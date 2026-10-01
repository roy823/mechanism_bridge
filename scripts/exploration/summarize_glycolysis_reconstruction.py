"""Audit glycolysis reconstruction and build one all-inventory TransitionNet."""
import hashlib
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.io import write
import networkx as nx
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from mechbridge.event_graph import geometry_mol
from mechbridge.molecular_visuals import entry, load_events
from mechbridge.report_layout import attach_checks, molecular_document, prepare_shared_assets
from mechbridge.species_network import layout_species

BASE = ROOT / 'reports/glycolysis_reconstruction_v14_final'
TARGETED = ROOT / 'reports/glycolysis_targeted_refine_v14'
RELAXED = ROOT / 'reports/glycolysis_relaxed_recovery_v15'
RELAXED_RESUME = ROOT / 'reports/glycolysis_relaxed_recovery_v15_resume'
STARTS = ROOT / 'data/processed/glycolysis_starts.jsonl'
DEFINITIONS = ROOT / 'data/processed/glycolysis_definitions.json'

PAIRS = [
    ('G6P ↔ F6P', 'hmp_g6p', 'hmp_f6p', 'full_atom'),
    ('FBP ↔ GAP + DHAP', 'cleavage_fbp', 'cleavage_gap_dhap', 'full_atom'),
    ('GAP ↔ DHAP', 'triose_gap', 'triose_dhap', 'full_atom'),
    ('3PG ↔ 2PG', 'pg_shift_3pg', 'pg_shift_2pg', 'full_atom'),
    ('2PG ↔ PEP + H2O', 'dehydration_2pg', 'dehydration_pep_water', 'full_atom'),
    ('PEP + H2O ↔ pyruvate + H3PO4', 'pep_hydrolysis_pep_water',
     'pep_hydrolysis_pyruvate_pi', 'phosphate_surrogate'),
    ('glucose + H3PO4 ↔ G6P + H2O', 'glucose_phosphorylation_glucose_pi',
     'glucose_phosphorylation_g6p_water', 'phosphate_surrogate'),
    ('F6P + H3PO4 ↔ FBP + H2O', 'f6p_phosphorylation_f6p_pi',
     'f6p_phosphorylation_fbp_water', 'phosphate_surrogate'),
    ('1,3-BPG + H2O ↔ 3PG + H3PO4', 'bpg_hydrolysis_bpg_water',
     'bpg_hydrolysis_3pg_pi', 'phosphate_surrogate'),
]


def constitution(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Cannot parse {smiles}')
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def network_paths():
    return (sorted((BASE / 'runs').rglob('network.json')) + sorted(TARGETED.rglob('network.json')) +
            sorted(RELAXED.rglob('network.json')) + sorted(RELAXED_RESUME.rglob('network.json')))


def main():
    RDLogger.DisableLog('rdApp.*')
    definitions = {row['id']: row for row in json.loads(DEFINITIONS.read_text(encoding='utf-8'))}
    networks = []
    for path in network_paths():
        data = json.loads(path.read_text(encoding='utf-8'))
        if data['status'] == 'completed':
            networks.append((path, data))
    edge_metrics = dict(chemical=0, conformational=0, source_connected=0,
                        proposal_endpoint_matches=0, carbon_skeleton_accepted=0)
    combined = nx.Graph()
    for _, network in networks:
        for edge in network.get('edges', []):
            edge_metrics[edge.get('kind', 'chemical')] += 1
            edge_metrics['source_connected'] += int(edge.get('source_connected', False))
            edge_metrics['carbon_skeleton_accepted'] += int(
                edge.get('endpoint_acceptance') == 'validated_core_descents')
            attempt = network['attempts'][edge['attempt']]
            proposal = attempt.get('proposal') or {}
            predicted = proposal.get('predicted_graph')
            source = network['nodes'][attempt['source_node']]['graph_smiles']
            observed = sorted(constitution(network['nodes'][i]['graph_smiles']) for i in edge['nodes'])
            if predicted and observed == sorted([constitution(source), constitution(predicted)]):
                edge_metrics['proposal_endpoint_matches'] += 1
            combined.add_edge(*observed)
    segments = []
    for label, left_id, right_id, evidence in PAIRS:
        target = {constitution(definitions[left_id]['canonical_smiles']),
                  constitution(definitions[right_id]['canonical_smiles'])}
        verified, unresolved = [], []
        for path, network in networks:
            for edge in network.get('edges', []):
                observed = {constitution(network['nodes'][i]['graph_smiles']) for i in edge['nodes']}
                if observed == target:
                    verified.append(dict(source=network['start']['id'], edge=edge['id'],
                        barriers_eV=edge['barriers_eV'],
                        acceptance=edge.get('endpoint_acceptance','legacy_full_system'),
                        network=path.relative_to(ROOT).as_posix()))
            for attempt in network.get('attempts', []):
                result_path = path.parent / attempt['artifact']
                if not result_path.exists():
                    continue
                result = json.loads(result_path.read_text(encoding='utf-8'))
                if len(result.get('endpoints', [])) != 2:
                    continue
                observed = {constitution(endpoint['graph_smiles']) for endpoint in result['endpoints']}
                if observed == target and attempt['status'] != 'new_connection':
                    unresolved.append(dict(source=network['start']['id'], attempt=attempt['id'],
                        status=attempt['status'], barriers_eV=[x['barrier_eV'] for x in result['endpoints']],
                        endpoint_force_max_eV_A=[x['force_max_eV_A'] for x in result['endpoints']],
                        result=result_path.relative_to(ROOT).as_posix()))
        target_nodes=sorted(target)
        path_states=(nx.shortest_path(combined,*target_nodes)
                     if len(target_nodes)==2 and all(combined.has_node(x) for x in target_nodes)
                        and nx.has_path(combined,*target_nodes) else [])
        segments.append(dict(name=label, evidence=evidence, recovered=bool(path_states),
                             direct_recovered=bool(verified),path=path_states,
                             verified_edges=verified, unresolved_exact_candidates=unresolved))
    summary = dict(model='aimnet2-2025 member0', reference='B97-3c',
        physical_scope='neutral fully protonated phosphate microstates; enzyme-free gas-phase surrogate',
        ATP_ADP_present=False, NAD_NADH_present=False, DFT_verified=False,
        networks=len(networks), completed_networks=sum(n['status']=='completed' for _, n in networks),
        starts=len(json.loads('[' + ','.join(STARTS.read_text(encoding='utf-8').splitlines()) + ']')),
        attempts=sum(len(n.get('attempts', [])) for _, n in networks),
        evaluations=sum(n.get('evaluations', 0) for _, n in networks),
        edges=sum(len(n.get('edges', [])) for _, n in networks), edge_metrics=edge_metrics,
        evaluations_per_registered_edge=(sum(n.get('evaluations', 0) for _, n in networks) /
            max(1, sum(len(n.get('edges', [])) for _, n in networks))), segments=segments,
        recovered_segments=sum(row['recovered'] for row in segments), total_segments=len(segments))
    BASE.mkdir(parents=True, exist_ok=True)
    (BASE / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    lines = ['# 无酶糖酵解代谢物网络重建：第一轮', '',
        f"共运行 {summary['networks']} 个起点网络、{summary['attempts']} 次 seed→TS 尝试、"
        f"{summary['evaluations']} 次势能/力评估，注册 {summary['edges']} 条 AIMNet2-2025 TS 边；"
        f"平均每条注册边 {summary['evaluations_per_registered_edge']:.0f} 次评估。", '',
        f"其中 {edge_metrics['chemical']} 条化学边、{edge_metrics['conformational']} 条构象边；"
        f"{edge_metrics['source_connected']} 条包含发起极小值，"
        f"{edge_metrics['proposal_endpoint_matches']} 条严格命中符号提议端点，"
        f"{edge_metrics['carbon_skeleton_accepted']} 条仅按含碳骨架验收。", '',
        '## 已知片段审计', '', '| 片段 | 证据类型 | 直接边 | 多步可达 | 未完全收敛的精确端点候选 |',
        '|---|---|---:|---:|---:|']
    for row in segments:
        lines.append(f"| {row['name']} | {row['evidence']} | {'是' if row['direct_recovered'] else '否'} | "
                     f"{'是' if row['recovered'] else '否'} | "
                     f"{len(row['unresolved_exact_candidates'])} |")
    lines += ['', f"严格恢复 {summary['recovered_segments']}/{summary['total_segments']} 个参考片段。", '',
              '## 严格恢复的物理边', '']
    for row in segments:
        for edge in row['verified_edges']:
            barriers = ' / '.join(f"{value:.3f}" for value in edge['barriers_eV'])
            lines.append(f"- **{row['name']}**：双向候选能垒 {barriers} eV；"
                         f"端点验收 `{edge['acceptance']}`；来源 `{edge['network']}`。")
    lines += ['', '## 证据边界', '',
        '- 所有物理边均为 AIMNet2-2025 鞍点与双侧下降，不是 DFT/IRC。',
        '- 磷酸基团采用完全质子化中性微观状态；结果不代表生理 pH。',
        '- H3PO4/H2O 只作为守恒原子的磷酸储库替代，不代表 ATP/ADP 酶促能垒。',
        '- 未加入 NAD+/NADH，因此 GAPDH 氧化还原步骤不在本轮物理搜索中。']
    (BASE / 'RESULTS_zh.md').write_text('\n'.join(lines), encoding='utf-8')

    output = BASE / 'global_transitionnet/molecules'
    output.mkdir(parents=True, exist_ok=True)
    events = []
    for path, network in networks:
        if network['status'] == 'completed' and network.get('edges'):
            required=[path.parent/f"attempt_{edge['attempt']:03d}"/'descent_-1.traj'
                      for edge in network['edges']]
            if not all(item.exists() for item in required):
                continue
            loaded, _ = load_events(path.parents[2], [path], layout_network=False)
            relative=path.relative_to(ROOT).as_posix()
            tag=('v15r' if 'glycolysis_relaxed_recovery_v15_resume' in relative else
                 'v15' if 'glycolysis_relaxed_recovery_v15' in relative else
                 'v14t' if 'glycolysis_targeted_refine_v14' in relative else 'v14')
            for event in loaded:
                event['id']=f"{tag}_{event['id']}"
            events.extend(loaded)
    nodes, node_by_key = [], {}

    def register(smiles, endpoint):
        if smiles not in node_by_key:
            node_by_key[smiles] = len(nodes)
            record = dict(endpoint, id=len(nodes), smiles=smiles, energy=0.,
                          energy_label='跨库存 · 不比较绝对能量')
            nodes.append(record)
        return node_by_key[smiles]

    starts = [json.loads(line) for line in STARTS.read_text(encoding='utf-8').splitlines()]
    for start in starts:
        mol = geometry_mol(start['atomic_numbers'], start['positions_A'], start['charge'])
        endpoint = entry(mol, start['positions_A'])
        endpoint['name'] = start['provenance']['system']
        register(endpoint['smiles'], endpoint)
    root_node_ids = list(range(len(nodes)))
    edges, pair_counts = [], {}
    for event in events:
        pair = [register(event[side]['smiles'], event[side]) for side in ('left', 'right')]
        event['physical_node_ids'] = event.get('physical_node_ids', event['node_ids'])
        event['node_ids'] = pair
        event['network_id'] = 'glycolysis-global'
        event['aggregate_edge_id'] = event['id']
        key = tuple(sorted(pair))
        pair_counts[key] = pair_counts.get(key, 0) + 1
        edges.append(dict(id=event['id'], label=len(edges), event_id=event['id'], nodes=pair,
                          barriers=[event['barrier_forward'], event['barrier_reverse']]))
    used = {}
    for edge in edges:
        key = tuple(sorted(edge['nodes']))
        edge['parallel_index'] = used.get(key, 0)
        edge['parallel_count'] = pair_counts[key]
        used[key] = edge['parallel_index'] + 1
    layout = layout_species(nodes, edges)
    for node in nodes:
        node['layout'] = layout[node['id']]
    for event in events:
        trajectory = []
        for frame in event['frames']:
            atoms = Atoms(numbers=event['numbers'], positions=frame['positions'])
            atoms.info['relative_energy_eV'] = frame['relative_energy']
            atoms.info['branch'] = frame['branch']
            trajectory.append(atoms)
        write(output / f"{event['id']}_path.xyz", trajectory, format='extxyz')
        write(output / f"{event['id']}_TS.xyz", trajectory[event['ts_index']], write_results=False)
        event['downloads'] = dict(path=f"{event['id']}_path.xyz", ts=f"{event['id']}_TS.xyz")
    prepare_shared_assets()
    attach_checks(events, output)
    network = dict(id='glycolysis-global', start='all known glycolysis metabolite inventories',
        strategy='aggregate', title='无酶糖酵解全代谢物 TransitionNet', nodes=nodes, edges=edges,
        root_nodes=root_node_ids)
    (output / 'global_network.json').write_text(
        json.dumps(network, ensure_ascii=False, indent=2), encoding='utf-8')
    payload = dict(events=events, networks=[network], figures=[], report_url='../../RESULTS_zh.md')
    (output / 'index.html').write_text(molecular_document(payload, output), encoding='utf-8')
    provenance = dict(events=len(events), species_nodes=len(nodes), edges=len(edges),
        source_hashes={path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path, _ in networks})
    (output / 'provenance.json').write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=True, indent=2))


if __name__ == '__main__':
    main()
