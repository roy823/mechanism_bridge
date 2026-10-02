"""Aggregate v17/v18 sparse-root searches and audit nine hidden glycolysis targets."""
import copy
import hashlib
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.io import write
import networkx as nx
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles
from mechbridge.molecular_visuals import entry, load_events
from mechbridge.network_aggregation import aggregate_network_records
from mechbridge.report_layout import attach_checks, molecular_document, prepare_shared_assets
from mechbridge.species_network import layout_species, project_species_network
from mechbridge.symbolic_library import parse_explicit
from scripts.data.prepare_glycolysis_sparse_endpoints import SPECIES

BASE = ROOT / 'reports/glycolysis_sparse_discovery_final_v18'
SOURCES = [ROOT / 'reports/glycolysis_sparse_endpoint_v17/runs',
           ROOT / 'reports/glycolysis_gap_recovery_v18/runs',
           ROOT / 'reports/glycolysis_gap_recovery_v18/protonation']


def skeleton_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def state(*names):
    return skeleton_key(graph_smiles(parse_explicit('.'.join(SPECIES[name] for name in names))))


def inventory(*names):
    mol = parse_explicit('.'.join(SPECIES[name] for name in names))
    return tuple(sorted(atom.GetAtomicNum() for atom in mol.GetAtoms())), Chem.GetFormalCharge(mol)


TARGETS = [
    ('G6P ↔ F6P', inventory('g6p'), state('g6p'), state('f6p')),
    ('FBP ↔ GAP + DHAP', inventory('fbp'), state('fbp'), state('gap','dhap')),
    ('GAP ↔ DHAP', inventory('gap'), state('gap'), state('dhap')),
    ('3PG ↔ 2PG', inventory('pg3'), state('pg3'), state('pg2')),
    ('2PG ↔ PEP + H2O', inventory('pg2'), state('pg2'), state('pep','water')),
    ('PEP + H2O ↔ pyruvate + Pi', inventory('pep','water'), state('pep','water'), state('pyruvate','pi')),
    ('glucose + Pi ↔ G6P + H2O', inventory('glucose','pi'), state('glucose','pi'), state('g6p','water')),
    ('FBP + H2O ↔ F6P + Pi', inventory('fbp','water'), state('fbp','water'), state('f6p','pi')),
    ('1,3-BPG + H2O ↔ 3PG + Pi', inventory('bpg13','water'), state('bpg13','water'), state('pg3','pi')),
]


def species_projection(group):
    projection = project_species_network(group['nodes'], group['edges'])
    physical = {node['id']: node for node in group['nodes']}
    nodes = []
    baseline = min(node['energy_eV'] for node in group['nodes'])
    for item in projection['nodes']:
        source = physical[item['representative_node']]
        mol = geometry_mol(source['atomic_numbers'], source['positions_A'], group['charge'])
        nodes.append(dict(id=item['id'], **entry(mol, source['positions_A']),
            energy=item['energy_eV']-baseline, conformer_count=len(item['physical_nodes']),
            physical_nodes=item['physical_nodes']))
    edges = projection['edges']
    layout = layout_species(nodes, edges)
    for node in nodes:
        node['layout'] = layout[node['id']]
    return projection, nodes, edges


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    paths = sorted({path.resolve() for source in SOURCES for path in source.rglob('network.json')})
    records = []
    for path in paths:
        network = json.loads(path.read_text(encoding='utf-8'))
        if network['status'] == 'completed':
            records.append(dict(id=path.relative_to(ROOT).as_posix(), network=network))
    groups = aggregate_network_records(records, group_by='inventory')
    group_for_record = {}
    projections = {}
    networks = []
    group_by_inventory = {}
    for key, group in groups.items():
        inv = tuple(group['composition']), group['charge']
        group_by_inventory[inv] = key
        for record in records:
            if (record['id'], 0) in group['node_map']:
                group_for_record[record['id']] = key
        projection, nodes, edges = species_projection(group)
        projections[key] = projection
        roots = sorted({projection['physical_to_species'][node] for node in group['root_nodes']})
        network_id = 'sparse:' + group['system'] + f':q{group["charge"]}'
        networks.append(dict(id=network_id, start=group['system'], strategy='aggregate',
            title=f"{group['system']} · 电荷 {group['charge']:+d} · 稀疏端点聚合网",
            nodes=nodes, edges=[], root_nodes=roots, physical_node_count=len(group['nodes'])))
    network_by_id = {network['id']: network for network in networks}
    events_by_edge = {}
    event_lookup = {}
    for path in paths:
        network = json.loads(path.read_text(encoding='utf-8'))
        if network['status'] != 'completed' or not network['edges']:
            continue
        loaded, _ = load_events(path.parents[2], [path], layout_network=False)
        record_id = path.relative_to(ROOT).as_posix()
        key = group_for_record[record_id]
        group = groups[key]
        projection = projections[key]
        network_id = 'sparse:' + group['system'] + f':q{group["charge"]}'
        for original in loaded:
            event = copy.deepcopy(original)
            aggregate_edge = group['edge_map'][(record_id, original['edge_id'])]
            event['id'] = 'v18_' + hashlib.sha1(
                f'{record_id}:{original["edge_id"]}'.encode()).hexdigest()[:12]
            physical = [group['node_map'][(record_id, node)] for node in original['physical_node_ids']]
            event['node_ids'] = [projection['physical_to_species'][node] for node in physical]
            event['network_id'] = network_id
            event['aggregate_edge_id'] = aggregate_edge
            events_by_edge[(key, aggregate_edge)] = event
            event_lookup[(key, aggregate_edge)] = event['id']
    events = list(events_by_edge.values())
    for key, group in groups.items():
        target = network_by_id['sparse:' + group['system'] + f':q{group["charge"]}']
        for edge in projections[key]['edges']:
            target['edges'].append(dict(edge, event_id=event_lookup[(key, edge['id'])], label=edge['id']))

    output = BASE / 'molecules'
    output.mkdir(exist_ok=True)
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
    (output / 'index.html').write_text(molecular_document(
        dict(events=events, networks=networks, figures=[], report_url='../RESULTS_zh.md'), output), encoding='utf-8')

    graph_by_inventory = {}
    evidence_by_inventory = {}
    for key, group in groups.items():
        inv = tuple(group['composition']), group['charge']
        graph = nx.Graph()
        evidence = {}
        for node in group['nodes']:
            graph.add_node(skeleton_key(node['graph_smiles']))
        for edge in group['edges']:
            a = skeleton_key(group['nodes'][edge['nodes'][0]]['graph_smiles'])
            b = skeleton_key(group['nodes'][edge['nodes'][1]]['graph_smiles'])
            graph.add_edge(a,b)
            evidence.setdefault(tuple(sorted((a,b))), []).append(dict(
                barriers_eV=edge['barriers_eV'], origins=edge['origins']))
        graph_by_inventory[inv] = graph
        evidence_by_inventory[inv] = evidence
    checks = []
    for label, inv, left, right in TARGETS:
        graph = graph_by_inventory[inv]
        present = left in graph and right in graph
        path = nx.shortest_path(graph,left,right) if present and nx.has_path(graph,left,right) else []
        steps = [evidence_by_inventory[inv][tuple(sorted((a,b)))] for a,b in zip(path,path[1:])]
        checks.append(dict(label=label,passed=bool(path),direct=graph.has_edge(left,right),
            path=path,path_edges=max(0,len(path)-1),step_evidence=steps))
    summary = dict(model='aimnet2-2025 member0', chemical_roots=9, initial_variants=18,
        completed_networks=len(records), evaluations=sum(r['network'].get('evaluations',0) for r in records),
        raw_edges=sum(len(r['network']['edges']) for r in records),
        aggregate_inventory_networks=len(networks), aggregate_species_nodes=sum(len(n['nodes']) for n in networks),
        aggregate_edges=sum(len(n['edges']) for n in networks), passed=sum(c['passed'] for c in checks),
        targets=len(checks), checks=checks, html='molecules/index.html')
    (BASE / 'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines = ['# 稀疏端点糖酵解网络重建', '',
        f"- 初始化学 roots：{summary['chemical_roots']}；构象/取向 starts：{summary['initial_variants']}",
        f"- 完成局部网络：{summary['completed_networks']}",
        f"- AIMNet2-2025 势能/力评估：{summary['evaluations']}",
        f"- 原始 TS 边：{summary['raw_edges']}",
        f"- 聚合固定库存网络：{summary['aggregate_inventory_networks']}；物种/边：{summary['aggregate_species_nodes']}/{summary['aggregate_edges']}",
        f"- 隐藏目标覆盖：{summary['passed']}/{summary['targets']}", '', '## 九段审计', '']
    for check in checks:
        lines.append(f"- {'通过' if check['passed'] else '未通过'}：{check['label']}；"
                     f"direct={check['direct']}；path={check['path_edges']} edges")
    lines += ['', '## 证据边界', '',
        '- 搜索只读取单侧 root、符号动作和实际发现节点；隐藏目标仅用于最终审计。',
        '- 物理证据为 AIMNet2-2025 气相一阶鞍点与双侧下降，不是酶催化自由能或 DFT IRC。',
        '- ATP/ADP、NAD+/NADH、Mg2+ 和蛋白环境尚未加入；磷酸化步骤是骨架替代反应。']
    (BASE / 'RESULTS_zh.md').write_text('\n'.join(lines), encoding='utf-8')
    provenance = dict(source_hashes={path.relative_to(ROOT).as_posix():hashlib.sha256(path.read_bytes()).hexdigest()
                                     for path in paths}, summary=summary)
    (output / 'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=True,indent=2))


if __name__ == '__main__':
    main()
