"""Audit and render the sparse charged glycolysis skeleton campaign."""
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
from mechbridge.intermolecular_actions import electron_actions
from mechbridge.molecular_visuals import entry, load_events
from mechbridge.network_aggregation import aggregate_network_records
from mechbridge.report_layout import attach_checks, molecular_document, prepare_shared_assets
from mechbridge.species_network import layout_species, project_species_network
from mechbridge.symbolic_library import parse_explicit, replay
from scripts.data.prepare_glycolysis_charged_sparse import SPECIES

BASE = ROOT / 'reports/glycolysis_charged_sparse_v16'
STARTS = ROOT / 'data/processed/glycolysis_charged_sparse_starts.jsonl'


def state(*names):
    return graph_smiles(parse_explicit('.'.join(SPECIES[name] for name in names)))


def enediols(name):
    mol = parse_explicit(SPECIES[name])
    result = set()
    for action_name, arrows in electron_actions(mol):
        if action_name != 'alpha_hydroxy_carbonyl_to_enediol':
            continue
        product, _ = replay(mol, arrows)
        result.add(graph_smiles(product))
    return result


def skeleton_key(smiles):
    """Constitutional key for benchmark scoring; display keeps stereochemistry."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Cannot parse endpoint graph: {smiles}')
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def species_projection(group):
    projection = project_species_network(group['nodes'], group['edges'])
    physical = {node['id']: node for node in group['nodes']}
    nodes = []
    for item in projection['nodes']:
        source = physical[item['representative_node']]
        mol = geometry_mol(source['atomic_numbers'], source['positions_A'], group['charge'])
        node = dict(id=item['id'], **entry(mol, source['positions_A']),
            energy=item['energy_eV'], conformer_count=len(item['physical_nodes']),
            physical_nodes=item['physical_nodes'])
        nodes.append(node)
    edges = projection['edges']
    layout = layout_species(nodes, edges)
    for node in nodes:
        node['layout'] = layout[node['id']]
    return projection, nodes, edges


def main():
    paths = sorted(BASE.rglob('network.json'))
    records = []
    for path in paths:
        network = json.loads(path.read_text(encoding='utf-8'))
        if network['status'] == 'completed':
            records.append(dict(id=path.relative_to(BASE).as_posix(), network=network))
    groups = aggregate_network_records(records, group_by='inventory')
    group_for_record = {}
    projections = {}
    networks = []
    for key, group in groups.items():
        for record in records:
            if (record['id'], 0) in group['node_map']:
                group_for_record[record['id']] = key
        projection, nodes, edges = species_projection(group)
        projections[key] = projection
        root_species = sorted({projection['physical_to_species'][node]
                               for node in group['root_nodes']})
        networks.append(dict(id='charged:' + group['system'], start=group['system'],
            strategy='aggregate', title=f"{group['system']} · 电荷 {group['charge']:+d} · 稀疏 root 聚合网",
            nodes=nodes, edges=[], root_nodes=root_species,
            physical_node_count=len(group['nodes'])))
    network_by_id = {network['id']: network for network in networks}
    events_by_edge = {}
    event_lookup = {}
    for path in paths:
        network = json.loads(path.read_text(encoding='utf-8'))
        if network['status'] != 'completed' or not network['edges']:
            continue
        run = path.parents[2]
        loaded, _ = load_events(run, [path], layout_network=False)
        record_id = path.relative_to(BASE).as_posix()
        key = group_for_record[record_id]
        group = groups[key]
        projection = projections[key]
        for original in loaded:
            event = copy.deepcopy(original)
            event['id'] = 'v16_' + hashlib.sha1(record_id.encode()).hexdigest()[:8] + '_' + event['id']
            aggregate_physical = [group['node_map'][(record_id, node)]
                                  for node in original['physical_node_ids']]
            event['node_ids'] = [projection['physical_to_species'][node]
                                 for node in aggregate_physical]
            event['network_id'] = 'charged:' + group['system']
            aggregate_edge = group['edge_map'][(record_id, original['edge_id'])]
            event['aggregate_edge_id'] = aggregate_edge
            event_lookup[(key, aggregate_edge)] = event['id']
            events_by_edge[(key, aggregate_edge)] = event
    events = list(events_by_edge.values())
    for key, group in groups.items():
        projection = projections[key]
        target = network_by_id['charged:' + group['system']]
        for edge in projection['edges']:
            event_id = event_lookup.get((key, edge['id']))
            if event_id is None:
                origin = group['edges'][edge['id']]['origins'][0]
                event_id = event_lookup[(key, group['edge_map'][(origin['record'], origin['edge'])])]
            target['edges'].append(dict(edge, event_id=event_id, label=edge['id']))
    out = BASE / 'aggregate' / 'molecules'
    out.mkdir(parents=True, exist_ok=True)
    for event in events:
        trajectory = []
        for frame in event['frames']:
            atoms = Atoms(numbers=event['numbers'], positions=frame['positions'])
            atoms.info['relative_energy_eV'] = frame['relative_energy']
            atoms.info['branch'] = frame['branch']
            trajectory.append(atoms)
        write(out / f"{event['id']}_path.xyz", trajectory, format='extxyz')
        write(out / f"{event['id']}_TS.xyz", trajectory[event['ts_index']], write_results=False)
        event['downloads'] = dict(path=f"{event['id']}_path.xyz", ts=f"{event['id']}_TS.xyz")
    prepare_shared_assets()
    attach_checks(events, out)
    payload = dict(events=events, networks=networks, figures=[], report_url='../../RESULTS_zh.md')
    (out / 'index.html').write_text(molecular_document(payload, out), encoding='utf-8')

    graph_by_group = {}
    for key, group in groups.items():
        graph = nx.Graph()
        for node in group['nodes']:
            graph.add_node(skeleton_key(node['graph_smiles']))
        for edge in group['edges']:
            left = skeleton_key(group['nodes'][edge['nodes'][0]]['graph_smiles'])
            right = skeleton_key(group['nodes'][edge['nodes'][1]]['graph_smiles'])
            graph.add_edge(left, right)
        graph_by_group[(tuple(group['composition']), group['charge'])] = graph
    starts = [json.loads(line) for line in STARTS.read_text(encoding='utf-8').splitlines()]
    inventory = {}
    for start in starts:
        inventory[start['provenance']['system']] = (tuple(sorted(start['atomic_numbers'])), start['charge'])

    shared_hexose = sorted(enediols('g6p') & enediols('f6p'))
    shared_triose = sorted(enediols('gap') & enediols('dhap'))
    checks = []
    primary = [
        ('G6P <-> F6P skeleton connection', 'g6p', state('g6p'), state('f6p')),
        ('GAP <-> DHAP skeleton connection', 'gap', state('gap'), state('dhap')),
        ('F6P + H2PO4- <-> FBP + H2O', 'f6p_pi', state('f6p', 'pi'), state('fbp', 'water')),
    ]
    for label, inventory_name, left, right in primary:
        graph = graph_by_group[inventory[inventory_name]]
        left, right = skeleton_key(left), skeleton_key(right)
        present = left in graph and right in graph
        path = nx.shortest_path(graph, left, right) if present and nx.has_path(graph, left, right) else []
        checks.append(dict(label=label, passed=bool(path), direct=graph.has_edge(left, right),
            shortest_path=path, nodes_present=present))
    for label, inventory_name, route in [
        ('stable G6P-enediol-F6P endpoint sequence', 'g6p', [state('g6p'), *(shared_hexose[:1]), state('f6p')]),
        ('stable GAP-enediol-DHAP endpoint sequence', 'gap', [state('gap'), *(shared_triose[:1]), state('dhap')]),
    ]:
        graph = graph_by_group[inventory[inventory_name]]
        route = [skeleton_key(node) for node in route]
        present = all(node in graph for node in route)
        step_edges = [graph.has_edge(a, b) for a, b in zip(route, route[1:])] if present else []
        checks.append(dict(label=label, passed=bool(present and all(step_edges)), route=route,
            nodes_present=present, step_edges=step_edges, auxiliary_mechanistic_check=True))
    summary = dict(model='aimnet2-2025 member0', starts=len(starts), completed_networks=len(records),
        evaluations=sum(record['network'].get('evaluations', 0) for record in records),
        raw_physical_nodes=sum(len(record['network']['nodes']) for record in records),
        raw_edges=sum(len(record['network']['edges']) for record in records),
        aggregate_species_nodes=sum(len(network['nodes']) for network in networks),
        aggregate_species_edges=sum(len(network['edges']) for network in networks),
        proposal_matched_events=sum(event.get('matched_proposal') is True for event in events),
        checks=checks, html='aggregate/molecules/index.html',
        scope='charged gas-phase skeleton benchmark with one stoichiometric Pi or water where required')
    (BASE / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    recovery_networks = summary['completed_networks'] - summary['starts']
    lines = ['# 带电稀疏糖酵解骨架搜索', '',
        f"- 初始 roots：{summary['starts']}（6 个化学起点，各 2 个构象/取向）",
        f"- 完成网络：{summary['completed_networks']}（含 {recovery_networks} 个 recovery/continuation 网络）",
        f"- AIMNet2-2025 势能/力评估：{summary['evaluations']}",
        f"- 原始物理极小值/TS 边：{summary['raw_physical_nodes']}/{summary['raw_edges']}",
        f"- 聚合物种节点/边：{summary['aggregate_species_nodes']}/{summary['aggregate_species_edges']}",
        f"- 实际端点匹配符号提议的事件：{summary['proposal_matched_events']}/{len(events)}", '',
        '## 目标通路']
    for check in checks:
        evidence = (f"direct={check['direct']}，shortest_path={len(check['shortest_path']) - 1} edges"
                    if 'direct' in check else f"逐步边={check['step_edges']}")
        lines.append(f"- {'通过' if check['passed'] else '未通过'}：{check['label']}；{evidence}")
    lines += ['', '## 解释边界', '',
        '- 每个双分子 seed 只含一个化学计量 Pi 或一个水分子，没有溶剂浴。',
        '- TS 要求全体系收敛且仅有一个显著虚频；端点允许按碳骨架力阈值验收。',
        '- 这是 AIMNet2-2025 气相 MLIP 双侧下降证据，不是酶催化机理或 DFT IRC。',
        '- 聚合 HTML 为每条聚合边保留最多 61 个采样构型帧，并单独保存 TS XYZ。',
        '- 旧版 F6P 输入的磷酸位置错误，旧 G6P–F6P 阴性结果不再作为证据。']
    (BASE / 'RESULTS_zh.md').write_text('\n'.join(lines), encoding='utf-8')
    provenance = dict(source_hashes={path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                                     for path in paths}, summary=summary)
    (out / 'provenance.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
