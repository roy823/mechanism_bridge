"""Summarize the parallel shared C3H6O3 formose/triose TransitionNet."""
from collections import Counter
import json
from pathlib import Path
import sys

from jinja2 import Environment,FileSystemLoader,select_autoescape
import networkx as nx

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.molecular_visuals import symbolic_graph_key
from mechbridge.report_layout import NAVIGATION,prepare_shared_assets
from mechbridge.species_network import project_species_network

BASE=ROOT/'reports/formose_shared_v11'
NETWORK=BASE/'formaldehyde_enediol_o1/hybrid/network.json'


def main():
    network=json.loads(NETWORK.read_text(encoding='utf-8'))
    projected=project_species_network(network['nodes'],network['edges'],network['root_component_nodes'])
    graph=nx.Graph();graph.add_nodes_from(range(len(projected['nodes'])))
    graph.add_edges_from(edge['nodes'] for edge in projected['edges'] if edge['nodes'][0]!=edge['nodes'][1])
    root=projected['physical_to_species'][0]
    distances=nx.single_source_shortest_path_length(graph,root)
    target=max(distances,key=distances.get);path=nx.shortest_path(graph,root,target)
    steps=[]
    for left,right in zip(path,path[1:]):
        candidates=[edge for edge in projected['edges'] if set(edge['nodes'])=={left,right}]
        edge=min(candidates,key=lambda item:max(item['barriers_eV']))
        attempt=network['attempts'][edge['attempt']]
        steps.append(dict(edge=edge['id'],left=projected['nodes'][left]['graph_smiles'],
            right=projected['nodes'][right]['graph_smiles'],barriers_eV=edge['barriers_eV'],
            seed_strategy=attempt['proposal'].get('seed_strategy'),
            template_id=attempt['proposal'].get('template_id'),
            predicted_graph=attempt['proposal'].get('predicted_graph')))
    proposal_matches=0;seed_edges=Counter()
    for edge in network['edges']:
        attempt=network['attempts'][edge['attempt']]
        pair=[network['nodes'][i]['graph_smiles'] for i in edge['nodes']]
        source=network['nodes'][attempt['source_node']]['graph_smiles']
        predicted=attempt['proposal'].get('predicted_graph')
        if predicted and sorted(map(symbolic_graph_key,pair))==sorted(map(symbolic_graph_key,[source,predicted])):
            proposal_matches+=1
        seed_edges[attempt['proposal'].get('seed_strategy')]+=1
    statuses=Counter(attempt['status'] for attempt in network['attempts'])
    summary=dict(model='aimnet2-2025 member0',reference='B97-3c',start='C=O.O/C=C/O',
        scope='fixed C3H6O3 gas-phase formose/triose network; not enzymatic glycolysis',
        workers=network['scheduler']['workers'],threads_per_worker=network['scheduler']['threads_per_worker'],
        attempts=len(network['attempts']),evaluations=network['evaluations'],reservations=len(network['reservations']),
        unique_reservations=len({item['id'] for item in network['reservations']}),
        physical_nodes=len(network['nodes']),conformer_clusters=len(network['conformer_clusters']),
        species=len(projected['nodes']),edges=len(network['edges']),chemical_edges=sum(e['kind']=='chemical' for e in network['edges']),
        conformational_edges=sum(e['kind']=='conformational' for e in network['edges']),
        species_components=nx.number_connected_components(graph),root_reachable_species=len(distances),
        max_species_depth=max(distances.values()),cycles=len(nx.cycle_basis(graph)),
        species_self_loops=sum(e['nodes'][0]==e['nodes'][1] for e in projected['edges']),
        proposal_matches=proposal_matches,source_connected_edges=sum(e['source_connected'] for e in network['edges']),
        arrow_edges=seed_edges['arrows'],geometry_edges=seed_edges['geometry'],attempt_statuses=dict(statuses),
        deepest_path=steps,deepest_species=projected['nodes'][target]['graph_smiles'],
        scheduler=network['scheduler'],network_file=NETWORK.relative_to(ROOT).as_posix())
    (BASE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 并行共享 C3H6O3 formose / 三碳糖 TransitionNet','',
        '完整细胞糖酵解包含十个酶促步骤以及 ATP/ADP、NAD⁺/NADH、无机磷酸、Mg²⁺和带电磷酸糖，不适合作为当前中性 CHNO 气相原型的直接测试。本轮选择 ReactionAtlas 同样关注的 formose 化学空间，在固定 C3H6O3 原子库存上检验较长连续网络。','',
        f"使用 {summary['workers']} 个 CPU worker（每个 {summary['threads_per_worker']} 线程），完成 {summary['attempts']} 个唯一 reservation 和 {summary['evaluations']} 次势能/力评估。得到 {summary['physical_nodes']} 个物理极小值、{summary['conformer_clusters']} 个构象簇、{summary['species']} 个物种和 {summary['edges']} 条 TS 边。",'',
        '## 网络拓扑','',
        f"物种图有 {summary['species_components']} 个连通分量；起点所在分量覆盖 {summary['root_reachable_species']}/{summary['species']} 个物种，最远最短路径为 {summary['max_species_depth']} 步，并出现 {summary['cycles']} 个简单回路。",'',
        f"30 条边中 {summary['arrow_edges']} 条由完整箭头 seed 找到、{summary['geometry_edges']} 条由几何 seed 找到；{summary['proposal_matches']} 条边的实际端点命中发起符号提议。这个计数没有同预算纯几何对照，不能单独解释为因果增益。",'',
        '## 一条最深物种路径','',
        '| 步 | A | B | 双向势垒/eV | seed |','|---:|---|---|---:|---|']
    for index,step in enumerate(steps,1):
        lines.append(f"| {index} | `{step['left']}` | `{step['right']}` | {step['barriers_eV'][0]:.3f} / {step['barriers_eV'][1]:.3f} | `{step['template_id'] or 'geometry'}` |")
    lines += ['', '第一步为甲醛＋烯二醇与甘油醛之间的羟醛成键事件。后续包含水合/脱水、羰基–烯醇互变、裂解和重组。','',
        '## 查看','', '- [真实 TS、双侧下降、物种合并 TransitionNet](molecules/index.html)','',
        '## 解释边界','',
        '- 这是无磷酸、无酶的 formose/三碳糖气相子网，不是糖酵解十步生物通路。',
        '- “10 步路径”位于物种合并图；同一物种的不同构象簇之间不一定已有显式构象转换 TS。',
        '- 所有新边均为 AIMNet2-2025 证据，没有进行 DFT/IRC 复核。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    prepare_shared_assets();env=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    html=env.get_template('formose_shared.html').render(title='并行三碳糖 TransitionNet',root='../',
        section='experiments',navigation=NAVIGATION,summary=summary)
    (BASE/'index.html').write_text(html,encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
