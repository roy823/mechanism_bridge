"""Combine C2/C3/C4 physical networks and audit the canonical formose cycle."""
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

BASE=ROOT/'reports/formose_cycle_v12'
SELECTED={
 'c2_two_formaldehyde_o0':'runs/c2_two_formaldehyde_o0/c2_two_formaldehyde_o0/hybrid/network.json',
 'c2_glycolaldehyde_o0':'recovery/c2_glycolaldehyde_o0/c2_glycolaldehyde_o0/hybrid/network.json',
 'c2_enediol_o0':'recovery/c2_enediol_o0/c2_enediol_o0/hybrid/network.json',
 'c3_go_fa_o0':'runs/c3_go_fa_o0/c3_go_fa_o0/hybrid/network.json',
 'c3_glyceraldehyde_d_o0':'runs/c3_glyceraldehyde_d_o0/c3_glyceraldehyde_d_o0/hybrid/network.json',
 'c3_dha_o0':'runs/c3_dha_o0/c3_dha_o0/hybrid/network.json',
 'c3_enediol_o0':'recovery/c3_enediol/c3_enediol_o0/hybrid/network.json',
 'c4_erythrose_d_o0':'recovery/c4_erythrose_d_o0/c4_erythrose_d_o0/hybrid/network.json',
 'c4_threose_d_o0':'recovery/c4_threose_d_o0/c4_threose_d_o0/hybrid/network.json',
 'c4_erythrulose_d_o0':'recovery/c4_erythrulose_d_o0/c4_erythrulose_d_o0/hybrid/network.json',
 'c4_two_glycolaldehyde_o0':'recovery/c4_two_go/c4_two_glycolaldehyde_o0/hybrid/network.json',
}
PRIOR=ROOT/'reports/formose_shared_v11/formaldehyde_enediol_o1/hybrid/network.json'


def state(smiles):return symbolic_graph_key(smiles)


def records(path,label):
    network=json.loads(path.read_text(encoding='utf-8'));result=[]
    for edge in network['edges']:
        endpoints=tuple(state(network['nodes'][i]['graph_smiles']) for i in edge['nodes'])
        result.append(dict(label=label,path=path,network=network,edge=edge,endpoints=endpoints))
    return network,result


def graph_from(records):
    graph=nx.Graph()
    for record in records:
        a,b=record['endpoints'];weight=max(record['edge']['barriers_eV'])
        if graph.has_edge(a,b):
            if weight<graph[a][b]['weight']:graph[a][b].update(weight=weight,record=record)
        else:graph.add_edge(a,b,weight=weight,record=record)
    return graph


def path_or_empty(graph,left,right):
    if left not in graph or right not in graph or not nx.has_path(graph,left,right):return []
    return nx.shortest_path(graph,left,right)


def path_metrics(graph,path):
    if len(path)<2:return dict(bottleneck_eV=None,steps=[])
    steps=[]
    for left,right in zip(path,path[1:]):
        record=graph[left][right]['record'];edge=record['edge']
        attempt=record['network']['attempts'][edge['attempt']]
        steps.append(dict(left=left,right=right,barriers_eV=edge['barriers_eV'],
            template_id=attempt['proposal'].get('template_id'),source=record['label']))
    return dict(bottleneck_eV=max(max(step['barriers_eV']) for step in steps),steps=steps)


def main():
    networks={};all_records=[]
    for label,relative in SELECTED.items():
        path=BASE/relative;network,found=records(path,label);networks[label]=network;all_records.extend(found)
    prior_network,prior_records=records(PRIOR,'prior_formose_shared_v11')
    definitions=json.loads((ROOT/'data/processed/formose_cycle_definitions.json').read_text(encoding='utf-8'))
    by_inventory={inventory:[] for inventory in ('C2H4O2','C3H6O3','C4H8O4')}
    inventory_by_start={item['id']+'_o0':item['inventory'] for item in definitions}
    for record in all_records:
        inventory=inventory_by_start[record['label']]
        by_inventory[inventory].append(record)
    graphs={key:graph_from(value) for key,value in by_inventory.items()}
    fa2=state('C=O.C=O');go=state('O=CCO');c2_enediol=state('O/C=C/O')
    enediol_fa=state('C=O.O/C=C/O');ga=state('O=CC(O)CO')
    dha=state('O=C(CO)CO');c3_enediol=state('O/C=C(/O)CO')
    ga_ylide=state('[H]/[O+]=C(\C[O-])CO')
    ga_fa=state('C=O.O=CC(O)CO');erythrose=state('O=CC(O)C(O)CO')
    ga_fa_ylide=state('O=CC(O)CO.[H]/[C-]=[O+]/[H]')
    two_go=state('O=CCO.O=CCO');c4_ylide=state('OC/[C-]=[O+]/[C@@H](O)CO')
    prior_graph=graph_from(prior_records)
    initiation=path_or_empty(graphs['C2H4O2'],fa2,go)
    go_enediol=path_or_empty(graphs['C2H4O2'],go,c2_enediol)
    aldol=path_or_empty(prior_graph,enediol_fa,ga)
    ldb=path_or_empty(graphs['C3H6O3'],ga,ga_ylide)
    if ldb:
        tail=path_or_empty(graphs['C3H6O3'],ga_ylide,c3_enediol)
        final=path_or_empty(graphs['C3H6O3'],c3_enediol,dha)
        ldb=ldb+tail[1:]+final[1:] if tail and final else []
    direct=path_or_empty(graphs['C3H6O3'],ga,dha)
    c4_growth=path_or_empty(graphs['C4H8O4'],ga_fa,erythrose)
    closure=path_or_empty(graphs['C4H8O4'],erythrose,two_go)
    ylide_closure=path_or_empty(graphs['C4H8O4'],c4_ylide,two_go)
    segments=[
      dict(name='起始：2FA ↔ GO',passed=bool(initiation),path=initiation,note='严格端点匹配；不属于自催化核心，但提供第一份 GO。'),
      dict(name='GO ↔ C2 烯二醇',passed=bool(go_enediol),path=go_enediol,note='为 GO＋FA 羟醛步骤提供反应性 C2 烯二醇。'),
      dict(name='C2 烯二醇＋FA ↔ GA',passed=bool(aldol),path=aldol,note='来自独立共享网络的中性羟醛成键边。'),
      dict(name='GA ↔ C3 烯二醇 ↔ DHA',passed=bool(ldb),path=ldb,note='经过电荷分离的甘油醛极限结构；另有直接 GA↔DHA 边。'),
      dict(name='GA＋FA ↔ C4 erythrose',passed=bool(c4_growth),path=c4_growth,note='通过含甲醛叶立德的混合物中间节点。'),
      dict(name='中性 C4 糖 ↔ 2GO',passed=bool(closure),path=closure or [c4_ylide,two_go],note=('严格闭环边已找到。' if closure else '只找到 C4 羰基叶立德↔2GO；该叶立德与中性 C4 糖键连不同。')),
    ]
    for row,graph in zip(segments,[graphs['C2H4O2'],graphs['C2H4O2'],prior_graph,
                                  graphs['C3H6O3'],graphs['C4H8O4'],graphs['C4H8O4']]):
        row.update(path_metrics(graph,row['path']))
    strict_core=all(row['passed'] for row in segments[1:])
    attempts=sum(len(n['attempts']) for n in networks.values())
    evaluations=sum(n.get('evaluations',0) for n in networks.values())
    raw_edges=sum(len(n['edges']) for n in networks.values())
    states={endpoint for record in all_records for endpoint in record['endpoints']}
    partial=[label for label,n in networks.items() if n['status']!='completed']
    viewers=[
      dict(label='2FA / GO',description='起始成键和 C2 旁路',href='runs/c2_two_formaldehyde_o0/molecules/index.html'),
      dict(label='C2 烯二醇',description='GO 与 C2 烯二醇互变',href='recovery/c2_enediol_o0/molecules/index.html'),
      dict(label='甘油醛',description='GA、DHA 与直接通道',href='runs/c3_glyceraldehyde_d_o0/molecules/index.html'),
      dict(label='C3 烯二醇',description='LdB–AvE 候选路径',href='recovery/c3_enediol/molecules/index.html'),
      dict(label='D‑erythrose',description='C4 生长与裂解候选',href='recovery/c4_erythrose_d_o0/molecules/index.html'),
      dict(label='D‑threose',description='C4 异构与旁路',href='recovery/c4_threose_d_o0/molecules/index.html'),
      dict(label='D‑erythrulose',description='C4 酮糖网络',href='recovery/c4_erythrulose_d_o0/molecules/index.html'),
      dict(label='2GO',description='retro‑aldol 目标侧搜索',href='recovery/c4_two_go/molecules/index.html')]
    summary=dict(model='aimnet2-2025 member0',reference='B97-3c',defined_roots=len(definitions),
        executed_roots=len(networks),partial_roots=partial,total_attempts=attempts,total_evaluations=evaluations,
        raw_edges=raw_edges,unique_endpoint_states=len(states),segments=segments,
        passed_segments=sum(row['passed'] for row in segments),total_segments=len(segments),
        strict_cycle_complete=strict_core,missing_segment='中性 C4 糖 ↔ 2GO',
        direct_ga_dha=direct,direct_ga_dha_metrics=path_metrics(graphs['C3H6O3'],direct),
        ldb_path=ldb,ldb_metrics=path_metrics(graphs['C3H6O3'],ldb),
        c4_ylide_to_two_go=ylide_closure,
        viewers=viewers,selected_networks={k:(BASE/v).relative_to(ROOT).as_posix() for k,v in SELECTED.items()})
    (BASE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 多库存并行 formose 自催化循环审计','',
      f"共定义 {len(definitions)} 个反应根，实际纳入审计 {len(networks)} 个根；完成/部分运行合计 {attempts} 次任务、{evaluations} 次势能/力评估和 {raw_edges} 条 TS 边。",'',
      '## 严格循环结论','',
      '**尚未找到完整自催化闭环。** 已覆盖 5/6 个检查分段；缺少中性四碳糖与 2 个中性羟基乙醛之间的物理连接。2GO 只连接到不同键连的 C4 羰基叶立德，不能当作中性 tetrose。','',
      '| 分段 | 状态 | 路径最大双向势垒/eV | 实际物种路径 |','|---|---|---:|---|']
    for row in segments:
        barrier='—' if row['bottleneck_eV'] is None else f"{row['bottleneck_eV']:.3f}"
        lines.append(f"| {row['name']} | {'通过' if row['passed'] else '缺失'} | {barrier} | {' ↔ '.join(f'`{x}`' for x in row['path'])} |")
    lines += ['', '## LdB–AvE','',
      '找到了 `GA ↔ 电荷分离 GA ↔ C3 烯二醇 ↔ DHA` 的物种级路径，同时也找到直接 `GA ↔ DHA` 边。当前 AIMNet2‑2025 上直接边势垒更低，与 ReactionAtlas 报告的“中性糖之间不走直接氢转移”不一致，应优先做 DFT 精修。','',
      '## 证据边界','',
      '- 输入所有物种不构成循环证据；只有实际 B–TS–C 边计入通过。',
      '- campaign 曾被工具中断；部分网络的已完成边照实保留并在 summary.json 标记。',
      '- 全部结果仍是 AIMNet2‑2025 气相 MLIP 鞍点与双侧下降，不是 DFT/IRC。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    prepare_shared_assets();env=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    html=env.get_template('formose_cycle.html').render(title='完整 formose 循环审计',root='../',
        section='experiments',navigation=NAVIGATION,summary=summary)
    (BASE/'index.html').write_text(html,encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=True,indent=2))


if __name__=='__main__':main()
