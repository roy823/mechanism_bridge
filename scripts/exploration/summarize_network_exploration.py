"""Summarize actual search artifacts, keeping conformers and chemical events apart."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.reaction_network import atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    a=p.parse_args()
    manifest=json.loads((a.run/'manifest.json').read_text())
    expected=[a.run/s/strategy/'network.json' for s in manifest['starts'] for strategy in manifest['strategies']]
    if any(not p.exists() or json.loads(p.read_text(encoding='utf-8'))['status']!='completed' for p in expected):
        raise RuntimeError('Wait for every declared run to complete before comparing strategies')
    rows=[]
    networks=[]
    for path in sorted(a.run.glob('*/*/network.json')):
        n=json.loads(path.read_text(encoding='utf-8'))
        connected=set(n.get('root_component_nodes',[]))
        active=[e for e in n['edges'] if set(e['nodes']) <= connected]
        pairs={tuple(sorted(n['nodes'][i]['graph_smiles'] for i in e['nodes'])) for e in active if e['kind']=='chemical'}
        new_species={n['nodes'][i]['graph_smiles'] for i in connected}
        if n['nodes']:
            new_species.discard(n['nodes'][0]['graph_smiles'])
        row=dict(start=n['start']['id'],strategy=n['strategy'],status=n['status'],
            attempts=len(n['attempts']),evaluations=n['evaluations'],
            validated_connections=len(n['edges']),root_connected_edges=len(active),
            chemical_graph_pairs=len(pairs),new_species=len(new_species),
            source_nodes_searched=sorted({t['source_node'] for t in n['attempts']}),
            nodes_considered=n.get('expanded_nodes',[]),
            continued_from_new_node=any(t['source_node']!=0 for t in n['attempts']),
            outcomes=dict(Counter(t['status'] for t in n['attempts'])))
        rows.append(row); networks.append((path,n))
    events=[]
    for path,n in networks:
        for edge in n['edges']:
            attempt=n['attempts'][edge['attempt']]
            observed=[n['nodes'][i]['graph_smiles'] for i in edge['nodes']]
            proposed=attempt['proposal'].get('predicted_graph')
            source=n['nodes'][attempt['source_node']]['graph_smiles']
            matched=(sorted([source,proposed])==sorted(observed)) if proposed else None
            events.append(dict(run=path.parent.relative_to(a.run).as_posix(),edge_id=edge['id'],
                atomic_numbers=n['start']['atomic_numbers'],charge=0,multiplicity=1,
                nodes=edge['nodes'],endpoint_graphs=observed,ts_positions_A=edge['ts_positions_A'],
                endpoint_positions_A=[n['nodes'][i]['positions_A'] for i in edge['nodes']],
                barriers_eV=edge['barriers_eV'],proposal=attempt['proposal'],
                symbolic_proposal_matches_observed_graph_pair=matched,
                evidence='AIMNet2-rxn_index1_two_sided_displacement_descent',
                root_connected=set(edge['nodes'])<=set(n.get('root_component_nodes',[])),
                is_IRC=False,DFT_verified=False,independently_reviewed_arrow_pair=False))
    (a.run/'discovered_events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in events),encoding='utf-8')
    aggregate=[]
    for strategy in sorted({r['strategy'] for r in rows}):
        selected=[r for r in rows if r['strategy']==strategy]
        aggregate.append(dict(strategy=strategy, runs=len(selected),
            **{key:sum(r[key] for r in selected) for key in ['attempts','evaluations','validated_connections','root_connected_edges','chemical_graph_pairs','new_species']},
            continuation_runs=sum(r['continued_from_new_node'] for r in selected)))
    result=dict(rows=rows,aggregate=aggregate,
        interpretation='Sum across runs is not a count of globally independent reactions; acetaldehyde/vinyl alcohol share a family.',
        evidence='One pretrained AIMNet2-rxn member; index-one Hessians and mode-displacement/BFGS minima; not DFT IRC.',
        claim='Pilot comparison only; source/sink coupling is a seed heuristic; no assertion of general superiority.')
    result['globally_distinct_root_connected_graph_pairs']=sorted({
        tuple(sorted(e['endpoint_graphs'])) for e in events if e['root_connected']})
    atomic_json(a.run/'summary.json',result)
    lines=['# 符号与几何引导的反应网络探索：实际结果','',
        '输入仅为 RGD1 起始构型与 SynEPD 可回放符号库；未输入参考 TS、产物三维坐标或参考能垒。',
        '物理后端：AIMNet2-rxn 单个预训练成员。连接证据为一阶鞍点检查与沿负模位移后的双侧 BFGS 下降；不是 DFT IRC。','',
        '| 起点 | 策略 | 尝试 | 势评估次数 | 根节点连通边 | 化学图对 | 新物种 | 从新节点继续 |',
        '|---|---|---:|---:|---:|---:|---:|---|']
    for r in rows:
        lines.append(f"| {r['start']} | {r['strategy']} | {r['attempts']} | {r['evaluations']} | {r['root_connected_edges']} | {r['chemical_graph_pairs']} | {r['new_species']} | {r['continued_from_new_node']} |")
    lines+=['','## 比较口径','',
        '- 各策略的总势评估上限、单次上限、求解器与收敛标准相同；实际消耗分别报告。',
        '- 边数包括不同 TS 与构象连接；化学图对与新物种另计。',
        '- 几何策略不读取箭头。净变键和箭头策略使用相同已发表符号提议；箭头策略额外使用源—汇进度耦合。',
        '- 符号规则完整体系匹配，覆盖不到的新物种停止该分支；这限制探索范围。',
        '- 乙醛与乙烯醇属于同一反应家族，不能把三次起点当作三个独立化学家族。',
        '- AIMNet2-rxn 训练包含 RGD1；本实验不证明未见反应泛化。',
        '- 搜索失败、同谷返回与找到意外有效事件分开记录。',
        '- 不同网络可能找到同一事件，汇总数量不代表全局去重后的反应数。',
        '- 普通下降建立探索级连接；更严格的机理结论需要同层级 IRC 或 DFT 复核。','']
    lines+=['## 实际发现了什么','',
        '根节点连通分量内的全局不同化学图对如下；同图对的不同原子映射/TS 不重复算化学反应：','']
    for pair in result['globally_distinct_root_connected_graph_pairs']:
        lines.append(f'- `{pair[0]}` ↔ `{pair[1]}`')
    mismatches=[e for e in events if e['root_connected'] and e['symbolic_proposal_matches_observed_graph_pair'] is False]
    if mismatches:
        lines+=['','存在与符号提议不同的根连通事件。这些属于 MLIP 上的意外连接候选，尚无 DFT 支持，不能称作新发现的真实化学机理。']
    lines+=['','三组在每个起点使用相同预算上限。geometry 是随机几何种子，结果不代表原 GPA 完整几何策略或其他强反应搜索基线。完整箭头的额外收益应与净变键组比较，不能仅与随机种子比较。','']
    (a.run/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    if not rows:
        return
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    labels=[r['strategy'] for r in aggregate]
    axes[0].bar(labels,[r['chemical_graph_pairs'] for r in aggregate])
    axes[0].set_title('Chemical graph pairs (sum across starts)')
    axes[1].bar(labels,[r['evaluations'] for r in aggregate])
    axes[1].set_title('Actual energy/force evaluations')
    fig.savefig(a.run/'comparison.png',dpi=160); plt.close(fig)
    for k,(path,n) in enumerate(networks):
        g=nx.MultiGraph()
        g.add_nodes_from(v['id'] for v in n['nodes'])
        for e in n['edges']: g.add_edge(*e['nodes'])
        fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
        pos=nx.spring_layout(g,seed=17)
        nx.draw_networkx(g,pos,ax=ax,labels={v['id']:f"{v['id']}: {v['graph_smiles']}" for v in n['nodes']},node_color='#b9def3',font_size=8)
        ax.set_title(f"{n['start']['id']} / {n['strategy']} (MLIP)");ax.axis('off')
        filename=f'network_{k}.svg'
        fig.savefig(a.run/filename);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
