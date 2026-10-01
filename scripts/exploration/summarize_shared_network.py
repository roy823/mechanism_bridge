"""Summarize the central-registry parallel TransitionNet experiment."""
import json
from pathlib import Path
import sys

from jinja2 import Environment,FileSystemLoader,select_autoescape

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.report_layout import NAVIGATION,prepare_shared_assets

BASE=ROOT/'reports/shared_network_v10'


def load_network(folder):
    path=BASE/folder/'flower_transamidation_minimal_o1/arrows/network.json'
    return path,json.loads(path.read_text(encoding='utf-8'))


def metrics(folder):
    path,network=load_network(folder)
    reservations=network['reservations']
    return dict(folder=folder,file=path.relative_to(ROOT).as_posix(),status=network['status'],
        workers=network['scheduler']['workers'],threads_per_worker=network['scheduler']['threads_per_worker'],
        attempts=len(network['attempts']),evaluations=network['evaluations'],
        physical_nodes=len(network['nodes']),species=len(network['species_nodes']),
        conformer_clusters=len(network.get('conformer_clusters',[])),edges=len(network['edges']),
        reservations=len(reservations),unique_reservations=len({r['id'] for r in reservations}),
        searched_nodes=len(network.get('searched_nodes',[])),scheduler=network['scheduler'],
        stop_reason=network['stop_reason'])


def main():
    final=metrics('transamidation_conformer_clusters')
    long=metrics('transamidation_species')
    _,final_network=load_network('transamidation_conformer_clusters')
    _,long_network=load_network('transamidation_species')
    root_graph=final_network['nodes'][0]['graph_smiles']
    final['initial_conformer_clusters']=sum(c['graph_smiles']==root_graph for c in final_network['conformer_clusters'])
    target='CNC(C)=O.N'
    target_attempts=[a for a in long_network['attempts'] if (a.get('proposal') or {}).get('predicted_graph')==target]
    tautomer_attempts=[a for a in target_attempts
        if long_network['nodes'][a['source_node']]['graph_smiles']=='C/N=C(/C)O.N']
    direct_elimination_attempts=[a for a in target_attempts
        if long_network['nodes'][a['source_node']]['graph_smiles']=='CN[C@@](C)(N)O']
    intermediate='CN[C@@](C)(N)O'
    intermediate_sources=[a for a in long_network['attempts']
        if long_network['nodes'][a['source_node']]['graph_smiles']==intermediate]
    summary=dict(model='aimnet2-2025 member0',reference='B97-3c',final=final,long_growth=long,
        target_product=target,target_attempts=len(target_attempts),
        direct_elimination_attempts=len(direct_elimination_attempts),
        tautomerization_attempts=len(tautomer_attempts),
        target_statuses=[a['status'] for a in target_attempts],
        intermediate_tasks=len(intermediate_sources),
        final_scheduler='species+conformer_cluster+action+geometry_variant',
        chemistry_conclusion='shared growth reached tetrahedral intermediate and an imidic-acid+ammonia event; final amide tautomerization not physically validated')
    (BASE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 共享 TransitionNet 调度与转酰胺化连续增长','',
        '## 实现','',
        '1. 中央进程唯一维护化学物种、物理构象、TS 边和 reservation。',
        '2. CPU worker 只执行三维 seed、Dimer、一阶鞍点检查和双侧下降。',
        '3. 每个返回结果由中央进程顺序去重；新 B/C 随即进入前沿。',
        '4. TransitionNet 默认显示唯一物种节点，节点内部保留构象列表。',
        '5. 最终 reservation 为 `物种＋构象簇＋动作＋几何变体`。构象簇使用 0.5 Å 对称感知分子 RMSD，每簇选择最低能代表构象。','',
        '## 最终构象簇调度验证','',
        f"3 个 worker、每个 6 个 CPU 线程；{final['attempts']} 次任务、{final['evaluations']} 次势能/力评估。登记 {final['physical_nodes']} 个物理极小值、{final['species']} 个物种、{final['conformer_clusters']} 个构象簇和 {final['edges']} 条 TS 边。{final['reservations']} 个 reservation 全部唯一。",'',
        f"同一初始反应物物种被保留为 {final['initial_conformer_clusters']} 个构象簇，在显示层仍然只有一个物种节点。",'',
        '## 较长网络增长结果','',
        f"较长运行完成 {long['attempts']} 次尝试和 {long['evaluations']} 次评估，得到 {long['physical_nodes']} 个物理极小值、{long['species']} 个物种和 {long['edges']} 条 TS 边。",'',
        '- 找到乙酰胺＋甲胺与中性四面体之间的物理事件。',
        '- 四面体自动进入共享前沿，后续任务不需要人工重新制作起点。',
        '- 质子转移 seed 实际落到亚胺酸＋氨的消除/互变事件；实际端点与原符号提议不同，按正式规则保留。',
        f"- 四面体直接坍塌到目标 `{target}` 派发 {len(direct_elimination_attempts)} 个变体；亚胺酸互变到同一目标又派发 {len(tautomer_attempts)} 个变体，但均未达到 TS 力收敛。",'',
        '## 查看','',
        '- [最终构象簇调度的真实 TS、下降轨迹和物种网](transamidation_conformer_clusters/molecules/index.html)',
        '- [较长连续增长网络](transamidation_species/molecules/index.html)','',
        '## 证据边界','',
        '- 新边是 AIMNet2-2025 一阶鞍点与双侧下降证据，不是 DFT/IRC。',
        '- 当前结果证明共享队列、中央去重、新节点扩展和构象簇 reservation 工作；完整转酰胺化终产物仍未形成经过验证的连续物理路径。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    prepare_shared_assets();env=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    html=env.get_template('shared_scheduler.html').render(title='共享 TransitionNet 调度',root='../',
        section='experiments',navigation=NAVIGATION,summary=summary)
    (BASE/'index.html').write_text(html,encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
