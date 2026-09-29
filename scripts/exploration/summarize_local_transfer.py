"""Summarize selected completed transfer runs without upgrading their evidence."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import networkx as nx
from ase.io import read

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.reaction_network import atomic_json,aligned_rmsd
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT/'reports/local_transfer_v3')
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    args=p.parse_args()
    args.root=args.root.resolve()
    rows=[]
    for run in args.runs:
        run=run.resolve()
        manifest=json.loads((run/'manifest.json').read_text())
        for start in manifest['starts']:
            for strategy in manifest['strategies']:
                path=run/start/strategy/'network.json'
                n=json.loads(path.read_text(encoding='utf-8'))
                if n['status']=='running':raise ValueError(f'Unfinished run: {path}')
                graph=nx.Graph();graph.add_nodes_from(range(len(n['nodes'])))
                graph.add_edges_from(e['nodes'] for e in n['edges'])
                connected=nx.node_connected_component(graph,0) if n['nodes'] else set()
                reachable=[e for e in n['edges'] if set(e['nodes'])<=connected]
                pairs=sorted({tuple(sorted(n['nodes'][i]['graph_smiles'] for i in e['nodes']))
                    for e in reachable if e['kind']=='chemical'})
                rows.append(dict(start=start,strategy=strategy,status=n['status'],
                    evaluations=n['evaluations'],initialization_evaluations=n['initialization_evaluations'],
                    attempts=len(n['attempts']),root_chemical_pairs=pairs,
                    root_edges=len(reachable),all_edges=len(n['edges']),
                    continued_from_new_node=any(a['source_node']!=0 for a in n['attempts']),
                    continued_from_new_species=any(n['nodes'][a['source_node']]['graph_smiles']!=
                        n['nodes'][0]['graph_smiles'] for a in n['attempts']),
                    outcomes=dict(Counter(a['status'] for a in n['attempts'])),
                    transferred_attempts=sum(bool((a['proposal'] or {}).get('transferred')) for a in n['attempts']),
                    file=path.relative_to(ROOT).as_posix()))
    coverage=json.loads((args.root/'coverage.json').read_text())
    qc=[];events=[];feedback=[]
    for path in sorted((args.root/'qc').glob('*/verification.json')):
        result=json.loads(path.read_text())
        if result['status'] in ('started','refining_source_ts','integrating_irc'):
            raise ValueError('QC still in progress')
        qc.append(dict(event_id=result['event_id'],status=result['status'],
            method=result['method'],basis=result['basis'],environment=result['environment'],
            gradients=result.get('gradient_evaluations'),seconds=result.get('elapsed_seconds'),
            physical_event_verified=result.get('physical_event_verified',False),
            expected_endpoint_match=result.get('expected_endpoint_match'),
            endpoint_graphs=[e['graph_smiles'] for e in result.get('endpoints',[])],
            barriers_eV=[e['barrier_electronic_eV'] for e in result.get('endpoints',[])],
            file=path.relative_to(ROOT).as_posix()))
        feedback.append(dict(event_id=result['event_id'],status=result['status'],
            original_MLIP_pair_preserved=result.get('expected_endpoint_match'),
            actual_endpoint_graphs=qc[-1]['endpoint_graphs'],
            source_verification=qc[-1]['file'],model_weights_updated=False,
            interpretation='Evidence for this attempted connection, not reaction impossibility'))
        if result.get('physical_event_verified'):
            ts=read(path.parent/'ts.xyz')
            qc[-1]['TS_RMSD_A']=aligned_rmsd(np.asarray(result['source']['positions_A']['ts']),ts.positions)
            events.append(dict(event_id=result['event_id'],atomic_numbers=ts.numbers.tolist(),
                positions_A=dict(ts=ts.positions.tolist(),endpoints=[
                    read(path.parent/f'minimum_{side}.xyz').positions.tolist() for side in ('forward','reverse')]),
                charge=result['charge'],multiplicity=result['multiplicity'],method=result['method'],
                basis=result['basis'],environment=result['environment'],endpoint_graphs=qc[-1]['endpoint_graphs'],
                barriers_eV=qc[-1]['barriers_eV'],physical_event_verified=True,is_IRC=True,
                verified_pair=False,symbolic_status='proposal_only_not_independently_reviewed',
                symbolic_proposal=result['source']['original_symbolic_proposal'],
                source_verification=qc[-1]['file']))
    for name,records in [('DFT_events',events),('DFT_feedback',feedback)]:
        (args.root/(name+'.jsonl')).write_text(''.join(json.dumps(r)+'\n' for r in records),encoding='utf-8')
    result=dict(rows=rows,DFT=qc,total_evaluations=sum(r['evaluations'] for r in rows),
        total_attempts=sum(r['attempts'] for r in rows),
        evidence='MLIP mode-displacement descents; DFT evidence applies only to listed QC events',
        family_holdout=False,trained_model=False,
        coverage_supported=sum(bool(r['proposals']) for r in coverage['panel']),
        coverage_panel_size=len(coverage['panel']))
    all_runs=list(args.root.glob('search*/*/*/*/network.json'))
    result['evaluations_including_initial_diagnostics']=sum(
        json.loads(p.read_text(encoding='utf-8'))['evaluations'] for p in all_runs)
    comparison=args.root/'historical_coverage_comparison.json'
    if comparison.exists():
        result['historical_coverage_comparison']=json.loads(comparison.read_text())
    continuation=args.root/'continuation/acetone_enol_observed/arrows/network.json'
    if continuation.exists():
        extra=json.loads(continuation.read_text(encoding='utf-8'))
        if extra['status']=='running':raise ValueError('Continuation probe still running')
        result['posthoc_continuation']=dict(status=extra['status'],evaluations=extra['evaluations'],
            attempts=len(extra['attempts']),source_graph=extra['start']['provenance']['reactant_smiles'],
            edges=[dict(kind=e['kind'],source_connected=e['source_connected'],
                graphs=[extra['nodes'][i]['graph_smiles'] for i in e['nodes']]) for e in extra['edges']],
            outcomes=dict(Counter(a['status'] for a in extra['attempts'])),
            file=continuation.relative_to(ROOT).as_posix(),included_in_primary_comparison=False)
        result['evaluations_including_continuation']=result['evaluations_including_initial_diagnostics']+extra['evaluations']
    diagnosis=args.root/'saddle_mode_diagnostic.json'
    if diagnosis.exists():
        result['posthoc_saddle_modes']=json.loads(diagnosis.read_text())
        result['all_MLIP_evaluations']=result.get('evaluations_including_continuation',
            result['evaluations_including_initial_diagnostics'])+result['posthoc_saddle_modes']['evaluations']
    refinement=args.root/'continuation_refinement/result.json'
    if refinement.exists():
        extra=json.loads(refinement.read_text())
        result['posthoc_tighter_refinement']=dict(status=extra['status'],
            evaluations=extra['total_evaluations'],ts=extra.get('ts'),
            endpoint_graphs=[e['graph_smiles'] for e in extra.get('endpoints',[])],
            source_connection=json.loads((refinement.parent/'connection_audit.json').read_text()),
            included_in_primary_comparison=False)
        result['all_MLIP_evaluations']+=extra['total_evaluations']
    atomic_json(args.root/'summary.json',result)
    lines=['# 局部电子动作迁移：v3 工程验证','',
        '固定起点为丙醛、丙酮，各一个随机种子。所有对照在同一起点使用相同势、位移幅度与预算。',
        '本轮改变了局部提议、seed 角度/异步进度和初始方向，不能与历史 v2 直接比较归因。','',
        '## 符号覆盖','',
        f"固定面板 {result['coverage_panel_size']} 个分子中 {result['coverage_supported']} 个有合法提议。"
        '该计数只表示图级回放合法，不证明存在相应 TS。',
        f"记录审计：{coverage['library_audit']}。来源文件和全部动作见 coverage.json。",'',
        '## 实际搜索','',
        '| 起点 | 策略 | 尝试 | 势评估 | 初始化评估 | 根连通化学图对 | 从新物种继续 |',
        '|---|---|---:|---:|---:|---:|---|']
    for r in rows:
        lines.append(f"| {r['start']} | {r['strategy']} | {r['attempts']} | {r['evaluations']} | "
            f"{r['initialization_evaluations']} | {len(r['root_chemical_pairs'])} | {r['continued_from_new_species']} |")
    lines+=['','图对按运行统计，反向事件和构象不算新的独立化学体系。候选数不是实际机理准确率。','',
        '## 已观察连接','']
    for r in rows:
        lines.append(f"- {r['start']} / {r['strategy']}："+
            ('；'.join(' ↔ '.join(pair) for pair in r['root_chemical_pairs']) or '无根连通化学连接')+
            f"；结果分类 {r['outcomes']}。")
    lines+=['','## DFT','']
    for r in qc:
        lines.append(f"- {r['event_id']}：{r['status']}；物理连接认证={r['physical_event_verified']}；"
            f"原 MLIP 端点对保留={r['expected_endpoint_match']}；实际端点={r['endpoint_graphs']}。")
        lines.append(f"  方法 {r['method']}/{r['basis']}，{r['environment']}；两侧电子能垒 {r['barriers_eV']} eV。")
    if qc:
        lines.append('DFT 与 AIMNet2-rxn 训练参考层级不同；能垒差不能全部归因于模型误差。DFT 复核没有产生独立电子箭头真值。')
    if not qc:lines.append('本报告尚无新增 DFT 事件；全部新增网络连接为 MLIP 候选。')
    if 'posthoc_continuation' in result:
        extra=result['posthoc_continuation']
        lines+=['','## 独立续探诊断','',
            f"直接使用新发现的烯醇端点，继续 {extra['attempts']} 次箭头搜索，消耗 {extra['evaluations']} 次势评估。",
            f"结果：{extra['outcomes']}；实际连接：{extra['edges']}。",
            '该探针在主试验结束后安排，不加入主对照，不作为独立反应体系或长期网络覆盖证据。']
    if 'posthoc_saddle_modes' in result:
        lines+=['','### 续探负模诊断','',
            '前三个拒绝候选重新计算 Hessian，共 186 次势评估，不改变原始判定。',
            '前两例的主负模约 −2110 cm⁻¹，与初始 Cartesian 方向的绝对重合约 0.83；'
            '额外负模约 −66～−79 cm⁻¹，主要位移在另一侧甲基氢。',
            '这是反应方向已接近、但横向柔性运动未稳定的线索；仅频率和方向重合不能独立完成模式归属。'
            '下一步可检验保持反应方向的横向曲率优化，不能直接忽略额外负模。',
            f"含该诊断的 MLIP 总评估数为 {result['all_MLIP_evaluations']}。"]
    if 'posthoc_tighter_refinement' in result:
        r=result['posthoc_tighter_refinement']
        lines+=['','### 收紧阈值的单独验证','',
            f"选择第一个被拒绝的续探候选，将力阈值收紧为 0.005 eV/Å，重新 Dimer 与双侧下降，消耗 {r['evaluations']} 次势评估。",
            f"结果 {r['status']}；实际端点 {r['endpoint_graphs']}；新端点与续探源极小值的几何/能量匹配见 continuation_refinement/connection_audit.json。",
            '该结果支持这一例的失败与收敛精度有关。它是原可逆图对的后续恢复，不是新的独立化学反应，也没有独立 DFT 验证这一个后处理候选。',
            '主对照 48 次和初次续探 3 次的原始判定均不修改。']
    lines+=['','## 起点诊断与证据边界','',
        '- 首轮 search 中两个起点均未通过极小值检查，原始输出保留。',
        '- search_refined 收紧起点力阈值；丙醛通过，丙酮仍有两个负模。',
        '- search_curvature 对丙酮沿负模微扰后重新最小化，额外计算计入初始化预算。',
        '- 两体系的最终运行路径及源码快照分别保存在 summary.json 指向的目录中。',
        '- 只有一个随机种子、两个相关体系；不作统计显著性或未见反应家族泛化声明。',
        '- 局部模式仅支持当前非芳香活动中心和中性闭壳层 CHNO 域；有限提议上限会截断动作空间。',
        '- 严格初始电子源检查排除了部分原库的中继式表达；不等于判定这些反应不可能。',
        '- 完整箭头仍是提议条件，未成为独立审核的电子机理真值；未更新模型权重。','']
    lines += [f"主对照共 {result['total_attempts']} 次尝试、{result['total_evaluations']} 次势评估；"
              f"连同起点诊断共 {result['evaluations_including_initial_diagnostics']} 次势评估。",'',
              '真实坐标可视化：[分子结构、TS 和下降轨迹](search_curvature/acetone_s17/molecules/index.html)。','']
    (args.root/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
