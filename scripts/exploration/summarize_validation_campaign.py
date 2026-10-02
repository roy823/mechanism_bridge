"""Report all planned repeats, budget curves and independently calculated QC results."""
import argparse
from collections import Counter,defaultdict
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
from mechbridge.reaction_network import atomic_json,aligned_rmsd


def graph_pairs(network,through_attempt=None):
    graph=nx.Graph();graph.add_nodes_from(range(len(network['nodes'])))
    edges=[e for e in network['edges'] if through_attempt is None or e['attempt']<=through_attempt]
    graph.add_edges_from(e['nodes'] for e in edges)
    if 0 not in graph:return set(),set(),[]
    connected=nx.node_connected_component(graph,0)
    reachable=[e for e in edges if set(e['nodes'])<=connected]
    pairs={tuple(sorted(network['nodes'][i]['graph_smiles'] for i in e['nodes']))
           for e in reachable if e['kind']=='chemical'}
    return pairs,connected,reachable


def summarize(root):
    root=Path(root).resolve()
    campaign=json.loads((root/'search/campaign.json').read_text())
    if 'finished_at' not in campaign:raise RuntimeError('Search campaign is still running')
    rows=[];curves=[];raw_outcomes=Counter()
    for job in campaign['plan']:
        for strategy in campaign['config']['strategies']:
            path=root/'search'/job['folder']/job['start']/strategy/'network.json'
            if not path.exists():raise RuntimeError(f'Missing planned run: {path}')
            n=json.loads(path.read_text(encoding='utf-8'))
            if n['status'] in ('running','aborted_error'):raise RuntimeError(f'Incomplete run: {path}')
            pairs,connected,edges=graph_pairs(n)
            cumulative=n.get('initialization_evaluations',n['evaluations'])
            attempts=[];first=None
            for trial in n['attempts']:
                cumulative+=trial['evaluations']
                observed,_,_=graph_pairs(n,trial['id'])
                attempts.append((cumulative,len(observed)))
                if observed and first is None:first=cumulative
            for budget in (1000,2000,4000,6000):
                count=0
                for cost,value in attempts:
                    if cost<=budget:count=value
                curves.append(dict(start=job['start'],seed=job['seed'],cohort=job['cohort'],
                    strategy=strategy,budget=budget,chemical_graph_pairs=count))
            statuses=Counter(t['status'] for t in n['attempts']);raw_outcomes.update(statuses)
            row=dict(start=job['start'],seed=job['seed'],cohort=job['cohort'],strategy=strategy,
                status=n['status'],evaluations=n['evaluations'],attempts=len(n['attempts']),
                chemical_graph_pairs=len(pairs),root_connected_edges=len(edges),
                all_discovered_edges=len(n['edges']),first_chemical_connection_evaluations=first,
                new_species=len({n['nodes'][i]['graph_smiles'] for i in connected}-
                                ({n['nodes'][0]['graph_smiles']} if n['nodes'] else set())),
                continued_from_new_node=any(t['source_node']!=0 for t in n['attempts']),
                unsupported_nodes=n['unsupported_nodes'],outcomes=dict(statuses),
                observed_graph_pairs=sorted(pairs),file=path.relative_to(ROOT).as_posix())
            row['candidate_quality_flags']=[dict(pair=pair,
                fragment_counts=[len(s.split('.')) for s in pair],
                contains_isolated_proton=any('[H+]' in s.split('.') for s in pair)) for pair in sorted(pairs)]
            rows.append(row)
    grouped=defaultdict(list)
    for r in rows:grouped[(r['cohort'],r['strategy'])].append(r)
    aggregates=[]
    for (cohort,strategy),group in grouped.items():
        aggregates.append(dict(cohort=cohort,strategy=strategy,runs=len(group),
            runs_with_chemical_connection=sum(r['chemical_graph_pairs']>0 for r in group),
            chemical_graph_pairs_sum=sum(r['chemical_graph_pairs'] for r in group),
            evaluations=sum(r['evaluations'] for r in group),
            continuation_runs=sum(r['continued_from_new_node'] for r in group)))
    paired=[]
    for start in campaign['config']['core_starts']:
        for seed in campaign['config']['core_random_seeds']:
            group={r['strategy']:r for r in rows if r['start']==start and r['seed']==seed}
            paired.append(dict(start=start,seed=seed,
                arrows_minus_net_edits=group['arrows']['chemical_graph_pairs']-group['bond_edits']['chemical_graph_pairs'],
                net_edits_minus_center=group['bond_edits']['chemical_graph_pairs']-group['center_random']['chemical_graph_pairs'],
                center_minus_geometry=group['center_random']['chemical_graph_pairs']-group['geometry']['chemical_graph_pairs']))
    qc=[]
    paired_events=[]
    for folder in ('acetaldehyde','nitromethane_unexpected'):
        path=root/'qc'/folder/'verification.json'
        if not path.exists():qc.append(dict(case=folder,status='missing'));continue
        result=json.loads(path.read_text())
        row=dict(case=folder,status=result['status'],
            physical_event_verified=result.get('physical_event_verified',False),
            MLIP_endpoint_pair_preserved=result.get('expected_endpoint_match'),
            method=result['method'],basis=result['basis'],
            observed_endpoints=[e['graph_smiles'] for e in result.get('endpoints',[])],
            barriers_eV=[e['barrier_electronic_eV'] for e in result.get('endpoints',[])],
            endpoint_checks=[dict(graph=e['graph_smiles'],force_converged=e['force_converged'],
                imaginary_count=e['imaginary_count'],lowest_frequency_cm=e['frequencies_cm-1'][0])
                for e in result.get('endpoints',[])],
            gradients=result.get('gradient_evaluations'),seconds=result.get('elapsed_seconds'),
            error=result.get('error'),file=path.relative_to(ROOT).as_posix())
        qc.append(row)
        recovery=path.parent/'curvature_recovery_forward/result.json'
        if recovery.exists():
            diagnosis=json.loads(recovery.read_text())
            row['posthoc_endpoint_diagnosis']=dict(status=diagnosis['status'],
                file=recovery.relative_to(ROOT).as_posix(),primary_verdict_unchanged=True,
                branches=diagnosis.get('branches',[]))
        if result.get('physical_event_verified'):
            source=result['source']
            from ase.io import read
            ts=read(path.parent/'ts.xyz')
            row['TS_rigid_aligned_RMSD_A']=aligned_rmsd(np.asarray(source['positions_A']['ts']),ts.positions)
            endfiles=['minimum_forward.xyz','minimum_reverse.xyz']
            paired_events.append(dict(event_id=result['event_id'],
                atomic_numbers=ts.numbers.tolist(),charge=result['charge'],multiplicity=result['multiplicity'],
                method=result['method'],basis=result['basis'],environment=result['environment'],
                positions_A={'ts':ts.positions.tolist(),
                    'endpoints':[read(path.parent/f).positions.tolist() for f in endfiles]},
                endpoint_graphs=row['observed_endpoints'],barriers_eV=row['barriers_eV'],
                physical_event_verified=True,is_IRC=True,
                MLIP_endpoint_pair_preserved=row['MLIP_endpoint_pair_preserved'],
                symbolic_proposal=source['original_symbolic_proposal'],
                symbolic_status='proposal_only_not_independently_reviewed',verified_pair=False,
                source_verification=path.relative_to(ROOT).as_posix()))
    (root/'DFT_events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in paired_events),encoding='utf-8')
    feedback=[dict(case=r['case'],source_verification=r['file'],
        primary_DFT_status=r['status'],original_MLIP_pair_preserved=r.get('MLIP_endpoint_pair_preserved'),
        observed_DFT_endpoint_graphs=r.get('observed_endpoints',[]),
        endpoint_checks=r.get('endpoint_checks',[]),
        posthoc_endpoint_diagnosis=r.get('posthoc_endpoint_diagnosis'),
        interpretation='Candidate compatibility evidence; not a label that the proposed reaction is impossible',
        model_weights_updated=False) for r in qc if 'file' in r]
    (root/'DFT_feedback.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in feedback),encoding='utf-8')
    summary=dict(aggregates=aggregates,rows=rows,paired_core=paired,budget_curves=curves,DFT=qc,
        outcomes=dict(raw_outcomes),independent_core_systems=2,random_repeats_per_core_system=3,
        source_model='AIMNet2-rxn member0, includes RGD1 pretraining',
        evidence_boundary='Search counts are MLIP descents; only listed representative events have separate DFT/IRC evidence.')
    atomic_json(root/'summary.json',summary)
    lines=['# 验证 v2：实际结果','',
        '所有策略统一内部位移范数；按固定配置运行，不使用参考 TS 或产品三维构型。统计单位为起点/随机种子的一次网络探索；反向和构象不是独立化学体系。','',
        '## 汇总','',
        '| 队列 | 方法 | 运行数 | 有根连通候选的运行 | 非恒等图对之和 | 势评估总次数 | 从新节点继续的运行 |',
        '|---|---|---:|---:|---:|---:|---:|']
    for r in aggregates:lines.append(f"| {r['cohort']} | {r['strategy']} | {r['runs']} | {r['runs_with_chemical_connection']} | {r['chemical_graph_pairs_sum']} | {r['evaluations']} | {r['continuation_runs']} |")
    lines+=['','core = 乙醛、硝基甲烷各三个随机种子；extension = 乙酰胺、环丁酮各一次。表内图对跨运行求和，并非全局不同反应数。这里统计的是 MLIP 势面上的非恒等图连接候选，不等于可信化学反应成功率。','',
        '## 原始逐次结果','',
        '| 起点 | 随机种子 | 方法 | 化学图对 | 首次发现所需势评估 | 总势评估 | 尝试 |',
        '|---|---:|---|---:|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['start']} | {r['seed']} | {r['strategy']} | {r['chemical_graph_pairs']} | {r['first_chemical_connection_evaluations']} | {r['evaluations']} | {r['attempts']} |")
    lines+=['','## 核心体系内配对差值','']
    for field,label in [('center_minus_geometry','中心随机 − 几何'),('net_edits_minus_center','净变键 − 中心随机'),('arrows_minus_net_edits','完整箭头 − 净变键')]:
        values=[r[field] for r in paired]
        lines.append(f"- {label}：化学图对差值 {values}；正/零/负 = {sum(v>0 for v in values)}/{sum(v==0 for v in values)}/{sum(v<0 for v in values)}。")
    lines+=['','仅两个核心化学体系，不作总体显著性或未见反应泛化声明。不同随机种子主要检验搜索稳定性。','',
        '## 独立 DFT/IRC','',
        '| 代表事件 | 状态 | DFT 连接验证 | 保留 MLIP 端点对 | DFT 实际端点 | 两侧电子势垒 / eV |',
        '|---|---|---|---|---|---|']
    for r in qc:
        lines.append(f"| {r['case']} | {r['status']} | {r.get('physical_event_verified')} | {r.get('MLIP_endpoint_pair_preserved')} | {' ↔ '.join(r.get('observed_endpoints',[]))} | {', '.join(f'{v:.4f}' for v in r.get('barriers_eV',[]))} |")
    for r in qc:
        for endpoint in r.get('endpoint_checks',[]):
            if endpoint['imaginary_count']:
                lines.append(f"\n- {r['case']} 的 `{endpoint['graph']}` 端点仍有 {endpoint['imaginary_count']} 个显著虚频，最低 {endpoint['lowest_frequency_cm']:.2f} cm⁻¹，不能认证为极小值。")
        if 'posthoc_endpoint_diagnosis' in r:
            diag=r['posthoc_endpoint_diagnosis']
            lines.append(f"\n- 单独的端点负模微扰诊断：{diag['status']}；原始 DFT 判定不改写。")
            for b in diag['branches']:
                lines.append(f"  - 方向 {b['sign']:+d}：`{b['graph_smiles']}`，极小值认证={b['minimum_certified']}，最低频率 {b['frequencies_cm-1'][0]:.2f} cm⁻¹。")
    lines+=['','DFT 标准：气相闭壳层 ωB97X/6-31G(d)。只对预先选定的代表事件复核，不能把全部 MLIP 网络升级为 DFT 网络。没有执行全局电子态稳定性认证，也未生成新的独立箭头真值。',
        '物理反馈保存在 `DFT_feedback.jsonl`：记录候选端点是否保留、实际 DFT 端点和补充诊断。未更新模型权重；候选未保留也不表示符号反应不可能。','',
        '## 候选质量','',
        '不能只凭候选数量判断方法能否预测真实机理。以下多碎片/离子结构仍计入原始对照，不因观察到结果而事后删掉，但必须单独检查：','']
    for r in rows:
        for flag in r['candidate_quality_flags']:
            if flag['contains_isolated_proton']:
                lines.append(f"- {r['start']} / seed {r['seed']} / {r['strategy']}：`{' ↔ '.join(flag['pair'])}`；包含孤立质子，未做 DFT 复核，不能视作可靠的化学发现。")
    lines+=['','完整箭头相对净变键的额外候选，需结合上述 DFT 结果判断；若其对应图对在 DFT 下未保留，不能将该额外计数解释为真实机理预测收益。','',
        '## 限制','',
        '- 所有方法有相同预算上限，但实际消耗随收敛与分支覆盖变化。',
        '- center_random 接收中心原子 ID，不接收形成/断裂方向；它是比纯随机几何更强的控制。',
        '- 源—汇进度耦合是提议启发式，不能保证适合异步反应。',
        '- 符号库完整体系匹配限制了新物种的后续扩展。',
        '- 新增起点使用反应物 ETKDG/UFF 初始化；未宣称重现催化、溶液或实验条件。',
        '- DFT 与 AIMNet2-rxn 的训练参考层级不同，能垒差异不能全部归因于模型预测误差。',
        '- 失败包括预算截断、非一阶驻点、未证实极小值、同谷返回、图感知失败，详见 JSON。','']
    (root/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for ax,cohort in zip(axes,('core','extension')):
        for strategy in campaign['config']['strategies']:
            ys=[]
            for budget in (1000,2000,4000,6000):
                values=[r['chemical_graph_pairs'] for r in curves if r['cohort']==cohort and r['strategy']==strategy and r['budget']==budget]
                ys.append(float(np.mean(values)))
            ax.plot([1000,2000,4000,6000],ys,marker='o',label=strategy)
        ax.set(title=cohort+' (pilot, descriptive only)',xlabel='Energy/force evaluation budget',ylabel='Mean root-connected candidate graph pairs')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(root/'budget_comparison.png',dpi=180);fig.savefig(root/'budget_comparison.pdf');plt.close(fig)
    print(json.dumps(dict(aggregates=aggregates,paired=paired,DFT=qc),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    summarize(parser.parse_args().root)
