"""Report larger-system search cost and physically connected multistep candidates."""
import argparse
from collections import Counter
import html
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.reaction_network import atomic_json
from mechbridge.network_metrics import growth_metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    a=p.parse_args();root=a.root.resolve()
    campaign=json.loads((root/'search/campaign.json').read_text())
    if 'finished_at' not in campaign:raise ValueError('Campaign is not finished')
    rows=[];selection=None
    for job in campaign['plan']:
        for strategy in campaign['config']['strategies']:
            run=root/'search'/job['folder']
            path=run/job['start']/strategy/'network.json'
            n=json.loads(path.read_text(encoding='utf-8'))
            if n['status']=='running':raise ValueError('Unfinished network')
            metric=growth_metrics(n)
            rows.append(dict(start=job['start'],strategy=strategy,status=n['status'],
                attempts=len(n['attempts']),evaluations=n['evaluations'],
                model_calls=n['model_calls'],model_seconds=n['model_seconds'],seconds=n['elapsed_seconds'],
                outcomes=dict(Counter(t['status'] for t in n['attempts'])),**metric,
                file=path.relative_to(ROOT).as_posix()))
            if selection is None and strategy in ('arrows','hybrid') and metric['chains']:
                chain=metric['chains'][0]
                selection=dict(run=run.relative_to(ROOT).as_posix(),start=job['start'],strategy=strategy,
                    chain=chain,edges_to_verify=chain['chemical_edges'][:2],
                    rule=campaign['config']['DFT_selection'])
    qc=[]
    for path in sorted((root/'qc').glob('*/verification.json')):
        q=json.loads(path.read_text())
        qc.append(dict(case=path.parent.name,status=q['status'],
            physical_event_verified=q.get('physical_event_verified',False),
            original_pair_preserved=q.get('expected_endpoint_match'),
            endpoint_graphs=[e['graph_smiles'] for e in q.get('endpoints',[])],
            file=path.relative_to(ROOT).as_posix()))
    summary=dict(rows=rows,DFT_selection=selection,DFT=qc,
        attempts=sum(r['attempts'] for r in rows),evaluations=sum(r['evaluations'] for r in rows),
        model_calls=sum(r['model_calls'] for r in rows),campaign_seconds=campaign['finished_at']-campaign['started_at'],
        multistep_runs=sum(bool(r['chains']) for r in rows),
        runs_expanding_new_species=sum(bool(r['expanded_new_species']) for r in rows),
        neural_model_weights_updated=False,
        evidence='All search chains use actual MLIP minimum/TS connectivity; only individual listed QC events have DFT evidence')
    summary['strategy_totals']={strategy:dict(
        chemical_pairs_sum=sum(len(r['chemical_pairs']) for r in rows if r['strategy']==strategy),
        multistep_runs=sum(bool(r['chains']) for r in rows if r['strategy']==strategy),
        evaluations=sum(r['evaluations'] for r in rows if r['strategy']==strategy))
        for strategy in campaign['config']['strategies']}
    supplementary={}
    for name,relative in [('hessian','hessian_validation.json'),('dimer_extrapolation','dimer_extrapolation/summary.json'),
                          ('convergence_contract','convergence_contract/summary.json'),('DFT_chain','DFT_chain_audit.json')]:
        path=root/relative
        if path.exists():supplementary[name]=json.loads(path.read_text())
    control_rows=[]
    for path in sorted((root/'geometry_budget_check').glob('*/*/geometry/network.json')):
        n=json.loads(path.read_text(encoding='utf-8'))
        old=next(r for r in rows if r['start']==n['start']['id'] and r['strategy']=='geometry')
        previous=json.loads((ROOT/old['file']).read_text(encoding='utf-8'))
        control_rows.append(dict(start=n['start']['id'],status=n['status'],evaluations=n['evaluations'],
            original_evaluations=old['evaluations'],metrics=growth_metrics(n),
            original_attempt_prefix_preserved=[(t['status'],t['evaluations']) for t in previous['attempts']]==
                [(t['status'],t['evaluations']) for t in n['attempts'][:len(previous['attempts'])]]))
    supplementary['geometry_budget_controls']=control_rows
    path=root/'convergence_network_check/glycolaldehyde/hybrid/network.json'
    if path.exists():
        n=json.loads(path.read_text(encoding='utf-8'))
        supplementary['convergence_network_check']=dict(status=n['status'],attempts=len(n['attempts']),
            evaluations=n['evaluations'],metrics=growth_metrics(n))
    summary['supplementary']=supplementary
    atomic_json(root/'summary.json',summary)
    atomic_json(root/'DFT_selection.json',selection)
    lines=['# 扩展体系与多步网络增长 v4','',
        '参照 ReactionAtlas 的双来源提议、验证后扩展和新极小值队列设计；这是净中性小体系的受限实验，不是完整 formose、离子或溶液动力学复现。',
        '参考：[ReactionAtlas 方法与补充实验](https://arxiv.org/html/2606.30778v1)。','',
        '## 实验设计','',
        '八个固定起点、几何/符号/混合三种策略，各一个随机种子。相同 16000 个几何评估预算，最多 24 次尝试；'
        '每个几何独立计费，批量调用不折算成一次评估。当前策略根据已观察网络分配动作，没有更新神经网络权重。','',
        '| 起点 | 策略 | 尝试 | 几何评估 | 后端调用 | 用时/s | 根连通物种 | 化学图对 | 最长已证路径/化学步 | 新物种续探 |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['start']} | {r['strategy']} | {r['attempts']} | {r['evaluations']} | {r['model_calls']} | "
            f"{r['seconds']:.1f} | {len(r['root_species'])} | {len(r['chemical_pairs'])} | "
            f"{r['longest_demonstrated_chemical_path']} | {len(r['expanded_new_species'])} |")
    lines+=['','物种包含立体信息；不同构象不重复算物种。两步增长要求实际极小值网络存在路径，且压缩连续构象后至少有三个不重复物种。'
            '不通过合并同 SMILES 的未连接构象制造通路，也不把 A→B→A 往返计为增长。','',
            '## 实际多步候选','']
    for r in rows:
        if r['chains']:
            c=r['chains'][0]
            lines.append(f"- {r['start']} / {r['strategy']}：{' → '.join(c['species'])}；实际节点 {c['nodes']}；边 {c['edges']}。")
    if not summary['multistep_runs']:lines.append('本轮没有达到上述严格定义的多步候选。')
    lines+=['','## 独立 DFT','']
    for q in qc:lines.append(f"- {q['case']}：{q['status']}；端点对保留={q['original_pair_preserved']}；{q['endpoint_graphs']}。")
    if not qc:lines.append('尚无本轮 DFT 复核结果。选择规则和入选边已保存为 DFT_selection.json。')
    if 'DFT_chain' in supplementary:
        verdict=supplementary['DFT_chain']
        lines.append(f"完整原链 DFT 认证：{verdict['original_MLIP_chain_preserved']}；实际两步 DFT 链：{verdict['actual_two_step_DFT_chain']}。")
        lines.append('DFT 未收敛属于该计算尝试未认证，不是反应不存在的证据；替代端点与原端点分别保存。')
    lines+=['','## 主实验后的独立检查','',
        '下述检查不改写主实验尝试、预算或成功计数；当前代码已包含真实力停止条件的修正。','']
    for r in control_rows:
        lines.append(f"- 几何预算补查 {r['start']}：{r['evaluations']} 个评估，根连通化学图对 {len(r['metrics']['chemical_pairs'])}；"
                     f"原尝试前缀保持={r['original_attempt_prefix_preserved']}。")
    if 'hessian' in supplementary:
        h=supplementary['hessian']
        speeds=[r['speedup'] for r in h['cases']]
        lines.append(f"- Hessian 三点检查通过={h['passed']}，虚频数量保持；批量单项计时比约 {min(speeds):.2f}–{max(speeds):.2f}。并发负载会影响计时。")
    if 'dimer_extrapolation' in supplementary:
        d=supplementary['dimer_extrapolation']
        for r in d.get('comparisons',[]):
            lines.append(f"- Dimer 力外推 {r['case']}：评估 {r['evaluations_off']}→{r['evaluations_on']}，状态 {r['status_off']}→{r['status_on']}，端点图对相同={r['same_endpoint_graphs']}。")
    if 'convergence_contract' in supplementary:
        for r in supplementary['convergence_contract']['rows']:
            lines.append(f"- 停止条件检查 {r['start']} / {r['attempt']}：{r['old_status']}→{r['new_status']}，实际端点 {r['endpoints']}；"
                         '同图端点不自动计为新化学连接。')
    if 'convergence_network_check' in supplementary:
        c=supplementary['convergence_network_check']
        lines.append(f"- 修正后端到端短运行：{c['attempts']} 次尝试、{c['evaluations']} 次评估、{c['metrics']['all_edges']} 条 MLIP 边。")
    lines+=['','## 效率与限制','',
        f"总尝试 {summary['attempts']}；总几何评估 {summary['evaluations']}；后端调用 {summary['model_calls']}；"
        f"整个并行活动耗时 {summary['campaign_seconds']:.1f} 秒。单个运行用时包含算法、图处理和写盘，不含进程导入。", 
        '- Hessian 分批求力，保持原有限差分步长和认证阈值。计时不能外推成整体或 DFT 加速倍数。',
        '- 使用预训练 AIMNet2-rxn member0 与 ASE Dimer；所有 TS/端点采用 0.005 eV/Å 力阈值和完整曲率检查。',
        '- 网络相关优先级是启发式；混合组每四次中安排一次几何提议，并允许符号无覆盖节点的几何探索。',
        '- 单随机种子不能支持统计显著性；新物种可能仍在预训练数据分布内。',
        '- 没有把电子能垒当成自由能或实验反应速率；没有把搜索次数当成反应概率。',
        '- 双向模型训练、独立箭头真值和反应物分子池的双分子增长仍属后续工作。','']
    lines+=['可视化：[总览](index.html)、[丙醛实际分子网络](search/propanal_s17/molecules/index.html)、'
            '[环丁酮实际分子网络](search/cyclobutanone_s17/molecules/index.html)。','']
    (root/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
