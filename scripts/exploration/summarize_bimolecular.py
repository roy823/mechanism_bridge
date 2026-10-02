"""Audit actual intermolecular events and compare held-out reference products."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from rdkit import Chem,RDLogger
from ase.io import read
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol,resonance_equivalent,graph_smiles
from mechbridge.symbolic_library import parse_explicit
from mechbridge.reaction_network import atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT/'reports/bimolecular_v5')
    a=p.parse_args();root=a.root.resolve();RDLogger.DisableLog('rdApp.*')
    campaign=json.loads((root/'search/campaign.json').read_text())
    if 'finished_at' not in campaign:raise ValueError('Wait for all planned runs')
    held=json.loads((root/'coley_held_out_reference.json').read_text())
    references={'coley_'+r['rxn_id']:parse_explicit(r['rxn_smiles'].split('>>')[1]) for r in held['records']}
    rows=[];events=[];selected=None
    for job in campaign['plan']:
        for strategy in campaign['config']['strategies']:
            run=root/'search'/job['folder']
            file=run/job['start']/strategy/'network.json'
            if not file.exists():raise ValueError('Missing run: '+str(file))
            n=json.loads(file.read_text(encoding='utf-8'))
            if n['status'] in ('running','aborted_error'):raise ValueError('Unfinished run')
            start=n['start'];system=start['provenance']['system']
            initial=geometry_mol(start['atomic_numbers'],start['positions_A'],0)
            connected=set(n.get('root_component_nodes',[]))
            these=[]
            for edge in n['edges']:
                endpoints=[n['nodes'][i] for i in edge['nodes']]
                mols=[geometry_mol(start['atomic_numbers'],e['positions_A'],0) for e in endpoints]
                chemistry=edge['endpoint_chemistry']
                domain_flags=[]
                if any(len(Chem.GetMolFrags(m))>2 for m in mols):domain_flags.append('more_than_two_endpoint_fragments')
                if any(len(f)==1 for m in mols for f in Chem.GetMolFrags(m)):domain_flags.append('isolated_atom_or_ion')
                if any(abs(atom.GetFormalCharge())>1 for m in mols for atom in m.GetAtoms()):domain_flags.append('formal_charge_magnitude_above_one')
                input_match=[resonance_equivalent(initial,m) for m in mols]
                reference_match=None
                if system in references:
                    reference_match=any(resonance_equivalent(references[system],m) for m in mols)
                event=dict(start=start['id'],system=system,strategy=strategy,edge=edge['id'],
                    attempt=edge['attempt'],nodes=edge['nodes'],**chemistry,
                    endpoint_graphs=[e['graph_smiles'] for e in endpoints],
                    endpoint_domain_flags=domain_flags,
                    input_pair_endpoint_match=input_match,
                    root_connected=set(edge['nodes'])<=connected,
                    reference_product_match=reference_match,
                    electronic_barriers_eV=edge['barriers_eV'],
                    pair_key=sorted(e['graph_smiles'] for e in endpoints),
                    network_file=file.relative_to(ROOT).as_posix(),
                    result_file=(file.parent/n['attempts'][edge['attempt']]['artifact']).relative_to(ROOT).as_posix())
                these.append(event);events.append(event)
                if (selected is None and strategy in ('arrows','hybrid') and
                    chemistry['classification']=='intermolecular_heavy_atom_bond'):
                    selected=dict(run=run.relative_to(ROOT).as_posix(),start=start['id'],strategy=strategy,
                        edge=edge['id'],rule=campaign['config']['DFT_selection'])
            heavy=[e for e in these if e['classification']=='intermolecular_heavy_atom_bond']
            rows.append(dict(start=start['id'],system=system,strategy=strategy,status=n['status'],
                attempts=len(n['attempts']),evaluations=n['evaluations'],seconds=n['elapsed_seconds'],
                nodes=len(n['nodes']),edges=len(n['edges']),heavy_events=len(heavy),
                unique_heavy_pairs=len({tuple(e['pair_key']) for e in heavy}),
                input_pair_heavy_events=sum(any(e['input_pair_endpoint_match']) for e in heavy),
                root_connected_heavy_events=sum(e['root_connected'] for e in heavy),
                hydrogen_transfer_events=sum(e['classification']=='intermolecular_hydrogen_transfer' for e in these),
                spectator_rearrangements=sum(e['classification']=='intramolecular_with_spectator' for e in these),
                outcomes=dict(Counter(t['status'] for t in n['attempts'])),
                initialized=n['status']=='completed',
                file=file.relative_to(ROOT).as_posix()))
    qc=[];verified_events=[];feedback=[]
    for path in sorted((root/'qc').glob('*/verification.json')):
        q=json.loads(path.read_text())
        qc.append(dict(file=path.relative_to(ROOT).as_posix(),status=q['status'],
            physical_event_verified=q.get('physical_event_verified',False),
            expected_endpoint_match=q.get('expected_endpoint_match'),
            endpoints=[e['graph_smiles'] for e in q.get('endpoints',[])],
            seconds=q.get('elapsed_seconds'),
            barriers_eV=[e['barrier_electronic_eV'] for e in q.get('endpoints',[])],
            imaginary_counts=[q.get('ts',{}).get('imaginary_count')]+
                             [e['imaginary_count'] for e in q.get('endpoints',[])]))
        feedback.append(dict(source=path.relative_to(ROOT).as_posix(),event_id=q['event_id'],
            status=q['status'],physical_event_verified=q.get('physical_event_verified',False),
            expected_endpoint_match=q.get('expected_endpoint_match'),model_weights_updated=False))
        if q.get('physical_event_verified'):
            endpoints=[]
            for direction,e in zip(('forward','reverse'),q['endpoints']):
                atoms=read(path.parent/f'minimum_{direction}.xyz')
                if graph_smiles(geometry_mol(atoms.numbers,atoms.positions,0))!=e['graph_smiles']:
                    raise ValueError('DFT endpoint export mismatch')
                endpoints.append(dict(**e,positions_A=atoms.positions.tolist()))
            verified_events.append(dict(event_id=q['event_id'],source=path.relative_to(ROOT).as_posix(),
                atomic_numbers=q['source']['atomic_numbers'],charge=q['charge'],multiplicity=q['multiplicity'],
                method=q['method'],basis=q['basis'],environment=q['environment'],
                ts_positions_A=read(path.parent/'ts.xyz').positions.tolist(),ts=q['ts'],endpoints=endpoints,
                physical_event_verified=True,independent_arrow_gold=False,
                arrow_status='seed_hypothesis_only_not_independently_validated',
                seed_hypothesis=q['source']['original_symbolic_proposal']))
    for name,records in [('DFT_events',verified_events),('DFT_feedback',feedback)]:
        (root/(name+'.jsonl')).write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in records),encoding='utf-8')
    totals={}
    for strategy in campaign['config']['strategies']:
        subset=[r for r in rows if r['strategy']==strategy]
        totals[strategy]={k:sum(r[k] for r in subset) for k in
            ['attempts','evaluations','heavy_events','unique_heavy_pairs','input_pair_heavy_events',
             'root_connected_heavy_events','hydrogen_transfer_events','spectator_rearrangements']}
        totals[strategy]['runs_with_heavy_event']=sum(r['heavy_events']>0 for r in subset)
    summary=dict(rows=rows,events=events,strategy_totals=totals,DFT_selection=selected,DFT=qc,
        attempts=sum(r['attempts'] for r in rows),evaluations=sum(r['evaluations'] for r in rows),
        campaign_seconds=campaign['finished_at']-campaign['started_at'],
        failed_jobs=campaign['failed'],neural_model_weights_updated=False,
        primary_evidence='MLIP stationary points and two-sided mode-displacement descent; not IRC',
        comparison='Full symbolic package including rigid orientation; not an isolated arrow ablation')
    summary['unique_heavy_graph_pairs']=sorted({tuple(e['pair_key']) for e in events
        if e['classification']=='intermolecular_heavy_atom_bond'})
    supplementary={}
    for name in ['coley_reference_ts_diagnostic','coley_reference_ts_strict']:
        path=root/name/'summary.json'
        if path.exists():
            supplementary[name]=[{k:r[k] for k in ['rxn_id','reference_TS_used','outcome',
                'reference_pair_preserved','endpoint_graphs','total_evaluations']} for r in json.loads(path.read_text())]
    for name in ['encounter_minima_diagnostic','encounter_minima_bfgs']:
        path=root/name/'summary.json'
        if path.exists():
            supplementary[name]=[dict(system=r['system'],rmsd_A=r['rmsd_A'],
                energy_difference_eV=r['energy_difference_eV'],evaluations=r['evaluations'],
                both_minima_converged=all(v['force_converged'] and v['imaginary_count']==0 for v in r['refined']),
                same_refined_minimum=r['same_refined_minimum_by_registry_criteria']) for r in json.loads(path.read_text())]
    summary['supplementary']=supplementary
    atomic_json(root/'summary.json',summary);atomic_json(root/'DFT_selection.json',selected)
    lines=['# 双分子探索边界验证 v5','',
        '## 方法和例子来源','',
        'ReactionAtlas 提供甲醛、水、小糖及其相遇复合物的研究背景；采用其随机刚体取向、接近与验证后扩展的思路。'
        '本实验使用 AIMNet2-rxn 与 Dimer，没有运行 MoreRed/MD-ET，也没有复现溶液 formose 动力学。'
        '[原文](https://arxiv.org/html/2606.30778v1)。','',
        'Coley 组实际数据选取 rxn_id 3217、3216：按中性闭壳层 CHNO 的显式原子数、反应 ID 排序，'
        '取前两个不同反应物组合。参考产物仅在搜索完成后比较。'
        '[数据与代码](https://github.com/coleygroup/dipolar_cycloaddition_dataset)。',
        '之前讨论的 Jin/Coley 强化学习论文案例为自由基氧化与自由基串联环化，超出本轮模型适用范围。'
        '[MIT 论文](https://hdl.handle.net/1721.1/151666)。','',
        '六组反应物、两种相对取向、三个策略；每次运行最多 12 次尝试和 12000 个几何评估。'
        '每个尝试最多 2000 次评估，物理力收敛阈值 0.005 eV/Å，TS 恰一个低于 -30 cm⁻¹ 的虚频，'
        '两端零显著虚频。初始极小值失败也保留。','',
        '符号动作包括公开 SynEPD 局部迁移与明确标注的人为电子动作规则；后者不是模型预测或独立金标准。'
        '符号组包含额外的反应位点取向采样，因此本轮比较整个符号引导方案，不单独归因于完整箭头。','',
        '## 实际结果','',
        f"主实验共 {summary['attempts']} 次尝试、{summary['evaluations']} 个几何评估，墙钟 {summary['campaign_seconds']/60:.1f} 分钟。",'',
        '| 起点 | 策略 | 状态 | 尝试 | 跨分子重原子事件 | 原始反应物匹配 | 严格根连通 | 评估 |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['start']} | {r['strategy']} | {r['status']} | {r['attempts']} | {r['heavy_events']} | "
                     f"{r['input_pair_heavy_events']} | {r['root_connected_heavy_events']} | {r['evaluations']} |")
    lines+=['','“原始反应物匹配”允许已枚举的共振等价，但不把同 SMILES 的不同相遇构型直接合并为同一极小值。'
        '“严格根连通”要求真实极小值图上的通路；脱离根的发现单独保留。重原子成键、氢转移、旁观分子存在时的单分子重排分开统计。',
        '事件数可包含同一化学图对的多个 TS，不能当作不同反应家族数。含孤立原子/离子、超过两个片段或高形式电荷的候选在结构化结果中另有标记；通过 MLIP 数值检查不代表它们具有可靠化学意义。','',
        '## 跨策略汇总','',
        '| 策略 | 尝试 | 几何评估 | 跨分子重原子事件 | 成功运行 / 12 | 严格根连通事件 |',
        '|---|---:|---:|---:|---:|---:|']
    for strategy,t in totals.items():
        lines.append(f"| {strategy} | {t['attempts']} | {t['evaluations']} | {t['heavy_events']} | {t['runs_with_heavy_event']} | {t['root_connected_heavy_events']} |")
    lines+=['',f"跨所有取向和策略去重后有 {len(summary['unique_heavy_graph_pairs'])} 个重原子成键图对；不等于相互独立的反应家族。",'',
        '## 实际反应图对','']
    for left,right in summary['unique_heavy_graph_pairs']:
        lines.append(f'- `{left}` ↔ `{right}`')
    lines+=['','## 独立 DFT 检查','',
        '| 检查 | 状态 | 物理事件通过 | 原图对保留 | 端点 | 气相电子能垒/eV |',
        '|---|---|---|---|---|---|']
    for q in qc:
        lines.append(f"| {Path(q['file']).parent.name} | {q['status']} | {q['physical_event_verified']} | "
                     f"{q['expected_endpoint_match']} | {' / '.join(q['endpoints'])} | {', '.join(f'{v:.4f}' for v in q['barriers_eV'])} |")
    lines+=['','水合为预先选择的主验证；二聚为事后补充，严格优化另存。能垒按上表端点顺序，参照实际相遇极小值，不是无限分离的单体能量之和。',
        '已验证事件坐标、Hessian 频率、端点与提议来源导出到 `DFT_events.jsonl`；失败记录保留于 `DFT_feedback.jsonl`。箭头仍是提议假设，没有独立金标准标签。','',
        '## Coley 参考 TS 诊断（不计入主实验）','',
        '| 条件 | 记录 | 结果 | 参考图对保留 | 评估 |',
        '|---|---|---|---|---:|']
    for name in ['coley_reference_ts_diagnostic','coley_reference_ts_strict']:
        for r in supplementary.get(name,[]):
            lines.append(f"| {name} | {r['rxn_id']} | {r['outcome']} | {r['reference_pair_preserved']} | {r['total_evaluations']} |")
    lines+=['','两组都使用作者参考 TS 和虚频方向；strict 将真实力阈值从 0.005 收紧到 0.001 eV/Å。'
        '3217 的通过证明该局部 MLIP 通路可以被保留，不能算作我们的反应物起始搜索命中。3216 仍有端点优化问题。','',
        '## 相遇盆地诊断','',
        '首次 BFGSLineSearch/0.0005 的五组检查均未同时收敛，作为失败诊断保留。随后 BFGS/0.001 的结果：','',
        '| 体系 | 两端均为合格极小值 | RMSD/Å | 能量差/eV | 合并条件满足 |',
        '|---|---|---:|---:|---|']
    for r in supplementary.get('encounter_minima_bfgs',[]):
        lines.append(f"| {r['system']} | {r['both_minima_converged']} | {r['rmsd_A']:.3f} | {r['energy_difference_eV']:.5f} | {r['same_refined_minimum']} |")
    lines+=['','这些诊断没有修改主实验网络。下一步需要实际探索相遇复合物之间的构象转换或分离/再相遇过程，不能用较宽 RMSD 阈值代替物理连接。','',
        '## 结论边界','',
        '气相 MLIP 电子能垒不能直接与 Coley 数据表中的水相活化自由能比较。'
        '下降轨迹不是 IRC；只有明确列出的独立 DFT 检查可升级证据。'
        '没有神经模型更新、浓度动力学或产率预测，也没有声称穷举全部符号空间。','',
        '每次运行仍固定原子库存；本轮验证两个分子之间的探索，没有动态加入第三个分子，也不把断裂片段自动拼接成已验证网络边。']
    (root/'RESULTS_zh.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(attempts=summary['attempts'],totals=totals,DFT_selection=selected),indent=2))


if __name__=='__main__':main()
