"""Summarize physical searches for reaction classes showcased by FlowER."""
from collections import Counter
import json
from pathlib import Path
import sys

from jinja2 import Environment, FileSystemLoader, select_autoescape

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))

from mechbridge.molecular_visuals import load_events, symbolic_graph_key
from mechbridge.report_layout import NAVIGATION, prepare_shared_assets

BASE=ROOT/'reports/flower_cases_aimnet2025'
CASE_DIRS={'flower_diels_alder_minimal':'diels_alder',
           'flower_prilezhaev_minimal':'prilezhaev',
           'flower_transamidation_minimal':'transamidation'}


def pair_has(event, left, right):
    observed=sorted(symbolic_graph_key(s) for s in (event['left']['smiles'],event['right']['smiles']))
    return observed==sorted(symbolic_graph_key(s) for s in (left,right))


def is_case_target(case_id,event):
    if case_id=='flower_diels_alder_minimal':
        return pair_has(event,'C=CC=C.C=C','C1=CCCCC1')
    if case_id=='flower_prilezhaev_minimal':
        return pair_has(event,'C=C.OOC=O','C1CO1.O=CO')
    if case_id=='flower_transamidation_minimal':
        return pair_has(event,'CC(=O)N.CN','CNC(C)(N)O')
    return False


def main():
    definitions={d['id']:d for d in json.loads(
        (ROOT/'data/processed/flower_case_definitions.json').read_text(encoding='utf-8'))}
    rows=[];event_rows=[];status_counts=Counter()
    for case_id,directory in CASE_DIRS.items():
        folder=BASE/directory
        events,_=load_events(folder.resolve())
        networks=[]
        for path in sorted(folder.glob('*/*/network.json')):
            network=json.loads(path.read_text(encoding='utf-8'));networks.append(network)
            status_counts.update(a['status'] for a in network['attempts'])
        targets=0
        for event in events:
            network_path=next(folder.glob(f"{event['start']}/*/network.json"))
            network=json.loads(network_path.read_text(encoding='utf-8'))
            edge=network['edges'][event['edge_id']]
            attempt=network['attempts'][edge['attempt']]
            target=is_case_target(case_id,event);targets+=target
            event_rows.append(dict(case=definitions[case_id]['reaction_class'],start=event['start'],
                edge=event['edge_id'],endpoints=[event['left']['smiles'],event['right']['smiles']],
                barriers_eV=[event['barrier_forward'],event['barrier_reverse']],
                template_id=attempt['proposal'].get('template_id'),predicted_graph=event['predicted_graph'],
                proposal_match=bool(event['matched_proposal']),case_target=target,
                result_file=event['source_result']))
        attempts=sum(len(n['attempts']) for n in networks)
        evaluations=sum(n['evaluations'] for n in networks)
        edges=sum(len(n['edges']) for n in networks)
        if case_id=='flower_diels_alder_minimal':
            conclusion=f'发现 {targets} 条反应物↔环己烯目标边，并保留开链与带电旁路。'
        elif case_id=='flower_prilezhaev_minimal':
            conclusion=f'两个取向发现 {targets} 条烯烃＋过酸↔环氧化物＋甲酸目标边。'
        else:
            full_key=symbolic_graph_key('CNC(C)=O.N')
            full_found=any(full_key in {symbolic_graph_key(s) for s in e['endpoints']} for e in event_rows
                           if e['case']==definitions[case_id]['reaction_class'])
            conclusion=('找到胺加成的中性四面体事件；' +
                        ('已' if full_found else '未') + '形成完整取代产物与三步网络。')
        definition=definitions[case_id]
        rows.append(dict(id=case_id,reaction_class=definition['reaction_class'],reactant=definition['reactant'],
            runs=len(networks),attempts=attempts,evaluations=evaluations,edges=edges,
            case_target_edges=targets,conclusion=conclusion,viewer=f'{directory}/molecules/index.html'))
    aggregate=json.loads((BASE/'aggregate_summary.json').read_text(encoding='utf-8'))
    summary=dict(model='aimnet2-2025 member0',reference='B97-3c',
        scope='minimal representatives of FlowER reaction classes; not exact licensed NameRxn records',
        completed_runs=sum(r['runs'] for r in rows),attempts=sum(r['attempts'] for r in rows),
        evaluations=sum(r['evaluations'] for r in rows),edges=len(event_rows),
        proposal_matched_edges=sum(e['proposal_match'] for e in event_rows),
        case_target_edges=sum(e['case_target'] for e in event_rows),
        attempt_statuses=dict(status_counts),cases=rows,events=event_rows,aggregate=aggregate)
    (BASE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# FlowER 论文反应类型的 GPA 式物理搜索','',
        '本轮选择 FlowER 论文明确展示的 Diels–Alder、Prilezhaev 环氧化和转酰胺化。使用的是三类反应的最小代表体系，不是论文受许可限制的逐条 NameRxn 测试分子。','',
        '流程为：完整电子箭头 → 净变键与活性原子 → 反应物相遇取向和三维 seed → AIMNet2-2025/Dimer → 一阶鞍点检查 → 双侧下降 → 登记实际 B–TS–C。搜索没有读取参考 TS 或产物几何。','',
        f"共完成 {summary['completed_runs']} 个取向运行、{summary['attempts']} 次尝试和 {summary['evaluations']} 次势能/力评估；登记 {summary['edges']} 条 TS 边，其中 {summary['proposal_matched_edges']} 条命中各自符号提议，{summary['case_target_edges']} 条属于三类案例的目标事件。",'',
        '## 结果','',
        '| 反应类型 | 代表反应 | 尝试/评估 | TS 边 | 结论 |','|---|---|---:|---:|---|']
    for row in rows:
        lines.append(f"| {row['reaction_class']} | `{row['reactant']}` | {row['attempts']} / {row['evaluations']} | {row['edges']} | {row['conclusion']} |")
    lines += ['', '## 目标物理事件','',
        '| 类型 | 实际 B ↔ C | 双向势垒/eV | 符号模板 |','|---|---|---:|---|']
    for event in event_rows:
        if event['case_target']:
            lines.append(f"| {event['case']} | `{event['endpoints'][0]}` ↔ `{event['endpoints'][1]}` | {event['barriers_eV'][0]:.3f} / {event['barriers_eV'][1]:.3f} | `{event['template_id']}` |")
    lines += ['', '## 查看真实构型和 TransitionNet','',
        '- [物种聚合总网](aggregate/molecules/index.html)',
        '- [Diels–Alder 全部 TS 与下降轨迹](diels_alder/molecules/index.html)',
        '- [Prilezhaev 全部 TS 与下降轨迹](prilezhaev/molecules/index.html)',
        '- [转酰胺化全部 TS 与下降轨迹](transamidation/molecules/index.html)','',
        '## 当前证据能说明什么','',
        '- 三种 FlowER 符号反应类型都至少有一个符号引导 seed 落到了化学相关的一阶鞍点事件。',
        '- Diels–Alder 与 Prilezhaev 的完整目标端点对已出现；说明多箭头协同动作可以转成有效三维 seed。',
        '- 转酰胺化只得到加成/互变相关局部事件，尚未自动形成完整的加成、质子转移和离去网络。它直接暴露了下一步应做的“新节点同步扩展”。',
        '- 所有能垒均为 AIMNet2-2025 气相电子能结果；双侧 BFGS 下降不是严格 IRC，也没有完成 DFT 复核。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    prepare_shared_assets()
    env=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    html=env.get_template('flower_cases.html').render(title='FlowER 反应类型物理搜索',root='../',
        section='experiments',navigation=NAVIGATION,summary=summary)
    (BASE/'index.html').write_text(html,encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('completed_runs','attempts','evaluations','edges',
        'proposal_matched_edges','case_target_edges')},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
