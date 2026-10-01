"""Summarize targeted tetrose <-> 2GO closure and the global TransitionNet."""
from collections import Counter
import json
from pathlib import Path
import sys

from jinja2 import Environment,FileSystemLoader,select_autoescape

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.molecular_visuals import symbolic_graph_key
from mechbridge.report_layout import NAVIGATION,prepare_shared_assets

BASE=ROOT/'reports/formose_closure_v13'


def main():
    campaign=json.loads((BASE/'campaign.json').read_text(encoding='utf-8'))
    networks=[];events=[]
    for path in sorted((BASE/'runs').rglob('network.json')):
        network=json.loads(path.read_text(encoding='utf-8'));networks.append(network)
        for edge in network['edges']:
            attempt=network['attempts'][edge['attempt']]
            endpoints=[network['nodes'][i]['graph_smiles'] for i in edge['nodes']]
            source=network['nodes'][attempt['source_node']]['graph_smiles']
            predicted=attempt['proposal'].get('predicted_graph')
            events.append(dict(start=network['start']['id'],endpoints=endpoints,barriers_eV=edge['barriers_eV'],
                template_id=attempt['proposal'].get('template_id'),proposal_match=bool(predicted and
                    sorted(map(symbolic_graph_key,endpoints))==sorted(map(symbolic_graph_key,[source,predicted]))),
                file=(path.parent/f"attempt_{edge['attempt']:03d}/result.json").relative_to(ROOT).as_posix()))
    tetrose=symbolic_graph_key('O=CC(O)C(O)CO');two_go=symbolic_graph_key('O=CCO.O=CCO')
    exact=[event for event in events if sorted(map(symbolic_graph_key,event['endpoints']))==sorted([tetrose,two_go])]
    if not exact:raise ValueError('Targeted campaign did not produce an exact closure edge')
    closure=min(exact,key=lambda event:max(event['barriers_eV']))
    global_provenance=json.loads((BASE/'global_transitionnet/molecules/provenance.json').read_text(encoding='utf-8'))
    v12=json.loads((ROOT/'reports/formose_cycle_v12/summary.json').read_text(encoding='utf-8'))
    statuses=Counter(attempt['status'] for network in networks for attempt in network['attempts'])
    summary=dict(model='aimnet2-2025 member0',reference='B97-3c',cross_screening_used=False,
        starts=len(networks),attempts=sum(len(n['attempts']) for n in networks),
        evaluations=sum(n['evaluations'] for n in networks),edges=len(events),exact_closure_edges=len(exact),
        closure=closure,attempt_statuses=dict(statuses),constitutional_cycle_complete=True,
        constitutional_segments_passed=6,constitutional_segments_total=6,
        stereo_resolved_cycle_complete=False,
        stereo_gap='growth branch reaches D-erythrose; exact closure edge starts from D-threose',
        prior_segments=v12['segments'][:-1],global_network=global_provenance)
    (BASE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# Formose 最后缺口：中性 tetrose ↔ 2GO','',
      '本轮只使用 AIMNet2‑2025。没有运行 AIMNet2‑rxn 或 ensemble 交叉筛选。','',
      f"在 {summary['starts']} 个构象起点上完成 {summary['attempts']} 次定向任务和 {summary['evaluations']} 次势能/力评估，登记 {summary['edges']} 条边；其中 {summary['exact_closure_edges']} 条严格连接中性 tetrose 与 2GO。",'',
      '## 严格命中','',
      f"`{closure['endpoints'][0]}` ↔ `{closure['endpoints'][1]}`，双向势垒 {closure['barriers_eV'][0]:.3f} / {closure['barriers_eV'][1]:.3f} eV。符号提议端点命中：{'是' if closure['proposal_match'] else '否'}。",'',
      '命中只出现在 4 个 threose 构象中的一个；4 个 erythrose 构象和 3 个 2GO 逆向构象均没有严格命中。','',
      '## 循环结论','',
      '- **组成级：** 结合 v12 的前五段，现在 6/6 分段均有物理边，得到完整的 GO 自催化增殖候选。',
      '- **立体分辨级：** 尚未闭合。C4 生长段到 D‑erythrose，本轮闭环边来自 D‑threose，二者之间没有物理互变边。','',
      '## 查看','',
      '- [全部 C2/C3/C4 分子和 92 条边的 TransitionNet](global_transitionnet/molecules/index.html)',
      '- [命中闭环 TS 与真实双侧下降](runs/threose_uff_c3/molecules/index.html)','',
      '## 证据边界','',
      '- 当前是 AIMNet2‑2025 气相 MLIP 一阶鞍点与双侧下降，不是 DFT/IRC。',
      '- 跨库存总网不比较节点绝对能量；各事件页面保留各自双向相对势垒。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    prepare_shared_assets();env=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    html=env.get_template('formose_closure.html').render(title='Formose 定向闭环',root='../',
        section='experiments',navigation=NAVIGATION,summary=summary)
    (BASE/'index.html').write_text(html,encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=True,indent=2))


if __name__=='__main__':main()
