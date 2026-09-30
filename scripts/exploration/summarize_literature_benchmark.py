"""Score literature targets only after physical search, requiring both endpoints."""
import json,sys
from collections import Counter
from pathlib import Path
from rdkit import RDLogger
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.symbolic_library import parse_explicit
from mechbridge.event_graph import geometry_mol,resonance_equivalent
from mechbridge.report_catalog import STRATEGIES
from mechbridge.reaction_network import atomic_json
from mechbridge.event_classification import matches_reference_pair


def main():
    RDLogger.DisableLog('rdApp.*');root=ROOT/'reports/literature_benchmark_v6'
    campaign=json.loads((root/'search/campaign.json').read_text())
    if 'finished_at' not in campaign:raise ValueError('Campaign still running')
    cases=json.loads((root/'cases.json').read_text());lookup={c['id']:c for c in cases['cases']}
    rows=[]
    for job in campaign['plan']:
        for strategy in campaign['config']['strategies']:
            path=root/'search'/job['folder']/job['start']/strategy/'network.json';n=json.loads(path.read_text())
            case=lookup[n['start']['provenance']['system']];reactant=parse_explicit(case['reactant'])
            targets=[parse_explicit(s) for s in case['targets']]
            molecules=[geometry_mol(n['start']['atomic_numbers'],v['positions_A'],0) for v in n['nodes']]
            hits=[];target_elsewhere=[]
            for edge in n['edges']:
                pair=[molecules[i] for i in edge['nodes']]
                if any(resonance_equivalent(m,t) for m in pair for t in targets):target_elsewhere.append(edge['id'])
                if matches_reference_pair(pair,reactant,targets):
                    hits.append(edge['id'])
            rows.append(dict(case=case['id'],start=n['start']['id'],strategy=strategy,status=n['status'],
                attempts=len(n['attempts']),evaluations=n['evaluations'],all_edges=len(n['edges']),
                outcomes=dict(Counter(t['status'] for t in n['attempts'])),
                target_hit=bool(hits),target_edges=hits,target_seen_in_other_pair=target_elsewhere,
                initial_graph_preserved=bool(molecules and resonance_equivalent(reactant,molecules[0])),
                file=path.relative_to(ROOT).as_posix()))
    aggregates=[]
    for case in cases['cases']:
        for strategy in campaign['config']['strategies']:
            subset=[r for r in rows if r['case']==case['id'] and r['strategy']==strategy]
            if subset:aggregates.append(dict(case=case['id'],strategy=strategy,strategy_label=STRATEGIES[strategy],
                runs=len(subset),hits=sum(r['target_hit'] for r in subset),attempts=sum(r['attempts'] for r in subset),
                evaluations=sum(r['evaluations'] for r in subset),
                outcomes=dict(sum((Counter(r['outcomes']) for r in subset),Counter()))))
    summary=dict(rows=rows,aggregates=aggregates,attempts=sum(r['attempts'] for r in rows),
        evaluations=sum(r['evaluations'] for r in rows),campaign_seconds=campaign['finished_at']-campaign['started_at'],
        reference_models_executed=False,evidence='MLIP index-one and actual two-sided descent; no new DFT',
        selection=cases,comparison='Matched geometry/net-edit/full-arrow controls; manual epoxide grammar was added before execution')
    atomic_json(root/'summary.json',summary)
    lines=['# 文献案例 benchmark v6','',
        '比较相同案例、相同模型与预算下的几何/净变键/完整箭头提议。FlowER 和 ReactionAtlas 原模型没有执行；本轮没有新增 DFT。',
        'FlowER 小分子开环是论文 Fig.1b 的表示示例，不是单例 DFT 认证；ReactionAtlas 两例有论文中相应反应的讨论。',
        '完整的 FlowER 十步示例含 Br/F/Cl，不能用当前势能；原测试集中符合严格大小/元素/组分筛选的一条离子对记录，在几何 Lewis 识别阶段受阻。失败与限制见 cases.json，未删原子简化。','',
        '| 案例 | 策略 | 命中 / 取向 | 尝试 | 几何评估 |','|---|---|---:|---:|---:|']
    for r in aggregates:lines.append(f"| {r['case']} | {r['strategy_label']} | {r['hits']}/{r['runs']} | {r['attempts']} | {r['evaluations']} |")
    lines+=['','命中要求同一条通过 MLIP 检查的边同时匹配文献反应物与目标产物。仅出现目标分子、另一端不同，不算命中。',
        '只有两个相遇取向，不能宣称统计泛化。完整箭头并未在本轮稳定优于净变键；环氧开环未通过体现当前 seed/搜索覆盖的不足。',
        '不同轮次使用不同相遇构型和预算；v5 的甲醛二聚通过不能改写成 v6 命中。',
        '',f"总尝试 {summary['attempts']}；总几何评估 {summary['evaluations']}；主搜索墙钟 {summary['campaign_seconds']/60:.1f} 分钟。",'',
        '统一入口：../index.html；文献与差距：../benchmark.html；实际分子视图：../events.html?stage=literature_benchmark_v6。']
    (root/'RESULTS_zh.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summary['aggregates'],indent=2))


if __name__=='__main__':main()
