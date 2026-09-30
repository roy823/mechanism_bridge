"""Explicit catalogue of recorded experiments; protocols are never pooled as a leaderboard."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
STAGES=[
 dict(id='network_exploration_v1',short='v1',title='初始网络探索',question='符号提议能否帮助找到鞍点连接？',
      conclusion='历史可行性结果；位移幅度存在混杂，不能作为最终公平对照。'),
 dict(id='validation_v2',short='v2',title='信息量与随机对照',question='统一位移与预算后，符号信息是否仍有价值？',
      conclusion='包含反应中心、净变键、完整箭头和几何对照；小样本及非物理碎裂候选单独审计。'),
 dict(id='local_transfer_v3',short='v3',title='局部符号迁移',question='局部电子动作能否迁移到不同取代基？',
      conclusion='支持的固定分子面板从 7/22 增至 16/22，但筛选规则也发生变化；丙酮互变已做 DFT/IRC。'),
 dict(id='network_growth_v4',short='v4',title='多步网络与效率',question='新极小值是否能继续扩展成多步网络？',
      conclusion='6/24 组运行出现严格的两步 MLIP 路径；所选两步链未获得完整 DFT 认证。'),
 dict(id='bimolecular_v5',short='v5',title='双分子探索边界',question='符号与相遇取向能否引导跨分子成键？',
      conclusion='跨分子重原子事件为几何 0、符号 8、混合 7；水合与甲醛二聚有独立 DFT/IRC 证据。'),
 dict(id='literature_benchmark_v6',short='v6',title='文献案例与匹配对照',question='面对作者使用的案例，差距落在提议、几何还是物理验证？',
      conclusion='相同案例与相同预算的本地测试；没有执行 FlowER 或 ReactionAtlas 原模型，不能称为原模型排名。')]
STRATEGIES={'geometry':'纯几何','center_random':'仅反应中心','bond_edits':'净变键','arrows':'完整箭头','hybrid':'符号＋几何'}


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def catalogue():
    stages=[]
    for definition in STAGES:
        stage=dict(definition);folder=ROOT/'reports'/stage['id'];summary_path=folder/'summary.json'
        if not summary_path.exists():
            if stage['short']!='v6':raise FileNotFoundError(summary_path)
            stage.update(state='运行中',runs=[],attempts=0,evaluations=0,summary_file=None)
            stages.append(stage);continue
        summary=read_json(summary_path);runs=[]
        for row in summary['rows']:
            if stage['short']=='v1':path=folder/row['start']/row['strategy']/'network.json'
            else:path=ROOT/row['file']
            n=read_json(path)
            roots=set(n.get('root_component_nodes',[]))
            runs.append(dict(start=row['start'],strategy=n['strategy'],strategy_label=STRATEGIES[n['strategy']],
                status=n['status'],attempts=len(n['attempts']),evaluations=n['evaluations'],
                all_edges=len(n['edges']),chemical_edges=sum(e['kind']=='chemical' for e in n['edges']),
                root_edges=sum(set(e['nodes'])<=roots for e in n['edges']),
                target_hit=row.get('target_hit'),
                file=path.relative_to(ROOT).as_posix(),node_count=len(n['nodes'])))
        stage.update(state='已完成',runs=runs,attempts=sum(r['attempts'] for r in runs),
            evaluations=sum(r['evaluations'] for r in runs),summary_file=summary_path.relative_to(ROOT).as_posix())
        stages.append(stage)
    return stages


def quantum_checks():
    checks=[]
    for path in sorted((ROOT/'reports').rglob('verification.json')):
        q=read_json(path)
        if 'physical_event_verified' not in q:continue
        passed=q['physical_event_verified'];matched=q.get('expected_endpoint_match')
        label=('原图对通过' if matched else '替代图对有效') if passed else '未通过 / 未完成'
        checks.append(dict(id=path.parent.name,stage=path.relative_to(ROOT/'reports').parts[0],
            file=path.relative_to(ROOT).as_posix(),status=q['status'],passed=passed,matched=matched,label=label,
            method=q['method']+'/'+q['basis'],endpoints=[e['graph_smiles'] for e in q.get('endpoints',[])],
            parent=q.get('source',{}).get('source_json'),event_id=q['event_id']))
    return checks
