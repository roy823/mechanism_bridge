"""Summarize actual AIMNet2-2025 B-TS-C events from the formal exploration algorithm."""
from collections import Counter
import hashlib,json,sys,zipfile
from pathlib import Path
from jinja2 import Environment,FileSystemLoader,select_autoescape
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.report_layout import NAVIGATION,prepare_shared_assets
from mechbridge.reaction_network import atomic_json

BASE=ROOT/'reports/aimnet2025_reaction_paths'
LABELS={'geometry':'纯几何 seed','arrows':'完整箭头引导三维 seed','hybrid':'符号 seed＋纯几何 seed'}
VIEWERS={'intramolecular':'intramolecular/molecules/index.html','targeted_bimolecular':'targeted_bimolecular/molecules/index.html',
         'epoxide_ammonia':'epoxide_ammonia/molecules/index.html','bimolecular':'bimolecular/molecules/index.html'}
VIEWERS['frontier_expansion']='frontier_expansion/molecules/index.html'


def constitution(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol,canonical=True)


def relevant(start,pair):
    pair=set(pair)
    glycer={'O=C[C@H](O)CO','O=C[C@@H](O)CO','O=CC(O)CO'}
    if start.startswith('formaldehyde_glycolaldehyde'):
        return 'C=O.O=CCO' in pair and bool(pair&glycer)
    if start.startswith('formaldehyde_enediol'):
        return 'C=O.O/C=C/O' in pair and bool(pair&glycer)
    if start.startswith('coley_3217'):
        return pair=={'C#[N+][N-]C.C=C','CN1CCC=N1'}
    if start.startswith('flower_epoxide_ammonia'):
        return 'C1CO1.N' in pair and bool(pair&{'NCCO','[NH3+]CC[O-]'})
    if start=='glycolaldehyde':
        return pair=={'O=CCO','C=O.C=O'} or ('O=CCO' in pair and bool(pair&{'O/C=C/O','O/C=C\\O','OC=CO'}))
    if start=='glyceraldehyde':
        return bool(pair&glycer) and bool(pair&{'O=C(CO)CO','OC=C(O)CO'})
    return False


def main():
    rows=[];events=[]
    for path in sorted(BASE.rglob('network.json')):
        data=json.loads(path.read_text(encoding='utf-8'));campaign=path.relative_to(BASE).parts[0]
        complete=data['status']=='completed';target_count=0
        for edge in data['edges']:
            attempt=data['attempts'][edge['attempt']];pair=[data['nodes'][i]['graph_smiles'] for i in edge['nodes']]
            system=data['start'].get('provenance',{}).get('aggregate_system',data['start']['id'])
            hit=relevant(system,pair);target_count+=hit
            predicted=attempt['proposal'].get('predicted_graph')
            events.append(dict(start=data['start']['id'],system=system,strategy=data['strategy'],strategy_label=LABELS[data['strategy']],
                edge=edge['id'],endpoints=pair,barriers_eV=edge['barriers_eV'],source_connected=edge['source_connected'],
                predicted_graph=predicted,proposal_endpoint_match=(predicted is not None and
                    constitution(predicted) in {constitution(s) for s in pair}),
                template_id=attempt['proposal'].get('template_id'),
                literature_relevant=hit,result_file=(path.parent/f"attempt_{edge['attempt']:03d}"/'result.json').relative_to(ROOT).as_posix()))
        rows.append(dict(start=data['start']['id'],strategy=data['strategy'],strategy_label=LABELS[data['strategy']],
            status=data['status'],attempts=len(data['attempts']),evaluations=data['evaluations'],edges=len(data['edges']),
            target_edges=target_count,viewer=VIEWERS.get(campaign) if complete else None,file=path.relative_to(ROOT).as_posix()))
    unique_pairs={(e['system'],tuple(sorted(e['endpoints']))) for e in events}
    unique_targets={(e['system'],tuple(sorted(e['endpoints']))) for e in events if e['literature_relevant']}
    aggregate_path=BASE/'aggregate_summary.json';aggregate=json.loads(aggregate_path.read_text(encoding='utf-8')) if aggregate_path.exists() else None
    summary=dict(model='aimnet2-2025 member0',reference='B97-3c',algorithm='formal v1-v6 event registration',
        registration='A proposes a seed; the observed B-TS-C is registered even when A is neither endpoint',
        completed_runs=sum(r['status']=='completed' for r in rows),attempts=sum(r['attempts'] for r in rows),
        evaluations=sum(r['evaluations'] for r in rows),edges=len(events),unique_graph_pairs=len(unique_pairs),
        literature_relevant_edges=sum(e['literature_relevant'] for e in events),unique_literature_pairs=len(unique_targets),
        runs=rows,events=events,outcomes=dict(Counter(e['system'] for e in events)),aggregate=aggregate)
    atomic_json(BASE/'summary.json',summary)
    lines=['# AIMNet2-2025 正式算法反应路径扩展','',
        '`arrows` 表示完整箭头引导的三维 seed，不是只输出箭头。seed 使用净变键、活性原子、相遇取向、反应进度、进攻角以及电子源/受体耦合；随后由 AIMNet2-2025、Dimer 和双侧下降获得实际 B/C。','',
        f"完成运行 {summary['completed_runs']} 组，{summary['attempts']} 次尝试，{summary['evaluations']} 个几何评估，登记 {summary['edges']} 条 B–TS–C 边、{summary['unique_graph_pairs']} 个按体系区分的端点图对；其中 {summary['literature_relevant_edges']} 条边、{summary['unique_literature_pairs']} 个图对与选定文献案例相符。",'',
        '## 查看保存的真实反应过程','',
        '- [羟基乙醛与甘油醛](intramolecular/molecules/index.html)',
        '- [甲醛＋烯二醇与 Coley 环加成](targeted_bimolecular/molecules/index.html)',
        '- [环氧乙烷＋氨](epoxide_ammonia/molecules/index.html)',
        '- [甲醛＋羟基乙醛](bimolecular/molecules/index.html)','',
        '- [聚合同体系总 TransitionNet](aggregate/molecules/index.html)',
        '- [新盆地前沿扩展轨迹](frontier_expansion/molecules/index.html)','',
        '每个页面均使用保存的 TS 和双侧 BFGS 下降帧，没有插值。','',
        '## 文献相关命中','',
        '| 体系 | 实际端点 | 势垒/eV | 符号模板 | seed 预测端点命中 |','|---|---|---|---|---|']
    for event in events:
        if event['literature_relevant']:
            lines.append(f"| {event['start']} | `{event['endpoints'][0]}` ↔ `{event['endpoints'][1]}` | {event['barriers_eV'][0]:.3f} / {event['barriers_eV'][1]:.3f} | `{event['template_id']}` | {'是' if event['proposal_endpoint_match'] else '否'} |")
    if aggregate:
        lines += ['', '## 聚合同体系 TransitionNet','',
            f"逐运行共有 {aggregate['raw_nodes']} 个节点、{aggregate['raw_edges']} 条边；按能量与置换对齐 RMSD 合并后为 {aggregate['merged_nodes']} 个物理极小值、{aggregate['unique_ts_edges']} 条不同 TS。",'',
            '| 体系 | 来源运行 | 原节点→合并节点 | 原边→不同 TS | 连通分量 |','|---|---:|---:|---:|---:|']
        for row in aggregate['systems']:
            lines.append(f"| {row['system']} | {row['source_runs']} | {row['raw_nodes']}→{row['merged_nodes']} | {row['raw_edges']}→{row['unique_ts_edges']} | {row['components']} |")
    lines += ['', '## 证据边界','',
        '- 新边是 AIMNet2-2025 势能面上的一阶鞍点和双侧极小值，不是 DFT/IRC。',
        '- `source_connected` 仅说明实际端点是否包含发起节点 A；不作为 B–TS–C 事件的拒绝条件。',
        '- 相同图对允许多条不同 TS；高能分解和替代通道继续保留。']
    lines += ['', '## 案例来源','',
        '- 甲醛、羟基乙醛、烯二醇和甘油醛来自 ReactionAtlas 讨论的 formose 化学空间。',
        '- 环氧乙烷＋氨来自 FlowER 的表示示例；该单例不是 FlowER 作者逐例 DFT 认证。',
        '- `coley_3217` 来自 Coley 组公开的偶极 [3+2] 环加成数据，并有参考 TS 可用于后续独立复核。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    prepare_shared_assets();env=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    html=env.get_template('aimnet2025_paths.html').render(title='AIMNet2-2025 反应路径扩展',root='../',section='experiments',navigation=NAVIGATION,summary=summary)
    (BASE/'index.html').write_text(html,encoding='utf-8')
    source_files=[ROOT/p for p in ['src/mechbridge/reaction_network.py','src/mechbridge/network_aggregation.py','src/mechbridge/search_seeds.py',
        'src/mechbridge/symbolic_library.py','src/mechbridge/intermolecular_actions.py','src/mechbridge/molecular_visuals.py',
        'scripts/exploration/run_network_exploration.py','scripts/exploration/build_aggregate_transitionnet.py',
        'scripts/exploration/render_aimnet2025_reaction_paths.py','scripts/exploration/summarize_aimnet2025_reaction_paths.py',
        'scripts/diagnostics/audit_aimnet2025_reaction_paths.py']]
    evidence=[p for p in BASE.rglob('*') if p.is_file() and 'molecules' not in p.parts and p.suffix.lower() in
              {'.json','.jsonl','.xyz','.md','.npy'} and p.name!='evidence_archive.json']
    files=sorted(set(source_files+evidence));manifest={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    archive=ROOT/'reports/repository_snapshots/aimnet2025_reaction_paths_evidence.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as output:
        for path in files:output.write(path,path.relative_to(ROOT).as_posix())
        output.writestr('MANIFEST.json',json.dumps(manifest,indent=2))
    receipt=dict(file=archive.relative_to(ROOT).as_posix(),files=len(files),bytes=archive.stat().st_size,
        sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    (BASE/'evidence_archive.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    result={k:summary[k] for k in ['completed_runs','attempts','evaluations','edges','unique_graph_pairs','literature_relevant_edges','unique_literature_pairs']};result['evidence_archive']=receipt
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
