"""Summarize strict AIMNet2-2025 TransitionNet experiments and literature matches."""
from collections import Counter
from itertools import permutations,product
import hashlib,json,re,zipfile
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol,bond_orders
from mechbridge.reaction_network import aligned_rmsd,atomic_json

BASE=ROOT/'reports/aimnet2025_transitionnet_v9'
REF=ROOT/'data/raw/landscape17/extracted/Landscape17/malonaldehyde/Conformations'


def pi_rmsd(numbers,x,y):
    numbers=np.asarray(numbers);groups=[np.flatnonzero(numbers==z) for z in sorted(set(numbers))];best=float('inf')
    for choices in product(*(permutations(group.tolist()) for group in groups)):
        order=np.empty(len(numbers),dtype=int)
        for group,choice in zip(groups,choices):order[group]=choice
        best=min(best,aligned_rmsd(x,y[order]),aligned_rmsd(-np.asarray(x),y[order]))
    return best


def reference_ts():
    rows=[]
    for index,path in enumerate(sorted(REF.glob('coordsts_*.xyz'))):
        match=re.match(r'coordsts_(\d+)_(\d+)(?:_\d+)?\.xyz',path.name)
        rows.append(dict(id=index,pair=[int(match.group(1)),int(match.group(2))],positions=read(path).positions))
    return rows


def draw_molecule(ax,numbers,positions,center,scale=.14):
    positions=np.asarray(positions);centered=positions-positions.mean(0)
    _,_,vt=np.linalg.svd(centered,full_matrices=False);xy=centered@vt[:2].T
    xy*=scale/max(np.ptp(xy[:,0]),np.ptp(xy[:,1]),1e-6);xy+=np.asarray(center)
    mol=geometry_mol(numbers,positions,0);colors={1:'#e8edf2',6:'#334155',8:'#dc4c4c'}
    for (i,j),order in bond_orders(mol).items():
        ax.plot(xy[[i,j],0],xy[[i,j],1],color='#7b8794',lw=1+order*.45,zorder=1)
    for i,z in enumerate(numbers):
        ax.scatter(*xy[i],s=18 if z==1 else 45,c=colors.get(z,'#888'),edgecolors='white',linewidths=.5,zorder=2)


def main():
    assisted=json.loads((BASE/'landscape17_malonaldehyde_pi/summary.json').read_text())
    pes=json.loads((BASE/'landscape17_malonaldehyde_pes_loop/summary.json').read_text())
    refs=reference_ts();numbers=json.loads((BASE/'landscape17_malonaldehyde_pi/autonomous_start.jsonl').read_text())['atomic_numbers']
    node_ref={r['node']:(r['reference_minimum'] if r['matched'] else None) for r in pes['node_reference_matches']}
    pes_edges=[]
    for connection in pes['connections']:
        if not connection['accepted']:continue
        rank=next(i for i,r in enumerate(connection['refinements']) if r['target_pair'])
        path=BASE/'landscape17_malonaldehyde_pes_loop'/f"pair_{connection['job']:03d}_refine_{rank:02d}"/'result.json'
        result=json.loads(path.read_text());pair=connection['target_nodes'];reference_pair=[node_ref[i] for i in pair]
        candidates=[]
        if None not in reference_pair:
            for ref in refs:
                if sorted(ref['pair'])==sorted(reference_pair):
                    candidates.append((pi_rmsd(numbers,np.asarray(result['ts_positions_A']),ref['positions']),ref['id']))
        match=min(candidates,default=(float('inf'),None))
        pes_edges.append(dict(job=connection['job'],nodes=pair,reference_pair=reference_pair,
            barriers_eV=[e['barrier_eV'] for e in result['endpoints']],reference_ts=match[1],
            reference_ts_rmsd_A=None if not np.isfinite(match[0]) else match[0],
            exact_reference_match=bool(match[0]<.3),result_file=path.relative_to(ROOT).as_posix()))
    autonomous=[]
    roots=[BASE/'landscape17_malonaldehyde_autonomous',BASE/'glyceraldehyde',BASE/'glycolaldehyde',
           BASE/'glycolaldehyde_hydration']
    for root in roots:
        if not root.exists():continue
        for path in root.rglob('network.json'):
            network=json.loads(path.read_text());outcomes=Counter(a['status'] for a in network['attempts'])
            root_graph=network['nodes'][0]['graph_smiles'] if network['nodes'] else None
            strict=[e for e in network['edges'] if e.get('source_connected')]
            physical=sum(all(network['nodes'][i]['graph_smiles']==root_graph for i in e['nodes']) for e in strict)
            autonomous.append(dict(case=network['start']['id'],strategy=network['strategy'],status=network['status'],
                attempts=len(network['attempts']),evaluations=network['evaluations'],nodes=len(network['nodes']),
                strict_edges=len(strict),same_graph_edges=physical,
                detached=outcomes['detached_connection']+len(network['edges'])-len(strict),
                file=path.relative_to(ROOT).as_posix()))
    summary=dict(model='aimnet2-2025',strict_contract='one descent endpoint must match the proposing physical minimum',
        landscape17=dict(reference_minima=2,reference_transition_states=4,
            assisted_recovered=assisted['recovered_transition_states'],assisted_complete=assisted['complete_network_recovered'],
            basin_trials=pes['basin_trials'],model_minima=len(pes['nodes']),matched_reference_minima=sum(r['matched'] for r in pes['node_reference_matches']),
            additional_model_minima=sum(not r['matched'] for r in pes['node_reference_matches']),
            accepted_model_edges=len(pes_edges),exact_reference_edges=sum(e['exact_reference_match'] for e in pes_edges),
            edges=pes_edges,node_matches=pes['node_reference_matches']),autonomous_runs=autonomous,
        conclusions=dict(complete_transitionnet_recovered=False,
            symbol_value='Malonaldehyde arrows found 3 strict edges versus geometry 0, but all 3 changed the reference bond graph and matched 0/4 reference TSs.',
            pes_loop_value='Torsion basin search recovered both reference minima plus one extra AIMNet minimum; NEB/Sella connected all three distinct minima.',
            limiting_layer='Both proposal recall and AIMNet2-2025 PES topology limit reference-network reproduction.'))
    atomic_json(BASE/'summary.json',summary)
    # Molecular network visualization using the actual optimized 3-D node coordinates.
    fig,axes=plt.subplots(1,2,figsize=(11,4.8));ax=axes[0]
    ref_pos={0:(-.65,0),1:(.65,0)}
    ax.set_title('Landscape17 DFT reference: 2 minima / 4 TS')
    for i,p in ref_pos.items():
        draw_molecule(ax,numbers,read(REF/f'coordsmin_{i}.xyz').positions,p,.28);ax.text(p[0],-.43,f'DFT min {i}',ha='center')
    for k,ref in enumerate(refs):
        rad=(-.35,-.16,.16,.35)[k]
        ax.annotate('',xy=ref_pos[1],xytext=ref_pos[0],arrowprops=dict(arrowstyle='-',connectionstyle=f'arc3,rad={rad}',color='#2463eb',lw=1.8))
    ax.set_xlim(-1.2,1.2);ax.set_ylim(-.8,.8);ax.axis('off')
    ax=axes[1];ax.set_title('AIMNet2-2025 PES loop: 3 minima / 3 validated edges')
    positions={0:(0,.62),1:(-.67,-.45),2:(.67,-.45)}
    for node in pes['nodes']:
        p=positions[node['id']];draw_molecule(ax,numbers,node['positions_A'],p,.22)
        match=node_ref[node['id']];ax.text(p[0],p[1]-.34,f"node {node['id']}"+(f" ≈ DFT {match}" if match is not None else ' extra'),ha='center',fontsize=9)
    for edge in pes_edges:
        a,b=[positions[i] for i in edge['nodes']];ax.plot([a[0],b[0]],[a[1],b[1]],color='#ef8b2c',lw=2,zorder=0)
        mid=((a[0]+b[0])/2,(a[1]+b[1])/2);ax.text(*mid,'/'.join(f'{v:.2f}' for v in edge['barriers_eV'])+' eV',fontsize=7,ha='center',bbox=dict(fc='white',ec='none',alpha=.8))
    ax.set_xlim(-1.25,1.25);ax.set_ylim(-1,1);ax.axis('off');fig.tight_layout()
    fig.savefig(BASE/'transitionnet_overview.png',dpi=180,bbox_inches='tight');fig.savefig(BASE/'transitionnet_overview.pdf',bbox_inches='tight');plt.close(fig)
    landscape=summary['landscape17'];lines=['# AIMNet2-2025 复杂体系 TransitionNet 验证 v9','',
        '本轮采用严格根连通判据：只有双侧下降至少一端回到发起物理极小值，事件才进入网络。旧报告中的脱根边不会继续计作增长。','',
        '## Landscape17 malonaldehyde 完整 KTN 对照','',
        f"- DFT 参考：{landscape['reference_minima']} 个极小值、{landscape['reference_transition_states']} 条 TS。",
        f"- 参考 TS 辅助上限：AIMNet2-2025 保留 {landscape['assisted_recovered']}/4 条连接；未得到完整参考网络。",
        f"- 自主 PES loop：32 次二面角盆地搜索得到 {landscape['model_minima']} 个极小值；匹配参考 {landscape['matched_reference_minima']}/2，额外盆地 {landscape['additional_model_minima']} 个。",
        f"- NEB→Sella→双侧下降得到 {landscape['accepted_model_edges']} 条严格模型边，其中精确匹配参考 TS {landscape['exact_reference_edges']}/4。",'',
        '![真实三维节点与网络](transitionnet_overview.png)','',
        '## 自主搜索对照','',
        '| 体系 | 策略 | 状态 | 尝试 | 评估 | 节点 | 严格边 | 同键图边 | 脱根候选 |','|---|---|---|---:|---:|---:|---:|---:|---:|']
    for r in autonomous:lines.append(f"| {r['case']} | {r['strategy']} | {r['status']} | {r['attempts']} | {r['evaluations']} | {r['nodes']} | {r['strict_edges']} | {r['same_graph_edges']} | {r['detached']} |")
    lines += ['', '## 结论','',
        '- 当前流程搭出了一个由 3 个 AIMNet2-2025 极小值和 3 条已验证边构成的局部 TransitionNet，但它不是完整 DFT 参考网络。',
        '- malonaldehyde 箭头搜索相对纯几何搜索增加了严格边数，但三条都是改变键图的高能额外通道，参考 TS 命中为 0/4，不能据此宣称符号机理提高了目标网络复现率。',
        '- PES loop 能恢复两个参考极小值，并暴露一个额外模型盆地；这直接表明势能面拓扑误差会限制网络完整性。',
        '- 下一轮应把“参考 KTN 的节点/边召回率、额外节点/边率”作为主指标，并对候选边做 ensemble/DFT 升级。','',
        '结构化结果：[summary.json](summary.json)。Landscape17：DOI 10.6084/m9.figshare.29949230.v1；参考层级 ωB97x/6-31G(d)。']
    (BASE/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    source_files=[ROOT/p for p in ['src/mechbridge/reaction_network.py','src/mechbridge/symbolic_endpoint.py',
        'src/mechbridge/exploration_actions.py','src/mechbridge/aimnetcentral_backend.py',
        'scripts/exploration/run_network_exploration.py','scripts/exploration/explore_landscape17_malonaldehyde.py',
        'scripts/exploration/summarize_transitionnet_v9.py','scripts/diagnostics/benchmark_landscape17_malonaldehyde.py',
        'scripts/diagnostics/audit_transitionnet_v9.py','scripts/data/fetch_landscape17.py']]
    evidence=[p for p in BASE.rglob('*') if p.is_file() and p.suffix.lower() in
              {'.json','.xyz','.npy','.png','.pdf','.md'} and 'sella_jax_cache' not in p.parts
              and p.name!='evidence_archive.json']
    files=sorted(set(source_files+evidence));manifest={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    archive=ROOT/'reports/repository_snapshots/transitionnet_v9_evidence.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as output:
        for path in files:output.write(path,path.relative_to(ROOT).as_posix())
        output.writestr('MANIFEST.json',json.dumps(manifest,indent=2))
    archive_receipt=dict(file=archive.relative_to(ROOT).as_posix(),files=len(files),
        bytes=archive.stat().st_size,sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    (BASE/'evidence_archive.json').write_text(json.dumps(archive_receipt,indent=2),encoding='utf-8')
    print(json.dumps(dict(landscape=landscape,evidence_archive=archive_receipt),indent=2))


if __name__=='__main__':main()
