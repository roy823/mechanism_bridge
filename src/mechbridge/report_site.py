"""Build one offline research site from primary records and actual trajectories."""
import copy,hashlib,json,zipfile
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read,write
from jinja2 import Environment,FileSystemLoader,select_autoescape
from rdkit import Chem,RDLogger
from .report_catalog import ROOT,STAGES,STRATEGIES,catalogue,quantum_checks,read_json
from .report_layout import NAVIGATION,relative_link,prepare_shared_assets,molecular_document,attach_checks
from .molecular_visuals import load_events,entry,display_rotation,render_visuals
from .event_graph import geometry_mol,bond_orders,graph_smiles

CASE_NAMES={'flower_epoxide_ammonia':'FlowER 环氧乙烷＋氨',
    'atlas_glycolaldehyde_hydration':'ReactionAtlas 羟基乙醛水合',
    'atlas_formaldehyde_dimer':'ReactionAtlas 甲醛二聚',
    'flower_test_tropylium_cyanide':'FlowER 测试集：环庚三烯正离子＋氰根'}
SOURCE_LABELS={'paper_illustration_not_individually_DFT_validated':'论文表示示例；无单例 DFT 认证',
    'paper_reports_gas_phase_DFT_barrier_for_same_reaction':'论文报告同一反应的气相 DFT 能垒',
    'paper_reports_alternative_TS_for_same_endpoint_pair':'同一端点对，但作者报告的是另一条 TS 通路',
    'exact_public_test_record_individual_prediction_unavailable':'原始测试记录；逐例预测成功未知'}


def event_from_dft(path):
    q=read_json(path);folder=path.parent;numbers=q['source']['atomic_numbers']
    endpoints=q['endpoints'];ends=[read(folder/f'minimum_{d}.xyz') for d in ('forward','reverse')]
    mols=[geometry_mol(numbers,a.positions,0) for a in ends]
    if [graph_smiles(m) for m in mols]!=[e['graph_smiles'] for e in endpoints]:
        raise ValueError('DFT endpoint coordinates and recorded graph order disagree')
    ts=read(folder/'ts.xyz').positions;frames=[]
    def add(x,energy,branch):frames.append(dict(positions=np.asarray(x).tolist(),energy=float(energy),branch=branch))
    add(ends[0].positions,endpoints[0]['energy_eV'],'left')
    for a in read(folder/'irc_forward.traj',index=':')[::-1]:add(a.positions,a.get_potential_energy(),'left')
    ts_index=len(frames);add(ts,q['ts']['energy_eV'],'ts')
    for a in read(folder/'irc_reverse.traj',index=':'):add(a.positions,a.get_potential_energy(),'right')
    add(ends[1].positions,endpoints[1]['energy_eV'],'right')
    rotation=display_rotation(ts);center=ts.mean(0)
    for f in frames:
        f['relative_energy']=f['energy']-endpoints[0]['energy_eV']
        f['display_positions']=((np.asarray(f['positions'])-center)@rotation).round(7).tolist()
    before,after=[bond_orders(m) for m in mols];symbols=Atoms(numbers=numbers).get_chemical_symbols()
    changes=[]
    for i,j in sorted(before.keys()|after.keys()):
        if before.get((i,j),0)==after.get((i,j),0):continue
        changes.append(dict(i=i,j=j,before=before.get((i,j),0),after=after.get((i,j),0),label=f'{symbols[i]}{i+1}–{symbols[j]}{j+1}',
            distances=[float(np.linalg.norm(np.asarray(f['positions'])[i]-np.asarray(f['positions'])[j])) for f in frames]))
    e=dict(start=q['event_id'],strategy='DFT',edge_id=0,node_ids=[0,1],symbols=symbols,numbers=numbers,
        left=entry(mols[0],ends[0].positions),right=entry(mols[1],ends[1].positions),frames=frames,ts_index=ts_index,
        changes=changes,ts_energy=q['ts']['energy_eV'],barrier_forward=endpoints[0]['barrier_electronic_eV'],
        barrier_reverse=endpoints[1]['barrier_electronic_eV'],imaginary_frequency=q['ts']['frequencies_cm-1'][0],
        matched_proposal=None,matched_reference_pair=q.get('expected_endpoint_match'),predicted_graph=None,root_connected=None,
        artifact=folder.relative_to(ROOT).as_posix(),source_result=path.relative_to(ROOT).as_posix(),
        evidence='DFT_IRC',method_label=q['method']+'/'+q['basis']+' · 气相',
        path_label='双向 IRC 轨迹与端点弛豫构型',is_IRC=True,DFT_verified=True,related_checks=[])
    network=dict(start=e['start'],strategy='DFT',root_nodes=[0,1],nodes=[dict(id=i,**entry(m,a.positions),
        energy=endpoints[i]['energy_eV']-endpoints[0]['energy_eV']) for i,(m,a) in enumerate(zip(mols,ends))],
        edges=[dict(id='',nodes=[0,1],barriers=[e['barrier_forward'],e['barrier_reverse']])])
    return e,network


def export_event(event,network,stage,checks):
    e=copy.deepcopy(event);n=copy.deepcopy(network)
    eid=('d' if e['evidence']=='DFT_IRC' else 'm')+hashlib.sha256((e['source_result']+e['evidence']).encode()).hexdigest()[:14]
    e['id']=eid;out=ROOT/'reports/_events'/eid;out.mkdir(parents=True,exist_ok=True)
    n['edges']=[dict(id=eid,nodes=e['node_ids'],barriers=[e['barrier_forward'],e['barrier_reverse']])]
    # Preserve the full observed node registry but only display this actual edge.
    e['downloads']=dict(path='path.xyz',ts='ts.xyz',left='left.mol',right='right.mol')
    trajectory=[]
    for f in e['frames']:
        a=Atoms(numbers=e['numbers'],positions=f['positions']);a.info['relative_energy_eV']=f['relative_energy'];a.info['branch']=f['branch'];trajectory.append(a)
    write(out/'path.xyz',trajectory,format='extxyz');write(out/'ts.xyz',trajectory[e['ts_index']],write_results=False)
    for side in ('left','right'):
        mol=geometry_mol(e['numbers'],e[side]['positions'],0)
        (out/(side+'.mol')).write_text(Chem.MolToMolBlock(mol),encoding='utf-8')
    attach_checks([e],out,checks)
    stage_url=ROOT/'reports'/stage/'index.html' if stage in {s['id'] for s in STAGES} else ROOT/'reports/data.html'
    payload=dict(events=[e],networks=[n],figures=[],report_url=relative_link(stage_url,out))
    (out/'index.html').write_text(molecular_document(payload,out),encoding='utf-8')
    (out/'provenance.json').write_text(json.dumps(dict(source=e['source_result'],
        sha256=hashlib.sha256((ROOT/e['source_result']).read_bytes()).hexdigest(),evidence=e['evidence'],
        path=e['path_label'],coordinates='Actual saved frames; no interpolation'),indent=2),encoding='utf-8')
    names=[e[s]['name'] if e[s]['name']!='已优化构型' else e[s]['smiles'] for s in ('left','right')]
    return dict(id=eid,title=' ↔ '.join(names),stage=stage,stage_short=next((s['short'] for s in STAGES if s['id']==stage),'历史原型'),
        strategy_label=STRATEGIES[e['strategy']] if e['strategy']!='DFT' else '量化复核',start=e['start'],
        evidence=e['evidence'],evidence_label='DFT / IRC' if e['evidence']=='DFT_IRC' else 'MLIP 双侧下降',
        left_smiles=e['left']['smiles'],right_smiles=e['right']['smiles'],source=e['source_result'],
        case_label=Path(e['source_result']).parent.name if e['evidence']=='DFT_IRC' else e['start']+f" / TS {e['edge_id']}",
        view=relative_link(out/'index.html',ROOT/'reports'),related_checks=e['related_checks'])


def build_site(refresh_legacy=True):
    RDLogger.DisableLog('rdApp.*');reports=ROOT/'reports';prepare_shared_assets()
    originals=[p for p in reports.rglob('index.html') if '_events' not in p.parts and '_site' not in p.parts]
    archive=reports/'repository_snapshots/pre_report_portal_html.zip'
    if not archive.exists():
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
            for p in originals:z.write(p,p.relative_to(ROOT).as_posix())
    stages=catalogue();checks=quantum_checks();events=[];sources=[]
    for stage in stages:
        for run in stage['runs']:
            path=ROOT/run['file'];sources.append(path)
            loaded,networks=load_events(ROOT,[path])
            for event in loaded:events.append(export_event(event,networks[0],stage['id'],checks))
    for q in checks:
        sources.append(ROOT/q['file'])
        if q['passed']:
            event,network=event_from_dft(ROOT/q['file']);events.append(export_event(event,network,q['stage'],checks))
    events.sort(key=lambda e:(e['evidence']!='DFT_IRC',e['stage']!='bimolecular_v5',e['id']))
    if len({e['id'] for e in events})!=len(events):raise ValueError('Duplicate event IDs')
    environment=Environment(loader=FileSystemLoader(ROOT/'assets/report_site'),autoescape=select_autoescape(['html']))
    common=dict(navigation=NAVIGATION,stages=stages,checks=checks,case_names=CASE_NAMES,source_labels=SOURCE_LABELS)
    def render(template,destination,section,**values):
        context=dict(common,root=relative_link(reports,destination.parent)+'/',section=section,**values)
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_text(environment.get_template(template).render(**context),encoding='utf-8')
    render('overview.html',reports/'index.html','overview',title='研究总览',run_count=sum(len(s['runs']) for s in stages),qc_passed=sum(q['passed'] for q in checks))
    plots={'network_exploration_v1':'comparison.png','validation_v2':'budget_comparison.png',
           'network_growth_v4':'growth_overview.png','bimolecular_v5':'bimolecular_overview.png'}
    for stage in stages:
        render('experiment.html',reports/stage['id']/'index.html','experiments',title=stage['title'],
            stage=stage,plot=plots.get(stage['id']),checks=[q for q in checks if q['stage']==stage['id']])
    cases=read_json(reports/'literature_benchmark_v6/cases.json')
    benchmark_path=reports/'literature_benchmark_v6/summary.json'
    aimnetcentral_path=reports/'aimnetcentral_v8/summary.json'
    reaction_paths_path=reports/'aimnet2025_reaction_paths/summary.json'
    render('benchmark.html',reports/'benchmark.html','benchmark',title='文献与基准',cases=cases,
        benchmark=read_json(benchmark_path) if benchmark_path.exists() else None,
        aimnetcentral=read_json(aimnetcentral_path) if aimnetcentral_path.exists() else None,
        reaction_paths=read_json(reaction_paths_path) if reaction_paths_path.exists() else None)
    render('events.html',reports/'events.html','events',title='物理事件',events_json=json.dumps(events,ensure_ascii=False).replace('</','<\\/'))
    render('data.html',reports/'data.html','data',title='数据与复现')
    if refresh_legacy:
        for page in originals:
            if page.parent.name=='molecules':
                sid=page.relative_to(reports).parts[0]
                if sid not in {s['id'] for s in stages}:continue
                render_visuals(page.parent.parent)
            elif page.parent!=reports and page.parent.name not in {s['id'] for s in stages}:
                sid=page.relative_to(reports).parts[0]
                base=next((s for s in stages if s['id']==sid),None)
                if base is None:continue
                subset=[r for r in base['runs'] if (ROOT/r['file']).is_relative_to(page.parent)]
                if not subset:continue
                detail=dict(base,title=base['title']+' / '+page.parent.name,runs=subset,
                    attempts=sum(r['attempts'] for r in subset),evaluations=sum(r['evaluations'] for r in subset))
                render('experiment.html',page,'experiments',title=detail['title'],stage=detail,plot=None,checks=[])
    manifest=dict(primary_runs=sum(len(s['runs']) for s in stages),events=len(events),
        DFT_IRC_events=sum(e['evidence']=='DFT_IRC' for e in events),event_ids=[e['id'] for e in events],
        independent_arrow_gold=0,raw_records_modified=False,legacy_pages_refreshed=refresh_legacy,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (reports/'site_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    (reports/'event_catalog.json').write_text(json.dumps(events,ensure_ascii=False,indent=2),encoding='utf-8')
    return {k:manifest[k] for k in ['primary_runs','events','DFT_IRC_events','raw_records_modified']}
