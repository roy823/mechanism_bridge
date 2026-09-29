"""Molecular depictions and actual saved geometry/energy trajectories for the pilot.

No new electronic-structure evaluations, coordinate interpolation or TS bond-order
assignment. Dashed contacts are changes between endpoint Lewis graphs, not arrows.
"""
from collections import Counter
from pathlib import Path
import base64
import hashlib
import io
import json
import shutil
import re
import numpy as np
from ase import Atoms
from ase.io import read, write
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from .event_graph import geometry_mol, graph_smiles, bond_orders

ROOT=Path(__file__).resolve().parents[2]
COLORS={1:'#dbe3ef',6:'#475569',7:'#3478d4',8:'#e95058'}
NAMES={'CC=O':'乙醛','C=CO':'乙烯醇','C[N+](=O)[O-]':'硝基甲烷',
       'C=[N+]([O-])O':'aci-nitro 互变体','O=[NH+]C[O-]':'意外连接候选',
       'CC(N)=O':'乙酰胺','O=C1CCC1':'环丁酮','O=CC1CC1':'环收缩产物候选'}
NAMES.update({'CC(C)=O':'丙酮','CCC=O':'丙醛','C=C(C)O':'丙酮烯醇'})
NAMES.update({'O=CCO':'羟基乙醛','O/C=C\\O':'乙烯二醇',
    'O=C[C@H](O)CO':'甘油醛','O=C[C@@H](O)CO':'甘油醛',
    'O=C(CO)CO':'二羟基丙酮','CC(=O)CC(C)=O':'乙酰丙酮',
    'CC[N+](=O)[O-]':'硝基乙烷','OC1=CCC1':'环丁酮烯醇',
    'CC.[C-]#[O+]':'乙烷 + 一氧化碳',
    '[H]/[C-]=[O+]\\CC':'局部电荷中间体候选'})


def depiction(mol, png=False):
    m=Chem.Mol(mol)
    for atom in m.GetAtoms(): atom.SetAtomMapNum(0)
    m=Chem.RemoveHs(m)
    rdDepictor.Compute2DCoords(m)
    drawer=(rdMolDraw2D.MolDraw2DCairo if png else rdMolDraw2D.MolDraw2DSVG)(460,200)
    options=drawer.drawOptions()
    options.bondLineWidth=2.5
    options.padding=.18
    drawer.DrawMolecule(m)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def clean_svg(svg):
    return svg[svg.index('<svg'):]


def entry(mol, coordinates):
    smiles=graph_smiles(mol)
    return dict(smiles=smiles,name=NAMES.get(smiles,'已优化构型'),
        formula=rdMolDescriptors.CalcMolFormula(mol),svg=clean_svg(depiction(mol)),
        positions=np.asarray(coordinates).tolist(),
        bonds=[dict(i=i,j=j,order=order) for (i,j),order in bond_orders(mol).items()])


def display_rotation(ts):
    centered=ts-ts.mean(0)
    _,_,vt=np.linalg.svd(centered)
    rotation=vt.T
    if np.linalg.det(rotation)<0: rotation[:,-1]*=-1
    return rotation


def load_events(run):
    events=[]
    networks=[]
    for path in sorted(run.glob('*/*/network.json')):
        n=json.loads(path.read_text(encoding='utf-8'))
        if n['status']!='completed': raise ValueError('Only completed runs can be visualized')
        numbers=n['start']['atomic_numbers']
        symbols=Atoms(numbers=numbers).get_chemical_symbols()
        network=dict(start=n['start']['id'],strategy=n['strategy'],nodes=[],edges=[],
                     root_nodes=n['root_component_nodes'])
        for node in n['nodes']:
            m=geometry_mol(numbers,node['positions_A'],0)
            network['nodes'].append(dict(id=node['id'],**entry(m,node['positions_A']),
                energy=node['energy_eV']-n['nodes'][0]['energy_eV']))
        for edge in n['edges']:
            attempt=n['attempts'][edge['attempt']]
            folder=path.parent/f"attempt_{edge['attempt']:03d}"
            data=json.loads((folder/'result.json').read_text())
            # Prefer the actual endpoint corresponding to the proposing node on the left.
            left=edge['nodes'].index(attempt['source_node']) if attempt['source_node'] in edge['nodes'] else 0
            right=1-left
            endpoints=data['endpoints']
            ts=np.array(data['ts_positions_A'])
            rotation=display_rotation(ts)
            center=ts.mean(0)
            mols=[geometry_mol(numbers,e['positions_A'],0) for e in endpoints]
            endpoint_entries=[entry(m,e['positions_A']) for m,e in zip(mols,endpoints)]
            branches=[read(folder/f'descent_{sign}.traj',index=':') for sign in (-1,1)]
            frames=[]
            for a in branches[left][::-1]:
                frames.append(dict(positions=a.positions.tolist(),energy=float(a.calc.results['energy']),branch='left'))
            ts_index=len(frames)
            frames.append(dict(positions=ts.tolist(),energy=data['ts']['energy_eV'],branch='ts'))
            for a in branches[right]:
                frames.append(dict(positions=a.positions.tolist(),energy=float(a.calc.results['energy']),branch='right'))
            baseline=endpoints[left]['energy_eV']
            for frame in frames:
                frame['display_positions']=((np.array(frame['positions'])-center)@rotation).round(7).tolist()
                frame['relative_energy']=frame['energy']-baseline
            before,after=bond_orders(mols[left]),bond_orders(mols[right])
            changes=[dict(i=i,j=j,before=before.get((i,j),0),after=after.get((i,j),0))
                     for i,j in sorted(before.keys()|after.keys()) if before.get((i,j),0)!=after.get((i,j),0)]
            for c in changes:
                c['label']=f"{symbols[c['i']]}{c['i']+1}–{symbols[c['j']]}{c['j']+1}"
                c['distances']=[float(np.linalg.norm(np.array(f['positions'])[c['i']]-np.array(f['positions'])[c['j']])) for f in frames]
            pair=[endpoint_entries[left]['smiles'],endpoint_entries[right]['smiles']]
            proposal=attempt['proposal'].get('predicted_graph')
            source=n['nodes'][attempt['source_node']]['graph_smiles']
            matched=sorted(pair)==sorted([source,proposal]) if proposal else None
            eid=f"{n['start']['id']}_{n['strategy']}_{edge['id']}"
            event=dict(id=eid,start=n['start']['id'],strategy=n['strategy'],edge_id=edge['id'],
                node_ids=[edge['nodes'][left],edge['nodes'][right]],symbols=symbols,numbers=numbers,
                left=endpoint_entries[left],right=endpoint_entries[right],frames=frames,ts_index=ts_index,
                changes=changes,ts_energy=data['ts']['energy_eV'],
                barrier_forward=endpoints[left]['barrier_eV'],barrier_reverse=endpoints[right]['barrier_eV'],
                imaginary_frequency=data['ts']['frequencies_cm-1'][0],
                matched_proposal=matched,predicted_graph=proposal,
                root_connected=set(edge['nodes'])<=set(n['root_component_nodes']),
                artifact=folder.relative_to(run).as_posix(),is_IRC=False,DFT_verified=False)
            events.append(event)
            network['edges'].append(dict(id=eid,nodes=edge['nodes'],barriers=edge['barriers_eV']))
        networks.append(network)
    # Start with the clean, familiar chemical example, while exposing every saved edge.
    events.sort(key=lambda e:(e['start']!='MR_8342_1',e['strategy']!='arrows',not e['root_connected'],e['edge_id']))
    return events,networks


def draw_structure(ax, frame, event, bonds, transition=False):
    x=np.array(frame['display_positions'])
    changed={(c['i'],c['j']) for c in event['changes']}
    for b in bonds:
        i,j=b['i'],b['j']
        if transition and tuple(sorted((i,j))) in changed: continue
        mid=(x[i]+x[j])/2
        ax.plot(*np.stack([x[i],mid]).T,color=COLORS[event['numbers'][i]],lw=5,solid_capstyle='round')
        ax.plot(*np.stack([mid,x[j]]).T,color=COLORS[event['numbers'][j]],lw=5,solid_capstyle='round')
    if transition:
        for c in event['changes']:
            p=x[[c['i'],c['j']]]
            ax.plot(*p.T,ls='--',lw=2,color='#d39a28')
    phi=np.linspace(0,np.pi,16);theta=np.linspace(0,2*np.pi,22)
    for i,(z,p) in enumerate(zip(event['numbers'],x)):
        r=.19 if z==1 else .28
        xx=p[0]+r*np.outer(np.sin(phi),np.cos(theta))
        yy=p[1]+r*np.outer(np.sin(phi),np.sin(theta))
        zz=p[2]+r*np.outer(np.cos(phi),np.ones_like(theta))
        ax.plot_surface(xx,yy,zz,color=COLORS[z],linewidth=0,shade=True,antialiased=True)
        ax.text(p[0],p[1],p[2]+r+.12,f"{event['symbols'][i]}{i+1}",fontsize=9,ha='center',color='#21354b')
    extent=max(np.ptp(np.array(f['display_positions']),axis=0).max() for f in event['frames'])/2+.65
    for setter in (ax.set_xlim,ax.set_ylim,ax.set_zlim):setter(-extent,extent)
    ax.set_box_aspect((1,1,1),zoom=1.38);ax.view_init(elev=65,azim=-90);ax.set_axis_off()


def static_figure(event,out):
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    plt.rcParams['pdf.fonttype']=42
    fig=plt.figure(figsize=(14,10),layout='constrained',facecolor='#f6f8fc')
    grid=fig.add_gridspec(4,3,height_ratios=[1.0,2.0,1.2,.18])
    left,right=event['left'],event['right']
    formula_math='$'+re.sub(r'(\d+)',r'_{\1}',left['formula'])+'$'
    fig.suptitle(f"{left['name']} → TS → {right['name']}   |   {formula_math}",fontsize=20,color='#152d48')
    for col,end,frame in [(0,left,event['frames'][0]),(2,right,event['frames'][-1])]:
        ax=fig.add_subplot(grid[0,col]);ax.axis('off')
        m=geometry_mol(event['numbers'],end['positions'],0)
        ax.imshow(plt.imread(io.BytesIO(depiction(m,png=True)),format='png'))
        formula_end='$'+re.sub(r'(\d+)',r'_{\1}',end['formula'])+'$'
        ax.set_title(f"{end['name']}  ·  {formula_end}",fontsize=13)
    ax=fig.add_subplot(grid[0,1]);ax.axis('off')
    ax.text(.5,.65,'→  TS  →',ha='center',fontsize=23,color='#b88724')
    ax.text(.5,.3,f"正向势垒 {event['barrier_forward']:.3f} eV\n逆向势垒 {event['barrier_reverse']:.3f} eV",ha='center',fontsize=12,linespacing=1.8)
    union={tuple(sorted((b['i'],b['j']))):b for b in left['bonds']+right['bonds']}
    for col,index,title,bonds,transition in [
        (0,0,'实际下降端点 A',left['bonds'],False),
        (1,event['ts_index'],'实际优化 TS · 一个虚频',list(union.values()),True),
        (2,len(event['frames'])-1,'实际下降端点 B',right['bonds'],False)]:
        ax=fig.add_subplot(grid[1,col],projection='3d',facecolor='#f6f8fc')
        draw_structure(ax,event['frames'][index],event,bonds,transition)
        ax.set_title(title,fontsize=12,color='#21354b')
    ax=fig.add_subplot(grid[2,:2])
    energies=[f['relative_energy'] for f in event['frames']]
    ax.plot(energies,color='#167d9a',lw=2)
    ax.scatter([event['ts_index']],[energies[event['ts_index']]],color='#d39a28',s=55,zorder=5)
    ax.axvline(event['ts_index'],ls='--',color='#d39a28',alpha=.5)
    ax.set(xlabel='保存的构型帧序号（左侧下降逆序 → TS → 右侧下降）',ylabel='相对端点 A 的势能 / eV')
    ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    ax=fig.add_subplot(grid[2,2]);ax.axis('off')
    changes=event['changes']
    text='关键原子对距离 / Å\n\n'
    for c in changes[:4]:
        d=c['distances']
        text+=f"{c['label']}:  {d[0]:.2f} → {d[event['ts_index']]:.2f} → {d[-1]:.2f}\n"
    text+=f"\nTS 虚频：{event['imaginary_frequency']:.1f} "+r'$\mathrm{cm}^{-1}$'
    ax.text(0,1,text,va='top',fontsize=10,linespacing=1.7)
    ax=fig.add_subplot(grid[3,:]);ax.axis('off')
    ax.text(0,.5,'坐标与能量来自实际 AIMNet2-rxn 计算；黄色虚线表示端点间变化的键。图示为双侧 BFGS 下降，非 IRC / 非动力学时间；独立 DFT 复核见对应报告。',fontsize=9,color='#58677d')
    stem='reaction_'+event['id']
    fig.savefig(out/(stem+'.png'),dpi=180)
    fig.savefig(out/(stem+'.pdf'))
    plt.close(fig)
    return stem


def render_visuals(run):
    run=Path(run).resolve()
    out=run/'molecules';out.mkdir(exist_ok=True)
    events,networks=load_events(run)
    # Keep actual source geometries in downloadable XYZ, SDF and multi-frame XYZ.
    for event in events:
        eid=event['id']
        path=[]
        for f in event['frames']:
            atoms=Atoms(numbers=event['numbers'],positions=f['positions'])
            atoms.info['relative_energy_eV']=f['relative_energy']
            atoms.info['branch']=f['branch'];path.append(atoms)
        write(out/(eid+'_path.xyz'),path,format='extxyz')
        write(out/(eid+'_TS.xyz'),path[event['ts_index']],write_results=False)
        for side in ['left','right']:
            m=geometry_mol(event['numbers'],event[side]['positions'],0)
            (out/(eid+'_'+side+'.mol')).write_text(Chem.MolToMolBlock(m),encoding='utf-8')
    selected=[];seen=set()
    for e in events:
        pair=tuple(sorted([e['left']['smiles'],e['right']['smiles']]))
        if e['root_connected'] and pair not in seen:
            seen.add(pair);selected.append(e)
    figures=[static_figure(e,out) for e in selected]
    assets=ROOT/'assets/molecular_viewer'
    shutil.copy2(assets/'3Dmol-min.js',out/'3Dmol-min.js')
    shutil.copy2(assets/'LICENSE',out/'3Dmol-LICENSE.txt')
    payload=dict(events=events,networks=networks,figures=figures)
    template=(assets/'viewer.html').read_text(encoding='utf-8')
    rendered=template.replace('__MOLECULAR_DATA__',json.dumps(payload,ensure_ascii=False).replace('</','<\\/'))
    (out/'index.html').write_text(rendered,encoding='utf-8')
    overview=run/'index.html'
    if overview.exists():
        source=overview.read_text(encoding='utf-8')
        if 'molecules/index.html' not in source:
            source=source.replace('<h1>Symbol-guided reaction network exploration</h1>',
                '<h1>Symbol-guided reaction network exploration</h1><p><a href="molecules/index.html">分子结构与真实三维反应过程</a></p>')
            overview.write_text(source,encoding='utf-8')
    provenance=dict(event_count=len(events),figures=figures,
        positions='Actual saved unmodified coordinates; display applies one common rigid rotation per event',
        energy='Saved trajectory energies and saved TS result; no new model calls',
        frames='Reversed left BFGS descent + saved TS + forward right BFGS descent; no interpolation',
        bond_display='Endpoint Lewis graphs; changing contacts dashed in TS/path, no TS bond-order claim',
        source_hashes={str(p.relative_to(run)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in run.glob('*/*/attempt_*/result.json') if json.loads(p.read_text())['status']=='validated_descents'})
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    return dict(page=str(out/'index.html'),events=len(events),figures=figures)
