"""Build one cross-inventory molecular TransitionNet from v12 and v13 evidence."""
import copy
import hashlib
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.io import write
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol
from mechbridge.molecular_visuals import load_events
from mechbridge.report_layout import attach_checks,molecular_document,prepare_shared_assets,relative_link
from mechbridge.species_network import layout_species

BASE=ROOT/'reports/formose_closure_v13'
V12=ROOT/'reports/formose_cycle_v12'


def selected_paths():
    summary=json.loads((V12/'summary.json').read_text(encoding='utf-8'))
    paths=[]
    for relative in summary['selected_networks'].values():
        path=ROOT/relative;network=json.loads(path.read_text(encoding='utf-8'))
        if network['status']=='completed':paths.append(path)
    for path in sorted((BASE/'runs').rglob('network.json')):
        network=json.loads(path.read_text(encoding='utf-8'))
        if network['edges']:paths.append(path)
    return paths


def main():
    out=BASE/'global_transitionnet/molecules';out.mkdir(parents=True,exist_ok=True)
    events=[]
    for path in selected_paths():
        run=path.parents[2]
        loaded,_=load_events(run,[path],layout_network=False)
        events.extend(loaded)
    nodes=[];node_by_key={}
    for event in events:
        for side in ('left','right'):
            endpoint=event[side];key=endpoint['smiles']
            if key not in node_by_key:
                node_by_key[key]=len(nodes)
                node=dict(endpoint,id=len(nodes),smiles=key,energy=0.,energy_label='跨库存 · 不比较绝对能量')
                nodes.append(node)
    edges=[];pair_counts={}
    for event in events:
        pair=[node_by_key[event[side]['smiles']] for side in ('left','right')]
        event['physical_node_ids']=event.get('physical_node_ids',event['node_ids'])
        event['node_ids']=pair;event['network_id']='formose-global'
        event['aggregate_edge_id']=event['id']
        key=tuple(sorted(pair));pair_counts[key]=pair_counts.get(key,0)+1
        edges.append(dict(id=event['id'],label=len(edges),event_id=event['id'],nodes=pair,
                          barriers=[event['barrier_forward'],event['barrier_reverse']]))
    used={}
    for edge in edges:
        key=tuple(sorted(edge['nodes']));edge['parallel_index']=used.get(key,0)
        edge['parallel_count']=pair_counts[key];used[key]=edge['parallel_index']+1
    layout=layout_species(nodes,edges)
    for node in nodes:node['layout']=layout[node['id']]
    for event in events:
        eid=event['id'];trajectory=[]
        for frame in event['frames']:
            atoms=Atoms(numbers=event['numbers'],positions=frame['positions'])
            atoms.info['relative_energy_eV']=frame['relative_energy'];atoms.info['branch']=frame['branch']
            trajectory.append(atoms)
        write(out/(eid+'_path.xyz'),trajectory,format='extxyz')
        write(out/(eid+'_TS.xyz'),trajectory[event['ts_index']],write_results=False)
        for side in ('left','right'):
            mol=geometry_mol(event['numbers'],event[side]['positions'],0)
            (out/(eid+'_'+side+'.mol')).write_text(Chem.MolToMolBlock(mol),encoding='utf-8')
        event['downloads']=dict(path=eid+'_path.xyz',ts=eid+'_TS.xyz',left=eid+'_left.mol',right=eid+'_right.mol')
    prepare_shared_assets();attach_checks(events,out)
    network=dict(id='formose-global',start='C2/C3/C4 formose',strategy='aggregate',
                 title='C2/C3/C4 全物种 formose TransitionNet',nodes=nodes,edges=edges,root_nodes=[])
    payload=dict(events=events,networks=[network],figures=[],report_url=relative_link(BASE/'index.html',out))
    (out/'index.html').write_text(molecular_document(payload,out),encoding='utf-8')
    provenance=dict(networks=len(selected_paths()),events=len(events),species_nodes=len(nodes),edges=len(edges),
        energy_policy='Cross-inventory nodes do not share an energy zero; node energies intentionally hidden',
        species_policy='Canonical stereo-resolved endpoint graph; constitutional cycle checks are reported separately',
        source_hashes={path.relative_to(ROOT).as_posix():hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in selected_paths()})
    (out/'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(provenance,ensure_ascii=True,indent=2))


if __name__=='__main__':main()
