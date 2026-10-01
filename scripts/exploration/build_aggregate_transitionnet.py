"""Build total same-system TransitionNets across strategies and orientations."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
import networkx as nx
import numpy as np

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol
from mechbridge.molecular_visuals import load_events,entry
from mechbridge.network_aggregation import aggregate_network_records
from mechbridge.report_layout import molecular_document,prepare_shared_assets,relative_link,attach_checks

DEFAULT_BASE=ROOT/'reports/aimnet2025_reaction_paths'


def layout_species(nodes,edges):
    count=len(nodes)
    if count<=2:return {i:[.2+.6*i/max(count-1,1),.5] for i in range(count)}
    cols=int(np.ceil(np.sqrt(count*1.4)));rows=int(np.ceil(count/cols));slots=[]
    for row in range(rows):
        present=min(cols,count-row*cols);xs=np.linspace(.08,.92,present) if present>1 else np.array([.5])
        y=.12+.76*row/max(rows-1,1)
        slots.extend([[float(x),float(y)] for x in xs])
    pairs=sorted({tuple(sorted(e['nodes'])) for e in edges if e['nodes'][0]!=e['nodes'][1]})
    def point_segment(p,a,b):
        delta=b-a;t=np.clip(np.dot(p-a,delta)/max(np.dot(delta,delta),1e-12),0,1)
        return np.linalg.norm(p-(a+t*delta))
    def score(order):
        pos={node:np.asarray(slots[slot]) for node,slot in enumerate(order)};cross=through=0
        for i,(a,b) in enumerate(pairs):
            p,q=pos[a],pos[b]
            for c,d in pairs[i+1:]:
                if len({a,b,c,d})<4:continue
                r,s=pos[c],pos[d]
                orient=lambda x,y,z:np.cross(y-x,z-x)
                if orient(p,q,r)*orient(p,q,s)<0 and orient(r,s,p)*orient(r,s,q)<0:cross+=1
            through+=sum(point_segment(pos[k],p,q)<.105 for k in range(count) if k not in (a,b))
        length=sum(np.sum((pos[a]-pos[b])**2) for a,b in pairs)
        return cross*10000+through*500+length
    rng=np.random.default_rng(20261001);best=None;best_score=float('inf')
    for restart in range(80):
        order=np.arange(count) if restart==0 else rng.permutation(count);current=score(order)
        for _ in range(500):
            a,b=rng.choice(count,2,replace=False);order[a],order[b]=order[b],order[a];candidate=score(order)
            if candidate<=current:current=candidate
            else:order[a],order[b]=order[b],order[a]
        if current<best_score:best,best_score=order.copy(),current
    return {node:slots[int(best[node])] for node in range(count)}


def species_projection(physical_nodes,physical_edges):
    species=[];lookup={};physical_to_species={}
    for node in physical_nodes:
        key=node['smiles']
        if key not in lookup:
            lookup[key]=len(species);species.append(dict(id=len(species),representative=node,physical_nodes=[],energy=node['energy']))
        sid=lookup[key];physical_to_species[node['id']]=sid;species[sid]['physical_nodes'].append(node['id'])
        if node['energy']<species[sid]['energy']:species[sid]['representative']=node;species[sid]['energy']=node['energy']
    nodes=[]
    for item in species:
        rep=dict(item['representative']);rep['id']=item['id'];rep['energy']=item['energy'];rep['conformer_count']=len(item['physical_nodes']);rep['physical_nodes']=item['physical_nodes'];nodes.append(rep)
    edges=[];pair_counts={}
    for edge in physical_edges:
        pair=[physical_to_species[i] for i in edge['nodes']];key=tuple(sorted(pair));pair_counts[key]=pair_counts.get(key,0)+1
        edges.append(dict(edge,nodes=pair))
    used={}
    for edge in edges:
        key=tuple(sorted(edge['nodes']));edge['parallel_index']=used.get(key,0);edge['parallel_count']=pair_counts[key];used[key]=edge['parallel_index']+1
    layout=layout_species(nodes,edges)
    for node in nodes:node['layout']=layout[node['id']]
    return nodes,edges


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=DEFAULT_BASE)
    BASE=parser.parse_args().base.resolve()
    paths=[]
    for path in sorted(BASE.rglob('network.json')):
        data=json.loads(path.read_text(encoding='utf-8'))
        if data['status']=='completed':paths.append((path,data))
    records=[dict(id=path.relative_to(BASE).as_posix(),network=data) for path,data in paths]
    groups=aggregate_network_records(records,group_by='inventory')
    # Load actual saved paths one network at a time and retain their original strategy labels.
    events=[];event_lookup={}
    for path,data in paths:
        run=BASE/path.relative_to(BASE).parts[0]
        loaded,_=load_events(run,[path])
        record=path.relative_to(BASE).as_posix();key=next(k for k,g in groups.items() if (record,0) in g['node_map'])
        group=groups[key]
        for event in loaded:
            item=copy.deepcopy(event);item['network_id']='aggregate:'+group['system']
            item['node_ids']=[group['node_map'][(record,i)] for i in event['node_ids']]
            item['aggregate_edge_id']=group['edge_map'][(record,event['edge_id'])]
            campaign=path.relative_to(BASE).parts[0]
            item['downloads']={kind:f'../../{campaign}/molecules/{event["id"]}_{suffix}' for kind,suffix in
                [('path','path.xyz'),('ts','TS.xyz'),('left','left.mol'),('right','right.mol')]}
            events.append(item);event_lookup[(record,event['edge_id'])]=item['id']
    networks=[];metrics=[]
    for key,group in sorted(groups.items(),key=lambda item:item[1]['system']):
        baseline=min(node['energy_eV'] for node in group['nodes']);nodes=[]
        for node in group['nodes']:
            mol=geometry_mol(node['atomic_numbers'],node['positions_A'],group['charge'])
            nodes.append(dict(id=node['id'],**entry(mol,node['positions_A']),energy=node['energy_eV']-baseline,
                              origins=len(node['origins'])))
        edges=[]
        for edge in group['edges']:
            first=edge['origins'][0];event_id=event_lookup[(first['record'],first['edge'])]
            edges.append(dict(id=edge['id'],label=edge['id'],event_id=event_id,nodes=edge['nodes'],
                barriers=edge['barriers_eV'],origins=edge['origins']))
        graph=nx.Graph();graph.add_nodes_from(range(len(nodes)));graph.add_edges_from(e['nodes'] for e in edges)
        components=nx.number_connected_components(graph) if nodes else 0
        species_nodes,species_edges=species_projection(nodes,edges);species_roots=sorted({next(n['id'] for n in species_nodes if root in n['physical_nodes']) for root in group['root_nodes']})
        networks.append(dict(id='aggregate:'+group['system'],start=group['system'],strategy='aggregate',
            title=f"{group['system']} 固定原子库存 · 物种聚合总网",nodes=species_nodes,edges=species_edges,root_nodes=species_roots))
        metrics.append(dict(system=group['system'],source_systems=group['source_systems'],source_runs=len({o['record'] for n in group['nodes'] for o in n['origins']}),
            raw_nodes=sum(len(record['network']['nodes']) for record in records if (record['id'],0) in group['node_map']),
            merged_nodes=len(nodes),raw_edges=sum(len(record['network']['edges']) for record in records if (record['id'],0) in group['node_map']),
            unique_ts_edges=len(edges),species_nodes=len(species_nodes),species_edges=len(species_edges),
            internal_conformer_ts=sum(e['nodes'][0]==e['nodes'][1] for e in species_edges),components=components,
            root_nodes=group['root_nodes'],parallel_pairs=len(edges)-len({tuple(sorted(e['nodes'])) for e in edges})))
    # Point every event to its total network; duplicate discoveries highlight the deduplicated edge.
    out=BASE/'aggregate/molecules';out.mkdir(parents=True,exist_ok=True);prepare_shared_assets();attach_checks(events,out)
    payload=dict(events=events,networks=networks,figures=[],report_url=relative_link(BASE/'index.html',out))
    (out/'index.html').write_text(molecular_document(payload,out),encoding='utf-8')
    provenance=dict(policy='same model, elemental inventory, charge and multiplicity; physical minima merge across different start species by graph isomorphism, energy and permutation-aware RMSD; TS merge by endpoints, energy and aligned RMSD when atom order agrees',
        geometry_tolerance_A=.15,energy_tolerance_eV=.03,ts_geometry_tolerance_A=.15,ts_energy_tolerance_eV=.03,
        source_hashes={path.relative_to(ROOT).as_posix():hashlib.sha256(path.read_bytes()).hexdigest() for path,_ in paths})
    (out/'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
    summary=dict(systems=metrics,total_source_runs=sum(m['source_runs'] for m in metrics),
        raw_nodes=sum(m['raw_nodes'] for m in metrics),merged_nodes=sum(m['merged_nodes'] for m in metrics),
        raw_edges=sum(m['raw_edges'] for m in metrics),unique_ts_edges=sum(m['unique_ts_edges'] for m in metrics),
        page=(out/'index.html').relative_to(ROOT).as_posix(),policy=provenance)
    continuation=[]
    for group in groups.values():
        for node in group['nodes']:
            if node['id'] in group['root_nodes']:continue
            continuation.append(dict(id=f"{group['system']}__aggregate_n{node['id']}",atomic_numbers=node['atomic_numbers'],
                positions_A=node['positions_A'],charge=group['charge'],multiplicity=group['multiplicity'],
                provenance=dict(aggregate_system=group['system'],aggregate_node=node['id'],
                    source='AIMNet2-2025 aggregate TransitionNet frontier',reference_TS_or_product_geometry_used=False)))
    (BASE/'continuation_starts.jsonl').write_text(''.join(json.dumps(start)+'\n' for start in continuation),encoding='utf-8')
    summary['continuation_starts']=len(continuation)
    (BASE/'aggregate_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
