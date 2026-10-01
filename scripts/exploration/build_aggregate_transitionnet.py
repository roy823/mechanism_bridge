"""Build total same-system TransitionNets across strategies and orientations."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
import networkx as nx
import numpy as np

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol
from mechbridge.molecular_visuals import load_events,entry
from mechbridge.network_aggregation import aggregate_network_records
from mechbridge.species_network import layout_species,project_species_network
from mechbridge.report_layout import molecular_document,prepare_shared_assets,relative_link,attach_checks

DEFAULT_BASE=ROOT/'reports/aimnet2025_reaction_paths'


def species_projection(physical_nodes,physical_edges):
    projection=project_species_network(physical_nodes,physical_edges)
    physical_by_id={node['id']:node for node in physical_nodes}
    nodes=[]
    for item in projection['nodes']:
        rep=dict(physical_by_id[item['representative_node']]);rep['id']=item['id']
        rep['energy']=item['energy_eV'];rep['conformer_count']=len(item['physical_nodes'])
        rep['physical_nodes']=item['physical_nodes'];nodes.append(rep)
    edges=projection['edges']
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
        loaded,_=load_events(run,[path],layout_network=False)
        record=path.relative_to(BASE).as_posix();key=next(k for k,g in groups.items() if (record,0) in g['node_map'])
        group=groups[key]
        for event in loaded:
            item=copy.deepcopy(event);item['network_id']='aggregate:'+group['system']
            physical_ids=event.get('physical_node_ids',event['node_ids'])
            item['node_ids']=[group['node_map'][(record,i)] for i in physical_ids]
            item['physical_node_ids']=physical_ids
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
