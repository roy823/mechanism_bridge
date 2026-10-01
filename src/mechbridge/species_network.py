"""Project physical minima and TS connections onto unique chemical species."""
from copy import deepcopy


def project_species_network(nodes, edges, root_nodes=()):
    """Return species nodes while retaining every conformer and every TS edge.

    Species identity is the canonical endpoint graph already recorded by the
    physical search. Energy/RMSD deduplication remains a separate physical-node
    operation; this projection only controls the chemical network view.
    """
    species=[]
    by_graph={}
    physical_to_species={}
    for physical in nodes:
        graph=physical.get('graph_smiles',physical.get('smiles'))
        if graph is None:
            raise ValueError('Physical node lacks a recorded chemical graph')
        if graph not in by_graph:
            sid=len(species)
            by_graph[graph]=sid
            species.append(dict(id=sid,graph_smiles=graph,physical_nodes=[],
                representative_node=physical['id'],energy_eV=physical.get('energy_eV',physical.get('energy')),
                depth_discovered=physical.get('depth_discovered')))
        sid=by_graph[graph]
        item=species[sid]
        item['physical_nodes'].append(physical['id'])
        energy=physical.get('energy_eV',physical.get('energy'))
        if energy is not None and (item['energy_eV'] is None or energy<item['energy_eV']):
            item['energy_eV']=energy
            item['representative_node']=physical['id']
        depth=physical.get('depth_discovered')
        if depth is not None and (item['depth_discovered'] is None or depth<item['depth_discovered']):
            item['depth_discovered']=depth
        physical_to_species[physical['id']]=sid
    projected=[]
    pair_counts={}
    for physical_edge in edges:
        edge=deepcopy(physical_edge)
        edge['physical_nodes']=list(physical_edge['nodes'])
        edge['nodes']=[physical_to_species[i] for i in physical_edge['nodes']]
        pair=tuple(sorted(edge['nodes']))
        pair_counts[pair]=pair_counts.get(pair,0)+1
        projected.append(edge)
    used={}
    for edge in projected:
        pair=tuple(sorted(edge['nodes']))
        edge['parallel_index']=used.get(pair,0)
        edge['parallel_count']=pair_counts[pair]
        used[pair]=edge['parallel_index']+1
    species_roots=sorted({physical_to_species[i] for i in root_nodes if i in physical_to_species})
    return dict(nodes=species,edges=projected,root_nodes=species_roots,
                physical_to_species=physical_to_species)
