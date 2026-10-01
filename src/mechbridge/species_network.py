"""Project physical minima and TS connections onto unique chemical species."""
from copy import deepcopy
import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment


def layout_species(nodes,edges,random_seed=20261001):
    """Assign nonoverlapping slots while minimizing crossings and edge occlusion."""
    count=len(nodes)
    if count<=2:return {i:[.2+.6*i/max(count-1,1),.5] for i in range(count)}
    cols=int(np.ceil(np.sqrt(count*1.4)));rows=int(np.ceil(count/cols));slots=[]
    for row in range(rows):
        present=min(cols,count-row*cols)
        xs=np.linspace(.08,.92,present) if present>1 else np.array([.5])
        y=.12+.76*row/max(rows-1,1)
        slots.extend([[float(x),float(y)] for x in xs])
    pairs=sorted({tuple(sorted(edge['nodes'])) for edge in edges if edge['nodes'][0]!=edge['nodes'][1]})
    def point_segment(point,a,b):
        delta=b-a;t=np.clip(np.dot(point-a,delta)/max(np.dot(delta,delta),1e-12),0,1)
        return np.linalg.norm(point-(a+t*delta))
    def orient(a,b,c):
        ab=b-a;ac=c-a
        return ab[0]*ac[1]-ab[1]*ac[0]
    def score(order):
        positions={node:np.asarray(slots[slot]) for node,slot in enumerate(order)}
        crossings=through=0
        for index,(a,b) in enumerate(pairs):
            p,q=positions[a],positions[b]
            for c,d in pairs[index+1:]:
                if len({a,b,c,d})<4:continue
                r,s=positions[c],positions[d]
                if orient(p,q,r)*orient(p,q,s)<0 and orient(r,s,p)*orient(r,s,q)<0:
                    crossings+=1
            through+=sum(point_segment(positions[k],p,q)<.105 for k in range(count) if k not in (a,b))
        length=sum(np.sum((positions[a]-positions[b])**2) for a,b in pairs)
        return crossings*10000+through*500+length
    rng=np.random.default_rng(random_seed);best=None;best_score=float('inf')
    if count>16:
        graph=nx.Graph();graph.add_nodes_from(range(count));graph.add_edges_from(pairs)
        spring=nx.spring_layout(graph,seed=random_seed,iterations=120)
        target=np.asarray([spring[i] for i in range(count)])
        span=np.ptp(target,axis=0);target=(target-target.min(0))/np.where(span>1e-12,span,1.)
        slot_array=np.asarray(slots)
        _,assignment=linear_sum_assignment(((target[:,None]-slot_array[None,:])**2).sum(2))
        order=assignment.astype(int);current=score(order)
        for _ in range(500):
            a,b=rng.choice(count,2,replace=False);order[a],order[b]=order[b],order[a]
            candidate=score(order)
            if candidate<=current:current=candidate
            else:order[a],order[b]=order[b],order[a]
        return {node:slots[int(order[node])] for node in range(count)}
    restarts=80;swaps=500
    for restart in range(restarts):
        order=np.arange(count) if restart==0 else rng.permutation(count);current=score(order)
        for _ in range(swaps):
            a,b=rng.choice(count,2,replace=False);order[a],order[b]=order[b],order[a]
            candidate=score(order)
            if candidate<=current:current=candidate
            else:order[a],order[b]=order[b],order[a]
        if current<best_score:best,best_score=order.copy(),current
    return {node:slots[int(best[node])] for node in range(count)}


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
