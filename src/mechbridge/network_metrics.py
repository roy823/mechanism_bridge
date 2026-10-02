"""Growth metrics on actual minimum/TS connectivity, never a merged-SMILES shortcut."""
import networkx as nx

MAX_CHEMICAL_STEPS=10
MAX_SPECIES_PATHS=20000


def basin_components(nodes,edges,connected):
    """Contract only observed non-chemical edges; equal SMILES alone never merges basins."""
    graph=nx.Graph();graph.add_nodes_from(connected)
    graph.add_edges_from(e['nodes'] for e in edges if e['kind']!='chemical' and set(e['nodes'])<=connected)
    components=sorted((sorted(c) for c in nx.connected_components(graph)),key=lambda c:c[0])
    member={i:k for k,c in enumerate(components) for i in c}
    labels=[frozenset(nodes[i]['graph_smiles'] for i in c) for c in components]
    return components,member,labels


def _species_paths(species,root,labels,max_steps,max_paths):
    """Simple component paths from root with no species repeated; >=2 chemical steps.

    Returns (paths, truncated, depth_limited): truncated means more than
    max_paths paths exist; depth_limited means some path could be extended
    beyond max_steps.
    """
    paths=[];stack=[[root]];depth_limited=False
    while stack:
        path=stack.pop()
        if len(path)>=3:
            if len(paths)==max_paths:return paths,True,depth_limited
            paths.append(path)
        seen=frozenset().union(*(labels[k] for k in path))
        following=[k for k in sorted(species[path[-1]],reverse=True) if not labels[k]&seen]
        if len(path)-1>=max_steps:
            depth_limited=depth_limited or bool(following)
            continue
        stack.extend(path+[k] for k in following)
    return paths,False,depth_limited


def _physical_route(path,member,edges,connected):
    """Earliest-completing physical route that follows the component sequence."""
    position={c:k for k,c in enumerate(path)}
    usable=[]
    for e in edges:
        if not set(e['nodes'])<=connected:continue
        a,b=e['nodes'];ca,cb=member[a],member[b]
        if ca not in position or cb not in position:continue
        if e['kind']!='chemical':
            if ca==cb:usable+=[(a,b,e),(b,a,e)]
        elif position[cb]==position[ca]+1:usable.append((a,b,e))
        elif position[ca]==position[cb]+1:usable.append((b,a,e))
    def route(limit):
        graph=nx.DiGraph()
        for a,b,e in usable:
            if e['attempt']>limit:continue
            if not graph.has_edge(a,b) or (e['attempt'],e['id'])<(graph[a][b]['edge']['attempt'],graph[a][b]['edge']['id']):
                graph.add_edge(a,b,edge=e)
        if 0 not in graph:return None
        reached=sorted((len(p),p) for target,p in nx.single_source_shortest_path(graph,0).items()
                       if member[target]==path[-1])
        if not reached:return None
        nodes=reached[0][1]
        return nodes,[graph[a][b]['edge'] for a,b in zip(nodes,nodes[1:])]
    # Reachability is monotone in the attempt limit: bisect for the bottleneck.
    limits=sorted({e['attempt'] for _,_,e in usable})
    low,high,best=0,len(limits)-1,None
    while low<=high:
        middle=(low+high)//2;found=route(limits[middle])
        if found is None:low=middle+1
        else:best,high=found,middle-1
    return best


def enumerate_chains(nodes,edges,connected,max_steps=MAX_CHEMICAL_STEPS,max_paths=MAX_SPECIES_PATHS):
    """Every simple multistep species path from the root, not only shortest physical paths.

    Differs from the legacy shortest-path chains in two ways: longer valid routes
    are counted, and physical minima joined by an observed resonance/conformer
    edge form one basin even when their graph SMILES differ. The second means
    the enumerated longest path is not guaranteed to be >= the legacy value.
    """
    components,member,labels=basin_components(nodes,edges,connected)
    species={k:set() for k in range(len(components))}
    for e in edges:
        if e['kind']=='chemical' and set(e['nodes'])<=connected:
            a,b=(member[i] for i in e['nodes'])
            if a!=b:species[a].add(b);species[b].add(a)
    paths,truncated,depth_limited=_species_paths(species,member[0],labels,max_steps,max_paths)
    chains=[]
    for path in paths:
        realized=_physical_route(path,member,edges,connected)
        if realized is None:continue
        route,route_edges=realized
        chains.append(dict(nodes=route,species=[min(labels[k]) for k in path],
            species_sets=[sorted(labels[k]) for k in path],edges=[e['id'] for e in route_edges],
            chemical_edges=[e['id'] for e in route_edges if e['kind']=='chemical'],
            completed_at_attempt=max(e['attempt'] for e in route_edges)))
    chains.sort(key=lambda c:(c['completed_at_attempt'],len(c['edges']),c['nodes']))
    return chains,enumeration_info(max_steps,max_paths,len(paths),truncated,depth_limited)


def enumeration_info(max_steps,max_paths,species_paths,truncated,depth_limited):
    return dict(policy='simple_species_paths_over_observed_conformer_basins',
        max_chemical_steps=max_steps,max_species_paths=max_paths,species_paths=species_paths,
        truncated=truncated,depth_limited=depth_limited)


def growth_metrics(network):
    nodes=network['nodes'];edges=network['edges']
    graph=nx.Graph();graph.add_nodes_from(range(len(nodes)))
    for edge in edges:graph.add_edge(*edge['nodes'])
    connected=nx.node_connected_component(graph,0) if nodes else set()
    root_graph=nodes[0]['graph_smiles'] if nodes else None
    reachable=[e for e in edges if set(e['nodes'])<=connected]
    pairs={tuple(sorted(nodes[i]['graph_smiles'] for i in e['nodes']))
           for e in reachable if e['kind']=='chemical'}
    # Legacy chains: one shortest physical path per target node. Kept unchanged
    # because historical DFT selections used chains[0]; it misses longer valid
    # routes, so prefer enumerated_chains for path-length statistics.
    chains=[]
    if nodes:
        for target,path in nx.single_source_shortest_path(graph,0).items():
            species=[]
            for i in path:
                smi=nodes[i]['graph_smiles']
                if not species or species[-1]!=smi:species.append(smi)
            if len(species)<3 or len(set(species))!=len(species):continue
            path_edges=[]
            for a,b in zip(path,path[1:]):
                edge=min((e for e in reachable if set(e['nodes'])=={a,b}),key=lambda e:e['attempt'])
                path_edges.append(edge)
            chemicals=[e for e in path_edges if e['kind']=='chemical']
            chains.append(dict(nodes=path,species=species,edges=[e['id'] for e in path_edges],
                chemical_edges=[e['id'] for e in chemicals],
                completed_at_attempt=max(e['attempt'] for e in path_edges)))
    chains.sort(key=lambda c:(c['completed_at_attempt'],len(c['edges']),c['nodes']))
    enumerated,enumeration=(enumerate_chains(nodes,edges,connected) if nodes else
        ([],enumeration_info(MAX_CHEMICAL_STEPS,MAX_SPECIES_PATHS,0,False,False)))
    attempted_species={nodes[a['source_node']]['graph_smiles'] for a in network['attempts']
                       if nodes[a['source_node']]['graph_smiles']!=root_graph}
    return dict(root_nodes=sorted(connected),root_species=sorted({nodes[i]['graph_smiles'] for i in connected}),
        chemical_pairs=sorted(pairs),root_edges=len(reachable),all_edges=len(edges),
        expanded_new_species=sorted(attempted_species),chains=chains,
        longest_demonstrated_chemical_path=max((len(c['species'])-1 for c in chains),
                                             default=1 if pairs else 0),
        enumerated_chains=enumerated,enumerated_chain_count=len(enumerated),chain_enumeration=enumeration,
        longest_enumerated_chemical_path=max((len(c['species'])-1 for c in enumerated),
                                             default=1 if pairs else 0))
