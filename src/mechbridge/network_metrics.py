"""Growth metrics on actual minimum/TS connectivity, never a merged-SMILES shortcut."""
import networkx as nx


def growth_metrics(network):
    nodes=network['nodes'];edges=network['edges']
    graph=nx.Graph();graph.add_nodes_from(range(len(nodes)))
    for edge in edges:graph.add_edge(*edge['nodes'])
    connected=nx.node_connected_component(graph,0) if nodes else set()
    root_graph=nodes[0]['graph_smiles'] if nodes else None
    reachable=[e for e in edges if set(e['nodes'])<=connected]
    pairs={tuple(sorted(nodes[i]['graph_smiles'] for i in e['nodes']))
           for e in reachable if e['kind']=='chemical'}
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
    attempted_species={nodes[a['source_node']]['graph_smiles'] for a in network['attempts']
                       if nodes[a['source_node']]['graph_smiles']!=root_graph}
    return dict(root_nodes=sorted(connected),root_species=sorted({nodes[i]['graph_smiles'] for i in connected}),
        chemical_pairs=sorted(pairs),root_edges=len(reachable),all_edges=len(edges),
        expanded_new_species=sorted(attempted_species),chains=chains,
        longest_demonstrated_chemical_path=max((len(c['species'])-1 for c in chains),
                                             default=1 if pairs else 0))
