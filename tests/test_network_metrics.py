from mechbridge.network_metrics import growth_metrics


def network(graphs,connections):
    return dict(nodes=[dict(id=i,graph_smiles=g) for i,g in enumerate(graphs)],attempts=[],
        edges=[dict(id=i,attempt=i,nodes=[a,b],kind='chemical' if graphs[a]!=graphs[b] else 'conformational')
               for i,(a,b) in enumerate(connections)])


def test_unconnected_same_species_does_not_create_a_path():
    n=network(['A','B','B','C'],[(0,1),(2,3)])
    assert growth_metrics(n)['chains']==[]
    assert growth_metrics(n)['root_species']==['A','B']


def test_actual_conformer_bridge_supports_multistep_path():
    n=network(['A','B','B','C'],[(0,1),(1,2),(2,3)])
    chain=growth_metrics(n)['chains'][0]
    assert chain['nodes']==[0,1,2,3]
    assert chain['chemical_edges']==[0,2]
    assert chain['species']==['A','B','C']


def test_return_to_original_species_is_not_counted_as_growth_chain():
    n=network(['A','B','A','C'],[(0,1),(1,2),(2,3)])
    assert growth_metrics(n)['chains']==[]


def test_enumerated_chains_agree_with_legacy_on_basic_cases():
    n=network(['A','B','B','C'],[(0,1),(1,2),(2,3)])
    metrics=growth_metrics(n)
    assert [c['species'] for c in metrics['enumerated_chains']]==[['A','B','C']]
    assert metrics['enumerated_chains'][0]['nodes']==[0,1,2,3]
    assert growth_metrics(network(['A','B','B','C'],[(0,1),(2,3)]))['enumerated_chains']==[]
    assert growth_metrics(network(['A','B','A','C'],[(0,1),(1,2),(2,3)]))['enumerated_chains']==[]


def test_non_shortest_chain_is_counted():
    # A->B->C->D with a direct A->C shortcut: the shortest path to C hides A-B-C.
    n=network(['A','B','C','D'],[(0,1),(1,2),(0,2),(2,3)])
    metrics=growth_metrics(n)
    species=[c['species'] for c in metrics['enumerated_chains']]
    assert ['A','B','C'] in species and ['A','B','C','D'] in species and ['A','C','D'] in species
    assert metrics['longest_enumerated_chemical_path']==3
    assert metrics['longest_demonstrated_chemical_path']==2


def test_repeated_species_shortcut_does_not_hide_valid_chain():
    # The shortest route to C revisits A; a longer route through D, E, F does not.
    n=network(['A','B','A','C','D','E','F'],[(0,1),(1,2),(2,3),(0,4),(4,5),(5,6),(6,3)])
    metrics=growth_metrics(n)
    assert ['A','D','E','F','C'] in [c['species'] for c in metrics['enumerated_chains']]
    assert metrics['longest_enumerated_chemical_path']==4
    assert not any(c['species'][-1]=='C' for c in metrics['chains'])


def test_chain_completion_uses_earliest_physical_route():
    # Two physical B basins joined by an observed conformer edge; the A-B-C species
    # path is completed by attempts {3,4} via node 2, earlier than {0,9} via node 1.
    n=dict(nodes=[dict(id=i,graph_smiles=g) for i,g in enumerate(['A','B','B','C'])],attempts=[],
           edges=[dict(id=0,attempt=0,nodes=[0,1],kind='chemical'),
                  dict(id=1,attempt=9,nodes=[1,3],kind='chemical'),
                  dict(id=2,attempt=3,nodes=[0,2],kind='chemical'),
                  dict(id=3,attempt=4,nodes=[2,3],kind='chemical'),
                  dict(id=4,attempt=5,nodes=[1,2],kind='conformational')])
    chains=growth_metrics(n)['enumerated_chains']
    chain=next(c for c in chains if c['species']==['A','B','C'])
    assert chain['completed_at_attempt']==4
    assert chain['nodes']==[0,2,3]


def test_enumeration_cap_reports_truncation():
    from mechbridge.network_metrics import enumerate_chains
    graphs=[chr(65+i) for i in range(8)]
    connections=[(i,j) for i in range(8) for j in range(i+1,8)]
    n=network(graphs,connections)
    chains,info=enumerate_chains(n['nodes'],n['edges'],set(range(8)),max_paths=10)
    assert info['truncated'] and info['species_paths']==10
    assert len(chains)==10
    # Exactly max_paths paths is not a truncation; a deeper path is a depth limit.
    line=network(['A','B','C','D'],[(0,1),(1,2),(2,3)])
    chains,info=enumerate_chains(line['nodes'],line['edges'],set(range(4)),max_paths=2)
    assert len(chains)==2 and not info['truncated'] and not info['depth_limited']
    chains,info=enumerate_chains(line['nodes'],line['edges'],set(range(4)),max_steps=2)
    assert [c['species'] for c in chains]==[['A','B','C']] and info['depth_limited']
