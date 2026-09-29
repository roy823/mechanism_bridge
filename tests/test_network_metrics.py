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
