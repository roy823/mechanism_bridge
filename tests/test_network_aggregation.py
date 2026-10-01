import copy
import numpy as np
from rdkit.Chem import AllChem

from mechbridge.network_aggregation import aggregate_network_records
from mechbridge.species_network import project_species_network
from mechbridge.symbolic_library import parse_explicit


def test_aggregate_merges_duplicate_minima_and_ts_but_keeps_parallel_ts():
    mol=parse_explicit('CC=O');assert AllChem.EmbedMolecule(mol,randomSeed=7)==0
    AllChem.UFFOptimizeMolecule(mol);x=np.asarray(mol.GetConformer().GetPositions());y=x.copy();y[0,2]+=.45
    start=dict(id='case_o0',atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],positions_A=x.tolist(),charge=0,multiplicity=1)
    nodes=[dict(id=0,positions_A=x.tolist(),energy_eV=-10.,graph_smiles='CC=O'),
           dict(id=1,positions_A=y.tolist(),energy_eV=-9.8,graph_smiles='CC=O')]
    edge=dict(id=0,nodes=[0,1],attempt=0,ts_energy_eV=-9.,ts_positions_A=((x+y)/2).tolist(),barriers_eV=[1.,.8])
    network=dict(start=start,nodes=nodes,edges=[edge],root_component_nodes=[0],attempts=[])
    records=[dict(id='a',network=network),dict(id='b',network=copy.deepcopy(network))]
    parallel=copy.deepcopy(network);parallel['edges'][0]['ts_energy_eV']=-8.5
    parallel['edges'][0]['ts_positions_A']=(np.asarray(edge['ts_positions_A'])+np.eye(len(x),3)*.3).tolist()
    records.append(dict(id='c',network=parallel))
    group=next(iter(aggregate_network_records(records).values()))
    assert len(group['nodes'])==2
    assert len(group['edges'])==2
    assert len(group['edges'][0]['origins'])==2
    assert group['root_nodes']==[0]


def test_species_projection_merges_conformers_and_preserves_all_ts_edges():
    nodes=[dict(id=0,graph_smiles='CC=O',energy_eV=-10.,depth_discovered=0),
           dict(id=1,graph_smiles='CC=O',energy_eV=-9.9,depth_discovered=1),
           dict(id=2,graph_smiles='C=CO',energy_eV=-9.8,depth_discovered=1)]
    edges=[dict(id=0,nodes=[0,1]),dict(id=1,nodes=[1,2]),dict(id=2,nodes=[0,2])]
    projected=project_species_network(nodes,edges,[0])
    assert len(projected['nodes'])==2
    assert projected['nodes'][0]['physical_nodes']==[0,1]
    assert projected['nodes'][0]['representative_node']==0
    assert projected['edges'][0]['nodes']==[0,0]
    assert projected['edges'][1]['nodes']==[0,1]
    assert projected['edges'][1]['parallel_count']==2
    assert projected['edges'][2]['parallel_index']==1
    assert projected['root_nodes']==[0]
