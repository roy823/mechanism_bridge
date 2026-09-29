import numpy as np
import pytest
from rdkit import Chem
from mechbridge.event_graph import geometry_mol, graph_smiles, pair_slots, arrow_hypotheses


def test_geometry_graph_keeps_identity_under_coordinate_atom_permutation():
    z=np.array([8,1,1]); xyz=np.array([[0.,0,0],[.96,0,0],[-.24,.93,0]])
    a=geometry_mol(z,xyz,0); perm=[2,0,1]
    b=geometry_mol(z[perm],xyz[perm],0)
    assert graph_smiles(a)==graph_smiles(b)=='O'
    assert [x.GetAtomMapNum() for x in b.GetAtoms()]==[1,2,3]
    assert len(pair_slots(a))==5


def test_unpaired_electron_and_fractional_graph_rejected():
    with pytest.raises(ValueError): pair_slots(Chem.AddHs(Chem.MolFromSmiles('[CH3]')))
    with pytest.raises(ValueError): pair_slots(Chem.AddHs(Chem.MolFromSmiles('c1ccccc1')))


def test_orbital_transport_preserves_electrons_without_claiming_independent_truth():
    mol=Chem.AddHs(Chem.MolFromSmiles('O'))
    slots=pair_slots(mol)
    populations=np.zeros((3,len(slots)))
    for i,site in enumerate(slots): populations[list(site),i]=2/len(site)
    changed=populations.copy(); changed[:,[0,1]]=changed[:,[1,0]]
    result=arrow_hypotheses(populations,changed,mol,mol)
    assert result['electron_pair_bookkeeping_passed']
    assert not result['independently_annotated']
    assert not result['unique_mechanism_certified']
    assert len(result['arrows'])==2


def test_pyscf_analytic_hessian_units_and_axis_order():
    pytest.importorskip('pyscf')
    from ase import Atoms
    from mechbridge.backends import PySCFCalculator
    from mechbridge.physics import finite_hessian
    atoms=Atoms('H2',positions=[[0,0,0],[0,0,.74]])
    atoms.calc=PySCFCalculator(0,1,'HF','sto-3g',1)
    analytic=atoms.calc.hessian(atoms)
    numeric=finite_hessian(atoms,step=.001)
    assert np.allclose(analytic,numeric,rtol=2e-3,atol=2e-3)


def test_unexpected_valid_event_is_preserved_and_failed_descent_is_unresolved():
    from mechbridge.verification import classify_connection
    ts={'force_converged':True,'imaginary_count':1}
    endpoints=[dict(graph_smiles=s,irc_converged=True,force_converged=True,
                    imaginary_count=0,barrier_electronic_eV=1.) for s in ['CCO','COC']]
    result=classify_connection(ts,endpoints,['CCO','CC=O'])
    assert result['physical_event_verified']
    assert result['event_outcome']=='alternative_valid_event'
    assert not result['expected_endpoint_match']
    endpoints[1]['irc_converged']=False
    result=classify_connection(ts,endpoints,['CCO','COC'])
    assert not result['physical_event_verified']
    assert result['event_outcome']=='unresolved'


def test_balanced_atom_relay_convention_retains_evidence():
    from mechbridge.event_graph import contract_atom_relays
    long=[dict(source=[0,3],sink=[0],electrons=2),dict(source=[0],sink=[0,1],electrons=2)]
    short=[dict(source=[0,3],sink=[0,1],electrons=2)]
    result=contract_atom_relays(long)
    assert result['normalized_flows']==contract_atom_relays(short)['normalized_flows']
    assert result['contracted_relays'][0]['relay_atom']==[0]
    assert not result['chemical_equivalence_certified']


def test_same_net_changes_do_not_erase_different_electron_pair_correspondences():
    from mechbridge.event_graph import contract_atom_relays
    a=[dict(source=[0,1],sink=[2,3],electrons=2),dict(source=[4,5],sink=[6,7],electrons=2)]
    b=[dict(source=[0,1],sink=[6,7],electrons=2),dict(source=[4,5],sink=[2,3],electrons=2)]
    assert contract_atom_relays(a)['normalized_flows']!=contract_atom_relays(b)['normalized_flows']
