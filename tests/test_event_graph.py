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


def test_resonance_form_that_rdkit_cannot_canonicalize_is_skipped(monkeypatch):
    import mechbridge.event_graph as eg
    from mechbridge.symbolic_library import parse_explicit
    left, right = parse_explicit('CC=O'), parse_explicit('C=CO')
    broken = object()
    real = eg.graph_smiles
    def graph(m):
        if m is broken:
            raise RuntimeError('Invariant Violation')
        return real(m)
    monkeypatch.setattr(eg, 'graph_smiles', graph)
    monkeypatch.setattr(eg.Chem, 'ResonanceMolSupplier', lambda *args, **kwargs: [broken, None])
    assert eg.resonance_equivalent(left, right) is False


def test_graph_smiles_falls_back_when_canonicalization_fails(monkeypatch):
    import mechbridge.event_graph as eg
    from mechbridge.symbolic_library import parse_explicit
    mol = parse_explicit('OC=CO')
    expected = eg.graph_smiles(mol)
    real = eg.Chem.MolToSmiles
    calls = []
    def flaky(m, *args, **kwargs):
        calls.append(kwargs.get('canonical', True))
        if len(calls) == 1:
            raise RuntimeError('Invariant Violation')
        return real(m, *args, **kwargs)
    monkeypatch.setattr(eg.Chem, 'MolToSmiles', flaky)
    assert eg.graph_smiles(mol) == expected
    assert calls == [True, False, True]


def test_resonance_check_skips_forms_reported_as_value_errors(monkeypatch):
    import mechbridge.event_graph as eg
    from mechbridge.symbolic_library import parse_explicit
    left, right = parse_explicit('CC=O'), parse_explicit('C=CO')
    broken = object()
    real = eg.graph_smiles
    def graph(m):
        if m is broken:
            raise ValueError('RDKit could not canonicalize the graph')
        return real(m)
    monkeypatch.setattr(eg, 'graph_smiles', graph)
    monkeypatch.setattr(eg.Chem, 'ResonanceMolSupplier', lambda *args, **kwargs: [broken])
    assert eg.resonance_equivalent(left, right) is False


def test_edge_chemistry_marks_unclassifiable_edges(monkeypatch):
    import mechbridge.reaction_network as net
    from mechbridge.symbolic_library import parse_explicit
    def fail(*args):
        raise ValueError('RDKit could not canonicalize the graph')
    monkeypatch.setattr(net, 'classify_event', fail)
    a, b = parse_explicit('CC=O'), parse_explicit('C=CO')
    fields = net.edge_chemistry([a, b], [a, b])
    assert fields['endpoint_chemistry']['classification'] == 'unclassified'
    assert fields['endpoint_chemistry']['resonance_equivalent'] is False
