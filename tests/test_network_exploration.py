import json
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from rdkit import Chem
from mechbridge.event_graph import graph_smiles
from mechbridge.symbolic_library import parse_explicit, replay, ArrowLibrary
from mechbridge.search_seeds import internal_direction, make_seed
from mechbridge.reaction_network import CountedCalculator, BudgetExceeded, molecular_rmsd


def example():
    m = parse_explicit('[O:1]=[CH:3][CH2:2][H:4]')
    ids = {a.GetAtomMapNum():a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()}
    arrows = [dict(source=[ids[i] for i in s],sink=[ids[i] for i in t],electrons=2)
              for s,t in [([1],[1,4]),([2,4],[2,3]),([1,3],[1])]]
    return m,arrows


def test_arrow_replay_and_reverse():
    m,a = example()
    p,edits = replay(m,a)
    assert graph_smiles(p) == 'C=CO'
    assert len(edits) == 4
    back,_ = replay(p,[dict(v,source=v['sink'],sink=v['source']) for v in a])
    assert graph_smiles(back) == graph_smiles(m)


def test_budget_counts_actual_calculations_and_cache():
    atoms=Atoms('H2', positions=[[0,0,0],[0,0,.75]])
    calc=CountedCalculator(EMT(), 2)
    atoms.calc=calc
    atoms.get_forces(); atoms.get_potential_energy()
    assert calc.calls==1
    atoms.positions[1,2]+=.01
    atoms.get_forces()
    atoms.positions[1,2]+=.01
    with pytest.raises(BudgetExceeded): atoms.get_forces()
    assert calc.calls==2


def test_internal_direction_removes_rigid_motion():
    x=np.array([[0.,0,0],[1.,0,0],[0.,1.,0]])
    d=internal_direction(x,np.array([[.1,.2,.3],[-.2,.1,0],[.3,0,.2]]))
    assert np.linalg.norm(d)==pytest.approx(1.)
    assert np.linalg.norm(d.sum(0))<1e-10
    assert np.linalg.norm(np.cross(x-x.mean(0),d).sum(0))<1e-10


def test_atom_permutation_and_rotation_node_match():
    m=parse_explicit('O')
    x=np.array([[0.,0,0],[.8,.6,0],[-.8,.6,0]])
    order=[0,2,1]
    other=Chem.RenumberAtoms(m,order)
    rotation=np.array([[0.,-1,0],[1.,0,0],[0.,0,1]])
    assert molecular_rmsd(m,x,other,x[order]@rotation+5)<1e-8


def test_library_uses_public_arrows(tmp_path):
    record=dict(id=17,reaction_name='tautomer',rsmi='[O:1]=[CH:3][CH2:2][H:4]>>[O:1]([CH:3]=[CH2:2])[H:4]',
                epd=[['LP-/Sigma+',[1],[1,4]],['Sigma-/Pi+',[2,4],[2,3]],['Pi-/LP+',[1,3],[1]]])
    path=tmp_path/'polar.json'
    path.write_text(json.dumps(dict(records=[record])))
    library=ArrowLibrary(path)
    m,_=example()
    proposals=library.propose(m)
    assert proposals and proposals[0]['predicted_graph']=='C=CO'
    assert library.propose(parse_explicit('CO'))==[]
    assert library.propose(parse_explicit('C=CO'))
    larger = parse_explicit('CCC=O')
    transferred = library.propose(larger)
    assert transferred and all(p['transferred'] for p in transferred)
    assert any(p['predicted_graph'] == 'CC=CO' for p in transferred)
    assert library.propose(parse_explicit('O=Cc1ccccc1')) == []
    for proposal in transferred:
        product, edits = replay(larger,proposal['arrows'])
        assert product.GetNumAtoms() == larger.GetNumAtoms()
        active = {i for e in edits for i in e['atoms']}
        old = {(b.GetBeginAtomIdx(),b.GetEndAtomIdx()):b.GetBondTypeAsDouble()
               for b in larger.GetBonds() if not ({b.GetBeginAtomIdx(),b.GetEndAtomIdx()} & active)}
        for (i,j),order in old.items():
            assert product.GetBondBetweenAtoms(i,j).GetBondTypeAsDouble() == order
    with pytest.raises(ValueError,match='explicit-H'):
        library.propose(Chem.MolFromSmiles('CCC=O'))


def test_single_electron_rejected():
    m,a=example()
    a[0]['electrons']=1
    with pytest.raises(ValueError): replay(m,a)


def test_missing_or_double_spent_electron_source_rejected():
    m,a=example()
    with pytest.raises(ValueError,match='occupied electron pairs'):
        replay(m,[a[1],a[1]])
    carbon=next(v.GetIdx() for v in m.GetAtoms() if v.GetAtomicNum()==6)
    with pytest.raises(ValueError,match='occupied electron pairs'):
        replay(m,[dict(source=[carbon],sink=[carbon,0],electrons=2)])


def test_all_information_conditions_have_equal_displacement():
    from rdkit.Chem import AllChem
    m,arrows=example()
    assert AllChem.EmbedMolecule(m,randomSeed=17)==0
    x=m.GetConformer().GetPositions()
    atoms=Atoms(numbers=[a.GetAtomicNum() for a in m.GetAtoms()],positions=x)
    product,edits=replay(m,arrows)
    proposal=dict(arrows=arrows,edits=edits,template_id='example',predicted_graph=graph_smiles(product))
    for sample,amplitude in enumerate((.6,1.,1.4)):
        for strategy in ('geometry','center_random','bond_edits','arrows'):
            y,d,meta=make_seed(atoms,m,strategy,proposal,sample,29)
            assert np.linalg.norm(y-x)==pytest.approx(amplitude)
            assert np.linalg.norm(d)==pytest.approx(1.)
            assert np.linalg.norm((y-x).sum(0))<1e-9
            assert np.linalg.norm(d.sum(0))<1e-9
            assert np.linalg.norm(np.cross(y-y.mean(0),d).sum(0))<1e-9
            if strategy in ('bond_edits','arrows'):
                assert meta['angle_targets']
                assert meta['mode_policy']=='damped_internal_coordinate_tangent_at_seed'
                if sample:
                    assert np.ptp(meta['progress_targets']) > 0
            if strategy=='center_random':
                assert 'predicted_graph' not in meta
                assert meta['information']=='active_atom_ids_only'


def test_network_registers_observed_endpoints_not_proposing_node(tmp_path,monkeypatch):
    """Topology-only mock: a seed from A discovering B/C must not invent A/B."""
    import mechbridge.reaction_network as net
    from rdkit.Chem import AllChem
    molecules=[]
    for smiles in ['CC=O','C=CO','C1CO1']:
        mol=parse_explicit(smiles)
        assert AllChem.EmbedMolecule(mol,randomSeed=13)==0
        AllChem.UFFOptimizeMolecule(mol)
        molecules.append((mol,mol.GetConformer().GetPositions()))
    m,x=molecules[0]
    class NoRelax:
        def __init__(self,*args,**kwargs):pass
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def run(self,**kwargs):return True
    monkeypatch.setattr(net,'BFGS',NoRelax)
    monkeypatch.setattr(net,'inspect_point',lambda atoms,protocol: (
        dict(energy_eV=0.,force_converged=True,imaginary_count=0),None))
    def observed(*args):
        ends=[dict(positions_A=xyz.tolist(),graph_smiles=graph_smiles(mol),energy_eV=.1,
                   barrier_eV=.9,force_converged=True,imaginary_count=0)
              for mol,xyz in molecules[1:]]
        return dict(status='validated_descents',evaluations=0,endpoints=ends,
                    ts={'energy_eV':1.},ts_positions_A=x.tolist())
    monkeypatch.setattr(net,'search_connection',observed)
    start=dict(id='topology_mock',atomic_numbers=[a.GetAtomicNum() for a in m.GetAtoms()],
               positions_A=x.tolist(),charge=0,multiplicity=1)
    report=net.explore(start,None,EMT(),'geometry',tmp_path/'network',
                       net.SearchProtocol(max_attempts=1,seeds_per_node=1))
    assert report['edges'][0]['nodes']==[1,2]
    assert report['edges'][0]['proposed_from']==0
    assert report['attempts'][0]['source_is_endpoint'] is False
    assert report['root_component_nodes']==[0]


def test_root_initialization_escapes_negative_curvature_and_charges_cost(tmp_path,monkeypatch):
    """Control-flow test: force convergence cannot hide a high-order stationary point."""
    import mechbridge.reaction_network as net
    class EvaluateOnly:
        def __init__(self,atoms,**kwargs):self.atoms=atoms
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def run(self,**kwargs):self.atoms.get_forces();return True
    checks=[]
    def point(atoms,protocol):
        checks.append(atoms.positions.copy())
        mode=np.zeros_like(atoms.positions);mode[1,2]=1.
        return dict(energy_eV=float(atoms.get_potential_energy()),force_converged=True,
                    imaginary_count=max(0,3-len(checks))),[mode]
    monkeypatch.setattr(net,'BFGS',EvaluateOnly)
    monkeypatch.setattr(net,'inspect_point',point)
    start=dict(id='water_control',atomic_numbers=[8,1,1],
               positions_A=[[0,0,0],[.8,.6,0],[-.8,.6,0]],charge=0,multiplicity=1)
    result=net.explore(start,None,EMT(),'geometry',tmp_path/'root_curvature',
                       net.SearchProtocol(max_attempts=0))
    assert result['status']=='completed'
    assert [r['imaginary_count'] for r in result['initial_relaxation_checks']]==[2,1,0]
    assert result['initialization_evaluations']==result['evaluations']==3
    assert result['initial_point']['imaginary_count']==0
