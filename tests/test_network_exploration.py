import json
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
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
        dict(energy_eV=0.,force_converged=True,imaginary_count=0,
             force_norms_eV_A=[0.]*len(atoms)),None))
    def observed(*args):
        ends=[dict(positions_A=xyz.tolist(),graph_smiles=graph_smiles(mol),energy_eV=.1,
                   barrier_eV=.9,force_converged=True,imaginary_count=0)
              for mol,xyz in molecules[1:]]
        return dict(status='validated_descents',evaluations=0,seconds=0.,endpoints=ends,
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
                    imaginary_count=max(0,3-len(checks)),
                    force_norms_eV_A=[0.]*len(atoms)),[mode]
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


def test_batched_hessian_counts_geometries_and_preserves_positions():
    from ase.calculators.calculator import Calculator,all_changes
    from mechbridge.physics import finite_hessian
    class Harmonic(Calculator):
        implemented_properties=['energy','forces']
        def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            x=self.atoms.positions
            self.results=dict(energy=float((x*x).sum()/2),forces=-x.copy())
        def evaluate_many(self,numbers,positions):
            return dict(energy=(positions**2).sum((1,2))/2,forces=-positions)
    atoms=Atoms('H2',positions=[[0.,0.,0.],[0.,0.,.8]])
    atoms.calc=CountedCalculator(Harmonic(),13)
    atoms.get_forces()
    before=atoms.positions.copy()
    h=finite_hessian(atoms,batch_forces=atoms.calc.batch_forces,batch_size=4)
    assert np.allclose(h,np.eye(6))
    assert np.array_equal(before,atoms.positions)
    assert atoms.calc.calls==13
    assert atoms.calc.model_calls==4
    with pytest.raises(BudgetExceeded):
        atoms.calc.batch_forces(atoms.numbers,np.array([before]))


def test_action_scheduler_diversity_and_finite_exhaustion():
    from mechbridge.exploration_actions import choose_action,action_key
    proposals=[dict(edits=[i],arrows=[i],predicted_graph=g) for i,g in enumerate(['A','A','B'])]
    used={}
    p,key,variant,strategy=choose_action(proposals,used,{'R'},'arrows')
    assert p['predicted_graph']=='A' and variant==0
    used[key]=1
    assert choose_action(proposals,used,{'R'},'arrows')[0]['predicted_graph']=='B'
    assert choose_action(proposals,{}, {'A'},'arrows')[0]['predicted_graph']=='B'
    assert choose_action(proposals,{action_key(p):3 for p in proposals},set(),'arrows') is None
    assert choose_action([],{},set(),'hybrid')[-1]=='geometry'


def test_parallel_reservation_identity_is_rigid_invariant_and_variant_specific():
    from mechbridge.parallel_network import reservation_identity
    node=dict(graph_smiles='CC=O',positions_A=[[0.,0.,0.],[1.2,0.,0.],[1.8,.8,0.]])
    proposal=dict(edits=[dict(atoms=[0,1],before=1,after=2)],
                  arrows=[dict(source=[0],sink=[0,1],electrons=2)])
    first,_=reservation_identity(node,proposal,0,'arrows','species-c0')
    rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    moved=dict(node,positions_A=(np.asarray(node['positions_A'])@rotation+4).tolist())
    same,_=reservation_identity(moved,proposal,0,'arrows','species-c0')
    distorted=dict(node,positions_A=[[0.,0.,0.],[1.4,0.,0.],[2.2,.9,0.]])
    same_cluster,_=reservation_identity(distorted,proposal,0,'arrows','species-c0')
    other_conformer,_=reservation_identity(distorted,proposal,0,'arrows','species-c1')
    other,_=reservation_identity(node,proposal,1,'arrows','species-c0')
    assert first==same
    assert first==same_cluster
    assert first!=other_conformer
    assert first!=other


def test_dimer_stops_on_physical_force_not_projected_norm():
    from mechbridge.saddle_optimization import StationaryDimerTranslate
    class Modes:
        force=np.array([[.0051,0.,0.]])
        curvature=-1.
        def get_forces(self,real=False):
            assert real
            return self.force
        def get_curvature(self):return self.curvature
    optimizer=object.__new__(StationaryDimerTranslate)
    optimizer.dimeratoms=Modes();optimizer.fmax=.005
    assert not optimizer.gradient_converged(np.array([.0049,0.,0.]))
    optimizer.dimeratoms.force=np.array([[.0049,0.,0.]])
    assert optimizer.gradient_converged(np.array([.0051,0.,0.]))
    optimizer.dimeratoms.curvature=1.
    assert not optimizer.gradient_converged(np.zeros(3))


def test_recoverable_failure_classification():
    from mechbridge.reaction_network import is_recoverable_failure
    from concurrent.futures.process import BrokenProcessPool
    for exc in (ValueError('valence'),np.linalg.LinAlgError('singular'),RuntimeError('diverged'),
                ZeroDivisionError('step'),BudgetExceeded('budget')):
        assert is_recoverable_failure(exc)
    for exc in (TypeError('bug'),KeyError('field'),AttributeError('bug'),IndexError('bug'),
                OSError('disk'),AssertionError('invariant'),NotImplementedError('missing'),
                RecursionError('loop'),BrokenProcessPool('worker died'),
                RuntimeError('CUDA out of memory. Tried to allocate 2 GiB')):
        assert not is_recoverable_failure(exc)


@pytest.mark.parametrize('error,status,raises',[(TypeError,'fatal_error',True),
                                                (ValueError,'calculation_failed',False)])
def test_search_connection_separates_bugs_from_failed_chemistry(tmp_path,monkeypatch,error,status,raises):
    import mechbridge.reaction_network as net
    class Broken:
        def __init__(self,*args,**kwargs):raise error('injected')
    monkeypatch.setattr(net,'DimerControl',Broken)
    seed=Atoms('H2',positions=[[0.,0.,0.],[0.,0.,.75]])
    calc=CountedCalculator(EMT(),100)
    run=lambda: net.search_connection(seed,np.zeros((2,3)),calc,tmp_path/'attempt',net.SearchProtocol())
    if raises:
        with pytest.raises(error):run()
    else:
        assert run()['status']==status
    result=json.loads((tmp_path/'attempt'/'result.json').read_text())
    assert result['status']==status
    assert result['error']==f'{error.__name__}: injected'
    assert (tmp_path/'attempt'/'error.log').exists()


def test_connection_status_barrier_gate_and_acceptance():
    from mechbridge.reaction_network import SearchProtocol,connection_status
    def end(barrier,full=True):
        return dict(full_system_minimum=full,carbon_skeleton_force_converged=True,barrier_eV=barrier)
    legacy,strict=SearchProtocol(),SearchProtocol(min_barrier_eV=1e-3)
    # The legacy gate accepts an endpoint 0.05 meV above the TS; the strict gate does not.
    assert connection_status([end(.5),end(-5e-5)],legacy)=='validated_descents'
    assert connection_status([end(.5),end(-5e-5)],strict)=='unresolved_minimum'
    assert connection_status([end(.5),end(5e-4)],strict)=='unresolved_minimum'
    assert connection_status([end(.5),end(2e-3)],strict)=='validated_descents'
    core=[end(.5,False),end(.4,False)]
    assert connection_status(core,legacy)=='unresolved_minimum'
    assert connection_status(core,SearchProtocol(endpoint_acceptance='carbon_skeleton'))=='validated_core_descents'
    with pytest.raises(ValueError):SearchProtocol(encounter_policy='unknown')
    with pytest.raises(ValueError):SearchProtocol(endpoint_acceptance='unknown')


def test_edge_chemistry_uses_attempt_atom_frame(tmp_path,monkeypatch):
    """A node matched only up to a fragment swap must not supply the bond indices."""
    import mechbridge.reaction_network as net
    from rdkit.Chem import AllChem
    monomer=parse_explicit('C=O')
    assert AllChem.EmbedMolecule(monomer,randomSeed=5)==0
    xyz=monomer.GetConformer().GetPositions()
    numbers=[a.GetAtomicNum() for a in monomer.GetAtoms()]*2
    shift=np.array([5.,0.,0.])
    x=np.vstack([xyz,xyz+shift])
    swapped=np.vstack([xyz+shift,xyz])     # same nuclei, fragment labels exchanged
    apart=np.vstack([xyz,xyz+2*shift])
    class NoRelax:
        def __init__(self,*args,**kwargs):pass
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def run(self,**kwargs):return True
    monkeypatch.setattr(net,'BFGS',NoRelax)
    monkeypatch.setattr(net,'inspect_point',lambda atoms,protocol: (
        dict(energy_eV=0.,force_converged=True,imaginary_count=0,
             force_norms_eV_A=[0.]*len(atoms)),None))
    def observed(*args):
        ends=[dict(positions_A=y.tolist(),graph_smiles='C=O.C=O',energy_eV=e,barrier_eV=1.-e,
                   force_converged=True,imaginary_count=0) for y,e in ((swapped,0.),(apart,.1))]
        return dict(status='validated_descents',evaluations=0,seconds=0.,endpoints=ends,
                    ts={'energy_eV':1.},ts_positions_A=x.tolist())
    monkeypatch.setattr(net,'search_connection',observed)
    seen=[]
    real=net.classify_event
    def capture(left,right):
        seen.append((left.GetConformer().GetPositions(),right.GetConformer().GetPositions()))
        return real(left,right)
    monkeypatch.setattr(net,'classify_event',capture)
    start=dict(id='fragment_swap',atomic_numbers=numbers,positions_A=x.tolist(),charge=0,multiplicity=1)
    report=net.explore(start,None,EMT(),'geometry',tmp_path/'network',
                       net.SearchProtocol(max_attempts=1,seeds_per_node=1))
    assert report['attempts'][0]['observed_nodes'][0]==0      # swapped copy matched the root
    assert report['edges'][0]['endpoint_chemistry_atom_frame']=='attempt'
    left,right=seen[0]
    assert np.allclose(left,swapped) and np.allclose(right,apart)


def test_worker_failure_is_charged_from_artifact_or_reservation(tmp_path):
    from concurrent.futures.process import BrokenProcessPool
    from mechbridge.parallel_network import worker_failure_evaluations
    artifact=tmp_path/'result.json'
    assert worker_failure_evaluations(artifact,TypeError('bug'),700)==(0,'no_artifact_before_search')
    assert worker_failure_evaluations(artifact,BrokenProcessPool('killed'),700)==(700,'reserved_upper_bound')
    artifact.write_text(json.dumps(dict(status='fatal_error',evaluations=42)))
    assert worker_failure_evaluations(artifact,TypeError('bug'),700)==(42,'artifact')


def test_explore_shared_uses_explicit_start_method_and_records_abort(tmp_path,monkeypatch):
    import mechbridge.parallel_network as par
    from mechbridge.reaction_network import SearchProtocol
    class Stop(Exception):pass
    captured={}
    class FakePool:
        def __init__(self,**kwargs):
            captured.update(kwargs);raise Stop('pool not started in this test')
    def root(start,calculator,outdir,protocol):
        info=dict(energy_eV=0.,force_converged=True,imaginary_count=0,full_system_minimum=True,
                  carbon_skeleton_force_converged=True,positions_A=start['positions_A'],graph_smiles='O')
        return info,[dict(info)],0
    monkeypatch.setattr(par,'ProcessPoolExecutor',FakePool)
    monkeypatch.setattr(par,'initialize_root',root)
    start=dict(id='water',atomic_numbers=[8,1,1],positions_A=[[0.,0.,0.],[.96,0.,0.],[-.24,.93,0.]],
               charge=0,multiplicity=1)
    with pytest.raises(ValueError,match='start method'):
        par.explore_shared(start,None,EMT(),'geometry',tmp_path/'bad','aimnet2-2025',tmp_path,
                           SearchProtocol(),workers=2,start_method='threads')
    with pytest.raises(Stop):
        par.explore_shared(start,None,EMT(),'geometry',tmp_path/'run','aimnet2-2025',tmp_path,
                           SearchProtocol(),workers=2)
    assert captured['mp_context'].get_start_method()=='spawn'
    report=json.loads((tmp_path/'run'/'network.json').read_text())
    assert report['status']=='aborted_error'
    assert report['error'].startswith('Stop')
    assert report['scheduler']['start_method']=='spawn'


class _NoOpContext:
    def __init__(self,*args,**kwargs):pass
    def __enter__(self):return self
    def __exit__(self,*args):return False


class _NoRelax(_NoOpContext):
    def run(self,**kwargs):return True


class _Harmonic(Calculator):
    """Analytic well centred on given positions."""
    implemented_properties=['energy','forces']
    def __init__(self,center):
        super().__init__();self.center=np.asarray(center,dtype=float)
    def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        d=self.atoms.positions-self.center
        self.results=dict(energy=float((d*d).sum()/2),forces=-d)


def _mock_search_engines(monkeypatch,calls,ts_imaginary=0):
    import mechbridge.reaction_network as net
    class Dimer(_NoOpContext):
        def run(self,fmax,steps):calls.append(('Dimer.run',fmax,steps));return True
    monkeypatch.setattr(net,'DimerControl',_NoOpContext)
    monkeypatch.setattr(net,'MinModeAtoms',lambda *args,**kwargs:None)
    monkeypatch.setattr(net,'StationaryDimerTranslate',Dimer)
    monkeypatch.setattr(net,'BFGS',_NoRelax)
    points=[]
    def inspect(atoms,protocol):
        points.append(atoms.positions.copy())
        first=len(points)==1
        mode=np.zeros_like(atoms.positions);mode[0,0]=1.
        return (dict(energy_eV=1. if first else .2,force_converged=True,
                     imaginary_count=ts_imaginary if first else 0,force_max_eV_A=0.,
                     force_norms_eV_A=[0.]*len(atoms)),[mode])
    monkeypatch.setattr(net,'inspect_point',inspect)
    return net


class _InnerLoopFailure(RuntimeError):
    pass


def _fake_sella(calls,fail_directions=()):
    import types
    class Sella(_NoOpContext):
        def __init__(self,atoms,**kwargs):calls.append(('Sella',kwargs))
        def run(self,fmax,steps):calls.append(('Sella.run',fmax,steps));return True
    class IRC(_NoOpContext):
        def __init__(self,atoms,**kwargs):calls.append(('IRC',kwargs));self.atoms=atoms
        def run(self,fmax,fmax_inner,steps,direction):
            calls.append(('IRC.run',direction,fmax,fmax_inner,steps));self.atoms.get_forces()
            if direction in fail_directions:raise _InnerLoopFailure
            return True
    return types.SimpleNamespace(Sella=Sella,IRC=IRC)


def _install_fake_sella(monkeypatch,calls,fail_directions=()):
    import sys,types
    monkeypatch.setitem(sys.modules,'sella',_fake_sella(calls,fail_directions))
    monkeypatch.setitem(sys.modules,'sella.optimize',types.SimpleNamespace())
    monkeypatch.setitem(sys.modules,'sella.optimize.irc',
                        types.SimpleNamespace(IRCInnerLoopConvergenceFailure=_InnerLoopFailure))


@pytest.mark.parametrize('optimizer',['dimer','dimer+sella'])
def test_ts_optimizer_stages(tmp_path,monkeypatch,optimizer):
    import sys
    calls=[]
    net=_mock_search_engines(monkeypatch,calls)
    if optimizer=='dimer+sella':_install_fake_sella(monkeypatch,calls)
    else:monkeypatch.setitem(sys.modules,'sella',None)
    x=np.array([[0.,0.,0.],[0.,0.,.75]])
    seed=Atoms('H2',positions=x)
    calc=CountedCalculator(_Harmonic(x),1000)
    protocol=net.SearchProtocol(ts_optimizer=optimizer)
    result=net.search_connection(seed,np.zeros((2,3)),calc,tmp_path/'attempt',protocol)
    assert result['status']=='not_index_one'
    if optimizer=='dimer':
        assert calls==[('Dimer.run',protocol.fmax,protocol.ts_steps)]
        assert 'ts_optimization' not in result
    else:
        assert calls[0]==('Dimer.run',protocol.dimer_fmax,protocol.ts_steps)
        assert calls[1]==('Sella',dict(order=1,internal=False,logfile=str(tmp_path/'attempt'/'sella.log'),
                                        trajectory=str(tmp_path/'attempt'/'sella.traj')))
        assert calls[2]==('Sella.run',protocol.fmax,protocol.sella_steps)
        assert result['ts_optimization']['sella_run'] is True


def test_sella_runs_only_after_the_dimer_gate(tmp_path,monkeypatch):
    import sys
    calls=[]
    net=_mock_search_engines(monkeypatch,calls)
    _install_fake_sella(monkeypatch,calls)
    x=np.array([[0.,0.,0.],[0.,0.,.75]])
    seed=Atoms('H2',positions=x+np.array([[.1,0.,0.],[0.,0.,0.]]))   # residual force 0.1 > 0.05
    result=net.search_connection(seed,np.zeros((2,3)),CountedCalculator(_Harmonic(x),1000),
                                 tmp_path/'attempt',net.SearchProtocol(ts_optimizer='dimer+sella'))
    assert result['status']=='ts_force_unconverged'
    assert [c[0] for c in calls]==['Dimer.run']
    assert result['ts_optimization']['sella_run'] is False


def test_irc_connection_protocol(tmp_path,monkeypatch):
    import sys
    calls=[]
    net=_mock_search_engines(monkeypatch,calls,ts_imaginary=1)
    _install_fake_sella(monkeypatch,calls)
    # Water, because the endpoints go through real RDKit bond perception.
    x=np.array([[0.,0.,0.],[.96,0.,0.],[-.24,.93,0.]])
    protocol=net.SearchProtocol(connection_protocol='irc')
    result=net.search_connection(Atoms('OH2',positions=x),np.zeros((3,3)),
                                 CountedCalculator(_Harmonic(x),1000),tmp_path/'attempt',protocol)
    runs=[c for c in calls if c[0]=='IRC.run']
    assert [c[1] for c in runs]==['reverse','forward']
    assert all(c[2:]==(protocol.irc_fmax,protocol.irc_inner_fmax,protocol.irc_steps) for c in runs)
    assert all(c[1]['dx']==protocol.irc_dx and c[1]['ninner_iter']==20 for c in calls if c[0]=='IRC')
    assert [e['irc_direction'] for e in result['endpoints']]==['reverse','forward']
    assert result['is_IRC'] is True and result['evidence']=='MLIP_bidirectional_IRC'
    assert result['status']=='validated_descents'
    with pytest.raises(ValueError):net.SearchProtocol(connection_protocol='neb')
    with pytest.raises(ValueError):net.SearchProtocol(ts_optimizer='sella')


def test_irc_inner_loop_failure_is_an_unconverged_branch(tmp_path,monkeypatch):
    calls=[]
    net=_mock_search_engines(monkeypatch,calls,ts_imaginary=1)
    _install_fake_sella(monkeypatch,calls,fail_directions=('forward',))
    x=np.array([[0.,0.,0.],[.96,0.,0.],[-.24,.93,0.]])
    result=net.search_connection(Atoms('OH2',positions=x),np.zeros((3,3)),
                                 CountedCalculator(_Harmonic(x),1000),tmp_path/'attempt',
                                 net.SearchProtocol(connection_protocol='irc'))
    assert [e['irc_converged'] for e in result['endpoints']]==[True,False]
    assert result['endpoints'][1]['irc_failure']=='inner_loop_convergence_failure'
    assert result['is_IRC'] is False and result['evidence']=='MLIP_IRC_not_converged_then_descent'
    assert result['status']=='validated_descents'
