"""GPA-style exploration: seeds -> index-one saddles -> observed minima -> graph.

Connections here are MLIP mode-displacement descents, explicitly not DFT IRC.
No expected product is used to accept or reject a physical connection.
"""
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import hashlib
import time
import traceback
import numpy as np
import networkx as nx
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
from ase.optimize import BFGS
from ase.mep import DimerControl, MinModeAtoms
from rdkit import Chem
from .physics import analyze_stationary
from .event_graph import geometry_mol, graph_smiles
from .search_seeds import ENCOUNTER_POLICIES, SEED_ABLATIONS, SEED_FEATURES, make_seed
from .exploration_actions import choose_action
from .saddle_optimization import StationaryDimerTranslate
from .event_classification import classify_event


class BudgetExceeded(RuntimeError):
    pass


ENDPOINT_ACCEPTANCE = ('full_system', 'carbon_skeleton')
TS_OPTIMIZERS = ('dimer', 'dimer+sella')
CONNECTION_PROTOCOLS = ('mode_displacement', 'irc')
FRONTIER_CONNECTIVITY = ('physical', 'species')


def join_same_species(graph, nodes):
    """Link minima that share a graph (conformers, encounter-complex orientations)."""
    first = {}
    for node in nodes:
        anchor = first.setdefault(node['graph_smiles'], node['id'])
        if anchor != node['id']:
            graph.add_edge(anchor, node['id'])
    return graph


def is_recoverable_failure(exc):
    """True for numerical or chemical-perception failures of one attempt.

    Examples are RDKit valence/charge perception errors, nonfinite model output,
    arithmetic errors and singular linear algebra. Programming and infrastructure
    errors (TypeError, KeyError, OSError, AssertionError, NotImplementedError,
    broken process pools, CUDA memory/device errors) must stop the run instead of
    being counted as failed chemistry.
    """
    if isinstance(exc, (NotImplementedError, RecursionError, BrokenProcessPool)):
        return False
    message = str(exc)
    if (type(exc).__name__ == 'OutOfMemoryError' or 'CUDA out of memory' in message
            or 'CUDA error' in message):
        return False
    return isinstance(exc, (ValueError, RuntimeError, ArithmeticError))


class CountedCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, backend, total_limit):
        super().__init__()
        self.backend = backend
        self.calls = 0
        self.model_calls = 0
        self.model_seconds = 0.
        self.total_limit = total_limit
        self.attempt_limit = total_limit

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if self.calls >= min(self.total_limit, self.attempt_limit):
            raise BudgetExceeded('Energy/force evaluation budget exhausted')
        self.calls += 1
        self.model_calls += 1
        started=time.perf_counter()
        try:
            self.backend.calculate(self.atoms, ['energy', 'forces'], all_changes)
        finally:
            self.model_seconds+=time.perf_counter()-started
        energy = float(self.backend.results['energy'])
        forces = np.asarray(self.backend.results['forces']).copy()
        if not np.isfinite(energy) or not np.isfinite(forces).all():
            raise ValueError('Nonfinite model prediction')
        self.results = {'energy': energy, 'forces': forces}

    def batch_forces(self,numbers,positions):
        """Budgets count geometries, not batched calls, preserving fair comparisons."""
        count=len(positions)
        if self.calls+count>min(self.total_limit,self.attempt_limit):
            raise BudgetExceeded('Insufficient geometry budget for Hessian batch')
        self.calls+=count;self.model_calls+=1
        started=time.perf_counter()
        try:
            result=self.backend.evaluate_many(numbers,positions)
        finally:
            self.model_seconds+=time.perf_counter()-started
        return result['forces']


@dataclass(frozen=True)
class SearchProtocol:
    seed_policy: str = 'reactant_graph_and_rigid_encounter_seeds_v6'
    max_attempts: int = 12
    seeds_per_node: int = 1
    symbolic_seed_scale: float = 1.0
    geometry_seeds_per_node: int = 9
    total_evaluations: int = 6000
    evaluations_per_attempt: int = 700
    ts_steps: int = 160
    initial_steps: int = 250
    descent_steps: int = 250
    initial_fmax: float = .003
    initial_curvature_steps: int = 3
    fmax: float = .005
    endpoint_acceptance: str = 'full_system'
    endpoint_core_fmax: float = .02
    hessian_step: float = .005
    hessian_batch_size: int = 1
    dimer_extrapolate_forces: bool = False
    mode_displacement: float = .15
    geometry_tolerance_A: float = .15
    energy_tolerance_eV: float = .03
    random_seed: int = 17
    # Each endpoint must satisfy E_TS - E_end >= min_barrier_eV. The legacy
    # value -1e-4 lets an endpoint lie up to 0.1 meV above the TS.
    min_barrier_eV: float = -1e-4
    # 'symbolic_only' (legacy): only bond_edits/arrows get the rigid encounter
    # search. 'matched_controls': geometry/center_random get the same search
    # with a proposal-free contact pair.
    encounter_policy: str = 'symbolic_only'
    # 'dimer' (legacy): Dimer to fmax. 'dimer+sella': Dimer to dimer_fmax, then
    # Sella P-RFO (order=1, Cartesian, as in DFT verification) to fmax.
    ts_optimizer: str = 'dimer'
    dimer_fmax: float = .05
    sella_steps: int = 100
    # 'mode_displacement' (legacy): +/- mode_displacement along the TS mode, then
    # BFGS. 'irc': Sella IRC on the MLIP in both directions, then BFGS to fmax.
    connection_protocol: str = 'mode_displacement'
    irc_dx: float = .08
    irc_steps: int = 300
    irc_fmax: float = .05
    irc_inner_fmax: float = .02
    # 'legacy': historical seeds. 'arrow_features_v1': arrow-only terms for
    # arrows, matched random-order progress axis for bond_edits (arrow_features).
    seed_features: str = 'legacy'
    # >0: symbolic proposals also come from up to N resonance forms of each node
    # (symbolic_library.ResonanceAwareLibrary); 0 is the historical behaviour.
    proposal_resonance_forms: int = 0
    # Iteration cap of the symbolic seed-geometry least-squares fit; 200 is the
    # historical value. A fit that converges earlier is unchanged by a larger cap.
    seed_fit_max_nfev: int = 200
    # Fig. 5c arms (search_seeds.SEED_ABLATIONS); 'none' keeps every arrow feature.
    seed_feature_ablation: str = 'none'
    # 'physical' (legacy): expand only minima joined to the root by observed edges.
    # 'species': minima with the same graph also count as joined, so a reaction
    # from a re-oriented encounter complex of the root species is root-connected.
    frontier_connectivity: str = 'physical'

    def __post_init__(self):
        if self.encounter_policy not in ENCOUNTER_POLICIES:
            raise ValueError(f'Unknown encounter_policy: {self.encounter_policy}')
        if self.endpoint_acceptance not in ENDPOINT_ACCEPTANCE:
            raise ValueError(f'Unknown endpoint_acceptance: {self.endpoint_acceptance}')
        if self.ts_optimizer not in TS_OPTIMIZERS:
            raise ValueError(f'Unknown ts_optimizer: {self.ts_optimizer}')
        if self.connection_protocol not in CONNECTION_PROTOCOLS:
            raise ValueError(f'Unknown connection_protocol: {self.connection_protocol}')
        if self.seed_features not in SEED_FEATURES:
            raise ValueError(f'Unknown seed_features: {self.seed_features}')
        if self.proposal_resonance_forms < 0:
            raise ValueError('proposal_resonance_forms must be >= 0')
        if self.seed_fit_max_nfev < 1:
            raise ValueError('seed_fit_max_nfev must be >= 1')
        if self.seed_feature_ablation not in SEED_ABLATIONS:
            raise ValueError(f'Unknown seed_feature_ablation: {self.seed_feature_ablation}')
        if self.frontier_connectivity not in FRONTIER_CONNECTIVITY:
            raise ValueError(f'Unknown frontier_connectivity: {self.frontier_connectivity}')
        if self.seed_feature_ablation != 'none' and self.seed_features != 'arrow_features_v1':
            raise ValueError("seed_feature_ablation needs seed_features='arrow_features_v1'")


def aligned_rmsd(x, y):
    x, y = x-x.mean(0), y-y.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
    return float(np.sqrt(np.mean(np.sum((x @ rotation-y)**2, axis=1))))


def molecular_rmsd(mol_a, x_a, mol_b, x_b):
    if graph_smiles(mol_a) != graph_smiles(mol_b):
        return float('inf')
    matches = mol_a.GetSubstructMatches(mol_b, uniquify=False, useChirality=True, maxMatches=256)
    return min((aligned_rmsd(x_a[list(m)], x_b) for m in matches), default=float('inf'))


def atomic_json(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    tmp.replace(path)


def inspect_point(atoms, protocol):
    force_norms=np.linalg.norm(atoms.get_forces(),axis=1)
    batch=atoms.calc.batch_forces if protocol.hessian_batch_size>1 else None
    analysis = analyze_stationary(atoms, protocol.fmax, protocol.hessian_step,
                                  batch_forces=batch,batch_size=protocol.hessian_batch_size)
    analysis['force_norms_eV_A']=force_norms.tolist()
    modes = analysis.pop('modes')
    return analysis, modes


def annotate_carbon_skeleton(analysis,mol,protocol):
    """Record a relaxed endpoint criterion without hiding full-system checks."""
    fragments=Chem.GetMolFrags(mol)
    core=sorted(i for fragment in fragments
        if any(mol.GetAtomWithIdx(j).GetAtomicNum()==6 for j in fragment) for i in fragment)
    if not core:
        core=list(range(mol.GetNumAtoms()))
    spectator=sorted(set(range(mol.GetNumAtoms()))-set(core))
    force=np.asarray(analysis['force_norms_eV_A'])
    analysis.update(carbon_skeleton_atoms=core,spectator_atoms=spectator,
        carbon_skeleton_force_max_eV_A=float(force[core].max()),
        carbon_skeleton_force_converged=bool(force[core].max()<=protocol.endpoint_core_fmax),
        spectator_force_max_eV_A=(float(force[spectator].max()) if spectator else 0.),
        full_system_minimum=bool(analysis['force_converged'] and analysis['imaginary_count']==0))
    return analysis


def edge_chemistry(attempt_mols, node_mols):
    """Edge chemistry in the attempt's atom numbering (both endpoints share the TS's).

    A registered node may match an endpoint only up to an atom permutation, so
    its stored molecule can give wrong bond indices. The old node-frame result is
    kept under node_frame_endpoint_chemistry when it differs, for comparison
    with runs made before this change.
    """
    try:
        chemistry = classify_event(*attempt_mols)
    except Exception as exc:
        if not is_recoverable_failure(exc):
            raise
        # RDKit could not perceive or canonicalize these endpoints; keep the edge,
        # mark it unclassified, and fall back to graph identity for its kind.
        same = graph_smiles(attempt_mols[0]) == graph_smiles(attempt_mols[1])
        chemistry = dict(classification='unclassified', error=f'{type(exc).__name__}: {exc}',
                         cross_fragment_bonds=[], resonance_equivalent=same)
    fields = dict(endpoint_chemistry=chemistry, endpoint_chemistry_atom_frame='attempt')
    try:
        legacy = classify_event(*node_mols)
    except Exception as exc:
        if not is_recoverable_failure(exc):
            raise
        return fields
    if ((legacy['classification'], legacy['cross_fragment_bonds']) !=
            (chemistry['classification'], chemistry['cross_fragment_bonds'])):
        fields['node_frame_endpoint_chemistry'] = legacy
    return fields


def connection_status(endpoints, protocol):
    """Accept or reject two descended endpoints; no expected product is used."""
    def downhill(e):
        return e['barrier_eV'] >= protocol.min_barrier_eV
    if all(e['full_system_minimum'] and downhill(e) for e in endpoints):
        return 'validated_descents'
    if (protocol.endpoint_acceptance == 'carbon_skeleton' and
            all(e['carbon_skeleton_force_converged'] and downhill(e) for e in endpoints)):
        return 'validated_core_descents'
    return 'unresolved_minimum'


def refine_with_sella(seed, outdir, protocol, result, initial_calls):
    """Sella P-RFO refinement of a Dimer point that passed the loose dimer_fmax gate.

    Every Sella evaluation, including its iterative Hessian estimates, goes
    through seed.calc and is billed like any other geometry.
    """
    calls = seed.calc.calls
    dimer_force = float(np.linalg.norm(seed.get_forces(), axis=1).max())
    stage = dict(ts_optimizer=protocol.ts_optimizer, dimer_fmax_eV_A=protocol.dimer_fmax,
                 dimer_converged=result['optimizer_converged'], dimer_force_max_eV_A=dimer_force,
                 dimer_evaluations=calls-initial_calls, sella_run=False)
    result['ts_optimization'] = stage
    if dimer_force > protocol.dimer_fmax:
        return
    from sella import Sella
    with Sella(seed, order=1, internal=False, logfile=str(outdir/'sella.log'),
               trajectory=str(outdir/'sella.traj')) as opt:
        stage['sella_converged'] = bool(opt.run(fmax=protocol.fmax, steps=protocol.sella_steps))
    stage.update(sella_run=True, sella_evaluations=seed.calc.calls-calls)
    result['optimizer_converged'] = stage['sella_converged']


def integrate_irc(atoms, outdir, protocol, direction):
    """One Sella IRC branch on the MLIP (DFT-verification settings); BFGS polish follows.

    An inner-loop convergence failure stops the branch at its last IRC point and
    is recorded as an unconverged IRC; the caller still polishes that point, so
    the endpoint is a descended minimum but the edge is not IRC evidence.
    """
    from sella import IRC
    from sella.optimize.irc import IRCInnerLoopConvergenceFailure
    calls = atoms.calc.calls
    failure = None
    with IRC(atoms, dx=protocol.irc_dx, ninner_iter=20, keep_going=False,
             logfile=str(outdir/f'irc_{direction}.log'),
             trajectory=str(outdir/f'irc_{direction}.traj')) as irc:
        try:
            converged = bool(irc.run(fmax=protocol.irc_fmax, fmax_inner=protocol.irc_inner_fmax,
                                     steps=protocol.irc_steps, direction=direction))
        except IRCInnerLoopConvergenceFailure:
            converged, failure = False, 'inner_loop_convergence_failure'
    return dict(irc_direction=direction, irc_converged=converged, irc_failure=failure,
                irc_evaluations=atoms.calc.calls-calls,
                irc_end_energy_eV=float(atoms.get_potential_energy()))


def descend_endpoints(seed, modes, ts_energy, calculator, outdir, protocol, charge, result):
    """Both endpoints of an index-one saddle: +/- mode displacement or MLIP IRC, then BFGS.

    Appends the analysed endpoints to result['endpoints'] and sets the IRC
    evidence fields, endpoint_acceptance and status. Shared by the search and by
    external baselines so that every arm uses the same post-processing chain.
    """
    np.save(outdir/'negative_mode.npy', modes[0])
    irc = protocol.connection_protocol == 'irc'
    for sign in (-1, 1):
        end = seed.copy()
        if irc:
            end.calc = calculator
            irc_info = integrate_irc(end, outdir, protocol, 'reverse' if sign < 0 else 'forward')
        else:
            end.positions += sign * protocol.mode_displacement * modes[0]
            end.calc = calculator
        with BFGS(end, maxstep=.1, logfile=str(outdir/f'descent_{sign}.log'),
                  trajectory=str(outdir/f'descent_{sign}.traj')) as opt:
            opt.run(fmax=protocol.fmax, steps=protocol.descent_steps)
        analysis, _ = inspect_point(end, protocol)
        if irc:
            irc_info['irc_polish_energy_change_eV'] = analysis['energy_eV']-irc_info['irc_end_energy_eV']
            analysis.update(irc_info)
        analysis['positions_A'] = end.positions.tolist()
        analysis['barrier_eV'] = ts_energy - analysis['energy_eV']
        mol = geometry_mol(end.numbers, end.positions, charge)
        annotate_carbon_skeleton(analysis,mol,protocol)
        analysis['graph_smiles'] = graph_smiles(mol)
        result['endpoints'].append(analysis)
        write(outdir/f'minimum_{sign}.xyz', end, write_results=False)
    if irc:
        converged = all(e['irc_converged'] for e in result['endpoints'])
        result.update(is_IRC=converged, evidence=('MLIP_bidirectional_IRC' if converged else
                                                  'MLIP_IRC_not_converged_then_descent'))
    result['endpoint_acceptance']=protocol.endpoint_acceptance
    result['status'] = connection_status(result['endpoints'], protocol)


def search_connection(seed, direction, calculator, outdir, protocol, charge=0):
    outdir.mkdir(parents=True, exist_ok=False)
    started, initial_calls = time.time(), calculator.calls
    initial_model_calls,initial_model_seconds=calculator.model_calls,calculator.model_seconds
    result = dict(status='started', evidence='MLIP_two_sided_mode_displacement_descent',
                  is_IRC=False, DFT_verified=False, endpoints=[],
                  optimizer_force_criterion='original_cartesian_per_atom_norm')
    try:
        seed.calc = calculator
        write(outdir / 'seed.xyz', seed, write_results=False)
        np.save(outdir / 'seed_positions_A.npy', seed.positions)   # full precision for exact replays
        np.save(outdir / 'seed_direction.npy', direction)
        with DimerControl(logfile=str(outdir/'dimer.log'), dimer_separation=.005,
                          maximum_translation=.1, max_num_rot=3,
                          extrapolate_forces=protocol.dimer_extrapolate_forces,
                          f_rot_min=.01, f_rot_max=.1) as control:
            mm = MinModeAtoms(seed, control=control, eigenmodes=[direction.copy()],
                              random_seed=protocol.random_seed)
            dimer_target = protocol.fmax if protocol.ts_optimizer == 'dimer' else protocol.dimer_fmax
            with StationaryDimerTranslate(mm, logfile=str(outdir/'opt.log'),
                                  trajectory=str(outdir/'search.traj')) as opt:
                result['optimizer_converged'] = bool(opt.run(fmax=dimer_target, steps=protocol.ts_steps))
        if protocol.ts_optimizer == 'dimer+sella':
            refine_with_sella(seed, outdir, protocol, result, initial_calls)
        residual_force=float(np.linalg.norm(seed.get_forces(),axis=1).max())
        if residual_force>protocol.fmax:
            result['status']='ts_force_unconverged'
            result['ts']=dict(force_max_eV_A=residual_force,force_converged=False)
            write(outdir/'ts.xyz',seed,write_results=False)
            return result
        ts, modes = inspect_point(seed, protocol)
        result['ts'] = ts
        write(outdir/'ts.xyz', seed, write_results=False)
        if not ts['force_converged'] or ts['imaginary_count'] != 1:
            result['status'] = 'not_index_one'
            return result
        descend_endpoints(seed, modes, ts['energy_eV'], calculator, outdir, protocol, charge, result)
        if result['status'] != 'unresolved_minimum':
            result['ts_positions_A'] = seed.positions.tolist()
    except BudgetExceeded:
        result['status'] = 'budget_exhausted'
    except Exception as exc:
        recoverable = is_recoverable_failure(exc)
        # Unrecoverable: a bug or an infrastructure fault. result.json is still
        # written below, then the run stops.
        result['status'] = 'calculation_failed' if recoverable else 'fatal_error'
        result['error'] = f'{type(exc).__name__}: {exc}'
        (outdir/'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        if not recoverable:
            raise
    finally:
        result['evaluations'] = calculator.calls - initial_calls
        result['model_calls'] = calculator.model_calls-initial_model_calls
        result['model_seconds'] = calculator.model_seconds-initial_model_seconds
        result['seconds'] = time.time() - started
        atomic_json(outdir/'result.json', result)
    return result


def initialize_root(start, calculator, outdir, protocol):
    """Relax and curvature-check the only user-supplied starting geometry."""
    outdir=Path(outdir)
    root=Atoms(numbers=start['atomic_numbers'],positions=start['positions_A'])
    root.calc=calculator
    with BFGS(root,maxstep=.1,logfile=str(outdir/'initial.log'),
              trajectory=str(outdir/'initial.traj')) as opt:
        opt.run(fmax=protocol.initial_fmax,steps=protocol.initial_steps)
    root_info,root_modes=inspect_point(root,protocol)
    checks=[root_info.copy()]
    for repair in range(protocol.initial_curvature_steps):
        if not root_info['force_converged'] or root_info['imaginary_count']==0:
            break
        root.positions+=protocol.mode_displacement*root_modes[0]
        with BFGS(root,maxstep=.1,logfile=str(outdir/f'initial_curvature_{repair}.log'),
                  trajectory=str(outdir/f'initial_curvature_{repair}.traj')) as opt:
            opt.run(fmax=protocol.initial_fmax,steps=protocol.initial_steps)
        root_info,root_modes=inspect_point(root,protocol)
        checks.append(root_info.copy())
    write(outdir/'initial.xyz',root,write_results=False)
    root_mol=geometry_mol(root.numbers,root.positions,start['charge'])
    annotate_carbon_skeleton(root_info,root_mol,protocol)
    root_info.update(positions_A=root.positions.tolist(),graph_smiles=graph_smiles(root_mol))
    return root_info,checks,calculator.calls


def explore(start, library, backend, strategy, outdir, protocol=SearchProtocol()):
    """Only an observed initial geometry and a reusable symbol library enter search."""
    outdir = Path(outdir)
    if strategy not in ('geometry','center_random','bond_edits','arrows','hybrid'):
        raise ValueError('Unknown search strategy')
    outdir.mkdir(parents=True, exist_ok=False)
    calculator = CountedCalculator(backend, protocol.total_evaluations)
    numbers = start['atomic_numbers']
    validator=getattr(backend,'validate_system',None)
    if validator is not None:
        validator(numbers,start['charge'],start['multiplicity'])
    if start['multiplicity'] != 1:
        raise ValueError('Current Lewis-graph registry supports closed-shell singlet exploration only')
    if strategy!='geometry' and not set(numbers)<={1,6,7,8,15}:
        raise ValueError('Current symbolic proposal layer only supports CHNOP; use geometry or add reviewed actions')
    nodes, edges, attempts, mols = [], [], [], []
    report = dict(start=start, strategy=strategy, protocol=asdict(protocol), nodes=nodes,
                  edges=edges, attempts=attempts, evidence='MLIP_descents_not_DFT_IRC',
                  library_policy=library.policy if library is not None else 'geometry_only',
                  status='running', unsupported_nodes=[],exhausted_nodes=[])
    started=time.perf_counter()
    def save():
        report['evaluations'] = calculator.calls
        report['model_calls'] = calculator.model_calls
        report['model_seconds'] = calculator.model_seconds
        report['elapsed_seconds'] = time.perf_counter()-started
        atomic_json(outdir/'network.json', report)
    def register(end, depth):
        """Return (node id, endpoint graph in this attempt's atom numbering)."""
        x = np.asarray(end['positions_A'])
        m = geometry_mol(numbers, x, start['charge'])
        for node, oldmol in zip(nodes, mols):
            if (abs(node['energy_eV']-end['energy_eV']) <= protocol.energy_tolerance_eV and
                molecular_rmsd(oldmol, np.array(node['positions_A']), m, x) < protocol.geometry_tolerance_A):
                return node['id'], m
        idx = len(nodes)
        nodes.append(dict(id=idx, depth_discovered=depth, **end))
        mols.append(m)
        return idx, m
    try:
        root_info,checks,initialization_evaluations=initialize_root(start,calculator,outdir,protocol)
        report['initial_relaxation_checks']=checks
        report['initial_point']=root_info.copy()
        report['initialization_evaluations']=initialization_evaluations
        root_valid=(root_info['full_system_minimum'] if protocol.endpoint_acceptance=='full_system'
                    else root_info['carbon_skeleton_force_converged'])
        if not root_valid:
            report['status'] = 'initial_minimum_unresolved'
            return report
        register(root_info, 0)
        expanded = set()
        visits = {}
        actions_used={}
        proposal_cache={}
        while len(attempts) < protocol.max_attempts and calculator.calls < protocol.total_evaluations:
            graph = nx.Graph()
            graph.add_nodes_from(range(len(nodes)))
            graph.add_edges_from(e['nodes'] for e in edges)
            if protocol.frontier_connectivity == 'species':
                join_same_species(graph, nodes)
            connected = nx.node_connected_component(graph, 0)
            available = sorted(connected - set(report['unsupported_nodes']) - set(report['exhausted_nodes']))
            if not available:
                break
            # Prefer unvisited basins, then spend remaining budget on new seeds
            # at covered basins. A failed first batch is not a proof of exhaustion.
            species_visits={}
            for i,count in visits.items():
                smi=nodes[i]['graph_smiles'];species_visits[smi]=species_visits.get(smi,0)+count
            node_id = min(available, key=lambda i: (species_visits.get(nodes[i]['graph_smiles'],0),visits.get(i,0),i))
            visit = visits.get(node_id, 0)
            visits[node_id] = visit + 1
            expanded.add(node_id)
            node = nodes[node_id]
            state = Atoms(numbers=numbers, positions=node['positions_A'])
            if node_id not in proposal_cache:
                try:
                    proposal_cache[node_id]=library.propose(mols[node_id]) if strategy!='geometry' else []
                except Exception as exc:
                    if not is_recoverable_failure(exc):
                        raise
                    # As in the parallel scheduler: no proposals for this node, recorded.
                    proposal_cache[node_id]=[]
                    report.setdefault('proposal_failures',[]).append(
                        dict(node=node_id,error=f'{type(exc).__name__}: {exc}'))
            proposals=proposal_cache[node_id]
            if strategy not in ('geometry','hybrid') and not proposals:
                report['unsupported_nodes'].append(node_id)
                continue
            for sample in range(protocol.seeds_per_node):
                if len(attempts) >= protocol.max_attempts or calculator.calls >= protocol.total_evaluations:
                    break
                aid = len(attempts)
                used=actions_used.setdefault(node_id,{})
                choice=choose_action(proposals,used,{n['graph_smiles'] for n in nodes},strategy,
                                     protocol.geometry_seeds_per_node,3*protocol.seeds_per_node)
                if choice is None:
                    report['exhausted_nodes'].append(node_id)
                    break
                proposal,key,variant,seed_strategy=choice
                used[key]=used.get(key,0)+1
                # Seed randomness follows state geometry + action + visit, not
                # discovery-order node IDs. Shared symbolic controls get the same RNG.
                identity = json.dumps(dict(graph=node['graph_smiles'],
                    positions=np.round(state.positions,6).tolist(),
                    edits=proposal['edits'] if proposal else None,
                    seed=protocol.random_seed,variant=variant),sort_keys=True)
                seed_rng = int.from_bytes(hashlib.sha256(identity.encode()).digest()[:4],'little')
                try:
                    seed_started=time.perf_counter()
                    x, direction, meta = make_seed(state, mols[node_id], seed_strategy, proposal,
                                                   variant, seed_rng, protocol.symbolic_seed_scale,
                                                   protocol.encounter_policy, protocol.seed_features,
                                                   protocol.seed_fit_max_nfev, protocol.seed_feature_ablation)
                except Exception as exc:
                    if not is_recoverable_failure(exc):
                        raise
                    dest = outdir/f'attempt_{aid:03d}'
                    dest.mkdir()
                    failure = dict(status='seed_generation_failed',error=str(exc),evaluations=0,
                                   is_IRC=False,DFT_verified=False,endpoints=[])
                    atomic_json(dest/'result.json', failure)
                    attempts.append(dict(id=aid,source_node=node_id,proposal=proposal,
                        status=failure['status'],evaluations=0,artifact=f'attempt_{aid:03d}/result.json'))
                    save()
                    continue
                meta['random_seed'] = seed_rng
                meta['seed_strategy']=seed_strategy
                meta['generation_seconds']=time.perf_counter()-seed_started
                trial = state.copy()
                trial.positions = x
                calculator.attempt_limit = min(calculator.total_limit,
                    calculator.calls + protocol.evaluations_per_attempt)
                dest = outdir/f'attempt_{aid:03d}'
                calls_before = calculator.calls
                try:
                    result = search_connection(trial,direction,calculator,dest,protocol,start['charge'])
                except Exception as exc:
                    # Keep the accounting identity evaluations = init + sum(attempts).
                    attempts.append(dict(id=aid, source_node=node_id, proposal=meta,
                        status='fatal_error', error=f'{type(exc).__name__}: {exc}',
                        evaluations=calculator.calls-calls_before,
                        artifact=f'attempt_{aid:03d}/result.json'))
                    raise
                attempt = dict(id=aid, source_node=node_id, proposal=meta,
                               status=result['status'], evaluations=result['evaluations'],
                               seconds=result['seconds'],
                               artifact=f'attempt_{aid:03d}/result.json')
                attempts.append(attempt)
                if result['status'] in ('validated_descents','validated_core_descents'):
                    registered = [register(e, node['depth_discovered']+1) for e in result['endpoints']]
                    ends = [item[0] for item in registered]
                    attempt['observed_nodes'] = ends
                    attempt['source_is_endpoint'] = node_id in ends
                    if ends[0] == ends[1]:
                        attempt['status'] = 'same_basin_return'
                    else:
                        duplicate = False
                        for edge in edges:
                            if sorted(edge['nodes']) != sorted(ends):
                                continue
                            # Atom identity is fixed throughout a run. Conservative TS
                            # comparison keeps uncertain symmetry duplicates separate.
                            if (abs(edge['ts_energy_eV']-result['ts']['energy_eV']) < .03 and
                                aligned_rmsd(np.array(edge['ts_positions_A']),
                                             np.array(result['ts_positions_A'])) < .15):
                                duplicate = True
                                break
                        attempt['status'] = 'duplicate_connection' if duplicate else 'new_connection'
                        if not duplicate:
                            chemistry=edge_chemistry([item[1] for item in registered],
                                                     [mols[i] for i in ends])
                            edges.append(dict(id=len(edges), nodes=ends, attempt=aid,
                                proposed_from=node_id, ts_energy_eV=result['ts']['energy_eV'],
                                ts_positions_A=result['ts_positions_A'],
                                barriers_eV=[e['barrier_eV'] for e in result['endpoints']],
                                source_connected=node_id in ends,
                                endpoint_acceptance=result['status'],
                                **chemistry,
                                kind=('conformational' if chemistry['endpoint_chemistry']['resonance_equivalent']
                                      else 'chemical')))
                save()
                print(json.dumps(dict(start=start['id'], strategy=strategy, attempt=aid,
                                      status=attempt['status'], evaluations=calculator.calls)), flush=True)
        report['expanded_nodes'] = sorted(expanded)
        report['searched_nodes'] = sorted({a['source_node'] for a in attempts})
        report['node_batches'] = visits
        report['status'] = 'completed'
        report['stop_reason'] = ('evaluation_budget' if calculator.calls >= protocol.total_evaluations else
                                 'attempt_budget' if len(attempts) >= protocol.max_attempts else 'frontier_exhausted')
        graph = nx.Graph()
        graph.add_nodes_from(range(len(nodes)))
        graph.add_edges_from(e['nodes'] for e in edges)
        report['root_component_nodes'] = sorted(nx.node_connected_component(graph, 0))
        if protocol.frontier_connectivity == 'species':
            report['species_root_component_nodes'] = sorted(
                nx.node_connected_component(join_same_species(graph, nodes), 0))
    except BudgetExceeded:
        report['status'] = 'initialization_budget_exhausted'
    except Exception as exc:
        report['status'] = 'aborted_error'
        report['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        save()
    return report
