"""GPA-style exploration: seeds -> index-one saddles -> observed minima -> graph.

Connections here are MLIP mode-displacement descents, explicitly not DFT IRC.
No expected product is used to accept or reject a physical connection.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import hashlib
import os
import time
import traceback
import numpy as np
import networkx as nx
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
from ase.optimize import BFGS, FIRE
from ase.mep import DimerControl, MinModeAtoms, NEB
from .physics import analyze_stationary, finite_hessian
from .event_graph import geometry_mol, graph_smiles
from .search_seeds import make_seed, internal_direction
from .symbolic_endpoint import product_geometry
from .exploration_actions import choose_action
from .saddle_optimization import StationaryDimerTranslate
from .event_classification import classify_event


class BudgetExceeded(RuntimeError):
    pass


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
    geometry_seeds_per_node: int = 9
    total_evaluations: int = 6000
    evaluations_per_attempt: int = 700
    ts_steps: int = 160
    sella_steps: int = 400
    descent_steps: int = 250
    endpoint_steps: int = 300
    neb_steps: int = 250
    neb_images: int = 7
    neb_fmax: float = .08
    initial_fmax: float = .003
    initial_curvature_steps: int = 3
    fmax: float = .005
    hessian_step: float = .005
    hessian_batch_size: int = 1
    dimer_extrapolate_forces: bool = False
    mode_displacement: float = .15
    geometry_tolerance_A: float = .15
    energy_tolerance_eV: float = .03
    random_seed: int = 17


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
    batch=atoms.calc.batch_forces if protocol.hessian_batch_size>1 else None
    analysis = analyze_stationary(atoms, protocol.fmax, protocol.hessian_step,
                                  batch_forces=batch,batch_size=protocol.hessian_batch_size)
    modes = analysis.pop('modes')
    return analysis, modes


def validate_saddle_descents(seed, calculator, outdir, protocol, charge, result):
    """Apply the common index-one and two-minimum acceptance contract."""
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
    np.save(outdir/'negative_mode.npy', modes[0])
    for sign in (-1, 1):
        end = seed.copy()
        end.positions += sign * protocol.mode_displacement * modes[0]
        end.calc = calculator
        with BFGS(end, maxstep=.1, logfile=str(outdir/f'descent_{sign}.log'),
                  trajectory=str(outdir/f'descent_{sign}.traj')) as opt:
            opt.run(fmax=protocol.fmax, steps=protocol.descent_steps)
        analysis, _ = inspect_point(end, protocol)
        analysis['positions_A'] = end.positions.tolist()
        analysis['barrier_eV'] = ts['energy_eV'] - analysis['energy_eV']
        mol = geometry_mol(end.numbers, end.positions, charge)
        analysis['graph_smiles'] = graph_smiles(mol)
        result['endpoints'].append(analysis)
        write(outdir/f'minimum_{sign}.xyz', end, write_results=False)
    if not all(e['force_converged'] and e['imaginary_count'] == 0 and e['barrier_eV'] >= -1e-4
               for e in result['endpoints']):
        result['status'] = 'unresolved_minimum'
    else:
        result['status'] = 'validated_descents'
        result['ts_positions_A'] = seed.positions.tolist()
    return result


def search_connection(seed, direction, calculator, outdir, protocol, charge=0):
    outdir.mkdir(parents=True, exist_ok=False)
    started, initial_calls = time.time(), calculator.calls
    result = dict(status='started', evidence='MLIP_two_sided_mode_displacement_descent',
                  is_IRC=False, DFT_verified=False, endpoints=[],
                  optimizer_force_criterion='original_cartesian_per_atom_norm')
    try:
        seed.calc = calculator
        write(outdir / 'seed.xyz', seed, write_results=False)
        np.save(outdir / 'seed_direction.npy', direction)
        with DimerControl(logfile=str(outdir/'dimer.log'), dimer_separation=.005,
                          maximum_translation=.1, max_num_rot=3,
                          extrapolate_forces=protocol.dimer_extrapolate_forces,
                          f_rot_min=.01, f_rot_max=.1) as control:
            mm = MinModeAtoms(seed, control=control, eigenmodes=[direction.copy()],
                              random_seed=protocol.random_seed)
            with StationaryDimerTranslate(mm, logfile=str(outdir/'opt.log'),
                                  trajectory=str(outdir/'search.traj')) as opt:
                result['optimizer_converged'] = bool(opt.run(fmax=protocol.fmax, steps=protocol.ts_steps))
        validate_saddle_descents(seed,calculator,outdir,protocol,charge,result)
    except BudgetExceeded:
        result['status'] = 'budget_exhausted'
    except Exception as exc:
        result['status'] = 'calculation_failed'
        result['error'] = f'{type(exc).__name__}: {exc}'
        (outdir/'error.log').write_text(traceback.format_exc(), encoding='utf-8')
    finally:
        result['evaluations'] = calculator.calls - initial_calls
        result['seconds'] = time.time() - started
        atomic_json(outdir/'result.json', result)
    return result


def search_sella_connection(seed, direction, calculator, outdir, protocol, charge=0):
    """Refine a path maximum with Sella P-RFO and apply the shared critic."""
    cache=outdir.parent/'sella_jax_cache';cache.mkdir(exist_ok=True)
    os.environ.setdefault('JAX_COMPILATION_CACHE_DIR',str(cache.resolve()))
    from sella import Sella
    outdir.mkdir(parents=True, exist_ok=False)
    started, initial_calls = time.time(), calculator.calls
    result = dict(status='started', evidence='Sella_PRFO_then_MLIP_two_sided_descents',
                  is_IRC=False, DFT_verified=False, endpoints=[])
    try:
        seed.calc=calculator
        write(outdir/'seed.xyz',seed,write_results=False)
        np.save(outdir/'seed_direction.npy',direction)
        def hessian_function(atoms):
            return finite_hessian(atoms,protocol.hessian_step,
                batch_forces=calculator.batch_forces,batch_size=protocol.hessian_batch_size)
        with Sella(seed,order=1,v0=direction.ravel(),internal=False,
                   hessian_function=hessian_function,diag_every_n=10,
                   logfile=str(outdir/'sella.log'),trajectory=str(outdir/'search.traj')) as opt:
            result['optimizer_converged']=bool(opt.run(fmax=protocol.fmax,steps=protocol.sella_steps))
        validate_saddle_descents(seed,calculator,outdir,protocol,charge,result)
    except BudgetExceeded:
        result['status']='budget_exhausted'
    except Exception as exc:
        result['status']='calculation_failed'
        result['error']=f'{type(exc).__name__}: {exc}'
        (outdir/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        result['evaluations']=calculator.calls-initial_calls
        result['seconds']=time.time()-started
        atomic_json(outdir/'result.json',result)
    return result


def search_neb_connection(source, source_mol, proposal, calculator, outdir, protocol,
                          random_seed, charge=0):
    """Construct a symbolic endpoint, optimize a NEB, then refine its saddle."""
    outdir.mkdir(parents=True, exist_ok=False)
    started, initial_calls = time.time(), calculator.calls
    result = dict(status='started', evidence='symbolic_endpoint_NEB_Sella_then_MLIP_descents',
                  is_IRC=False, DFT_verified=False, endpoints=[])
    try:
        guess, endpoint_meta = product_geometry(source.numbers, source.positions, source_mol,
                                                 proposal, random_seed)
        endpoint = Atoms(numbers=source.numbers, positions=guess, calculator=calculator)
        write(outdir/'proposed_endpoint.xyz', endpoint, write_results=False)
        with BFGS(endpoint, maxstep=.1, logfile=str(outdir/'endpoint.log'),
                  trajectory=str(outdir/'endpoint.traj')) as opt:
            endpoint_converged = bool(opt.run(fmax=protocol.initial_fmax,
                                               steps=protocol.endpoint_steps))
        endpoint_analysis, _ = inspect_point(endpoint, protocol)
        endpoint_graph = graph_smiles(geometry_mol(endpoint.numbers, endpoint.positions, charge))
        result['symbolic_endpoint'] = dict(**endpoint_meta, optimized_graph=endpoint_graph,
            optimizer_converged=endpoint_converged, stationary=endpoint_analysis,
            positions_A=endpoint.positions.tolist())
        write(outdir/'optimized_endpoint.xyz', endpoint, write_results=False)
        if (not endpoint_analysis['force_converged'] or endpoint_analysis['imaginary_count'] != 0
                or endpoint_graph != proposal['predicted_graph']):
            result['status'] = 'symbolic_endpoint_unresolved'
            return result
        images = [source.copy()]
        images += [source.copy() for _ in range(protocol.neb_images-2)]
        images += [endpoint.copy()]
        for image in images:
            image.calc = calculator
        neb = NEB(images, climb=False, allow_shared_calculator=True,
                  remove_rotation_and_translation=True, method='improvedtangent')
        neb.interpolate(method='idpp')
        with FIRE(neb, maxstep=.05, dt=.05, logfile=str(outdir/'neb.log'),
                  trajectory=str(outdir/'neb.traj')) as opt:
            neb_converged = bool(opt.run(fmax=protocol.neb_fmax, steps=protocol.neb_steps))
        energies = [float(image.get_potential_energy()) for image in images]
        maxima=[i for i in range(1,len(images)-1)
                if energies[i]>=energies[i-1] and energies[i]>=energies[i+1]]
        if not maxima:maxima=[1+int(np.argmax(energies[1:-1]))]
        maxima.sort(key=lambda i:energies[i],reverse=True)
        top=maxima[0]
        result['neb'] = dict(images=protocol.neb_images, optimizer_converged=neb_converged,
            fmax_eV_A=protocol.neb_fmax, energies_eV=energies, highest_image=top,
            local_maxima=maxima,
            endpoint_graphs=[graph_smiles(geometry_mol(i.numbers, i.positions, charge))
                             for i in (images[0], images[-1])])
        for index, image in enumerate(images):
            write(outdir/f'neb_image_{index:02d}.xyz', image, write_results=False)
        source_graph=graph_smiles(source_mol);target_graph=proposal['predicted_graph']
        result['refinements']=[];selected=None
        for rank,index in enumerate(maxima):
            direction = internal_direction(images[index].positions,
                images[index+1].positions-images[index-1].positions)
            refined = search_sella_connection(images[index].copy(),direction,calculator,
                outdir/f'saddle_refinement_{rank:02d}_image_{index:02d}',protocol,charge)
            observed=[e.get('graph_smiles') for e in refined.get('endpoints',[])]
            target_pair=(refined['status']=='validated_descents' and
                         sorted(observed)==sorted([source_graph,target_graph]))
            result['refinements'].append(dict(rank=rank,image=index,status=refined['status'],
                endpoint_graphs=observed,target_graph_pair=target_pair,
                evaluations=refined['evaluations']))
            if target_pair:
                selected=refined;break
        if selected is None:
            result['status']='path_not_recovered'
        else:
            result.update({k:v for k,v in selected.items() if k not in ('evaluations','seconds')})
    except BudgetExceeded:
        result['status'] = 'budget_exhausted'
    except Exception as exc:
        result['status'] = 'calculation_failed'
        result['error'] = f'{type(exc).__name__}: {exc}'
        (outdir/'error.log').write_text(traceback.format_exc(), encoding='utf-8')
    finally:
        result['evaluations'] = calculator.calls - initial_calls
        result['seconds'] = time.time() - started
        atomic_json(outdir/'result.json', result)
    return result


def explore(start, library, backend, strategy, outdir, protocol=SearchProtocol()):
    """Only an observed initial geometry and a reusable symbol library enter search."""
    outdir = Path(outdir)
    if strategy not in ('geometry','center_random','bond_edits','arrows','hybrid','neb_arrows'):
        raise ValueError('Unknown search strategy')
    outdir.mkdir(parents=True, exist_ok=False)
    calculator = CountedCalculator(backend, protocol.total_evaluations)
    numbers = start['atomic_numbers']
    validator=getattr(backend,'validate_system',None)
    if validator is not None:
        validator(numbers,start['charge'],start['multiplicity'])
    if start['multiplicity'] != 1:
        raise ValueError('Current Lewis-graph registry supports closed-shell singlet exploration only')
    if strategy!='geometry' and (start['charge']!=0 or not set(numbers)<={1,6,7,8}):
        raise ValueError('Current symbolic proposal layer only supports neutral CHNO; use geometry or add reviewed actions')
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
        x = np.asarray(end['positions_A'])
        m = geometry_mol(numbers, x, start['charge'])
        for node, oldmol in zip(nodes, mols):
            if (abs(node['energy_eV']-end['energy_eV']) <= protocol.energy_tolerance_eV and
                molecular_rmsd(oldmol, np.array(node['positions_A']), m, x) < protocol.geometry_tolerance_A):
                return node['id']
        idx = len(nodes)
        nodes.append(dict(id=idx, depth_discovered=depth, **end))
        mols.append(m)
        return idx
    try:
        root = Atoms(numbers=numbers, positions=start['positions_A'])
        root.calc = calculator
        with BFGS(root, maxstep=.1, logfile=str(outdir/'initial.log'),
                  trajectory=str(outdir/'initial.traj')) as opt:
            opt.run(fmax=protocol.initial_fmax, steps=protocol.descent_steps)
        root_info, root_modes = inspect_point(root, protocol)
        report['initial_relaxation_checks'] = [root_info.copy()]
        # Symmetric UFF structures can have zero force at a torsional maximum.
        # A minimum initializer must use curvature, not force convergence alone.
        for repair in range(protocol.initial_curvature_steps):
            if not root_info['force_converged'] or root_info['imaginary_count'] == 0:
                break
            root.positions += protocol.mode_displacement * root_modes[0]
            with BFGS(root,maxstep=.1,logfile=str(outdir/f'initial_curvature_{repair}.log'),
                      trajectory=str(outdir/f'initial_curvature_{repair}.traj')) as opt:
                opt.run(fmax=protocol.initial_fmax,steps=protocol.descent_steps)
            root_info, root_modes = inspect_point(root,protocol)
            report['initial_relaxation_checks'].append(root_info.copy())
        report['initial_point'] = root_info.copy()
        report['initialization_evaluations'] = calculator.calls
        write(outdir/'initial.xyz', root, write_results=False)
        if not root_info['force_converged'] or root_info['imaginary_count'] != 0:
            report['status'] = 'initial_minimum_unresolved'
            return report
        root_info.update(positions_A=root.positions.tolist(), graph_smiles=graph_smiles(
            geometry_mol(numbers, root.positions, start['charge'])))
        register(root_info, 0)
        expanded = set()
        visits = {}
        actions_used={}
        proposal_cache={}
        while len(attempts) < protocol.max_attempts and calculator.calls < protocol.total_evaluations:
            graph = nx.Graph()
            graph.add_nodes_from(range(len(nodes)))
            graph.add_edges_from(e['nodes'] for e in edges)
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
                proposal_cache[node_id]=library.propose(mols[node_id]) if strategy!='geometry' else []
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
                                     protocol.geometry_seeds_per_node)
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
                    if strategy == 'neb_arrows':
                        x, direction = None, None
                        meta = dict(template_id=proposal['template_id'], origin=proposal.get('origin'),
                            predicted_graph=proposal['predicted_graph'], arrows=proposal['arrows'],
                            net_edits=proposal['edits'], reference_product_geometry_used=False)
                    else:
                        x, direction, meta = make_seed(state, mols[node_id], seed_strategy, proposal,
                                                       variant, seed_rng)
                except (ValueError, RuntimeError) as exc:
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
                calculator.attempt_limit = min(calculator.total_limit,
                    calculator.calls + protocol.evaluations_per_attempt)
                dest = outdir/f'attempt_{aid:03d}'
                if strategy == 'neb_arrows':
                    result = search_neb_connection(state, mols[node_id], proposal, calculator, dest,
                                                   protocol, seed_rng, start['charge'])
                else:
                    trial = state.copy()
                    trial.positions = x
                    result = search_connection(trial,direction,calculator,dest,protocol,start['charge'])
                attempt = dict(id=aid, source_node=node_id, proposal=meta,
                               status=result['status'], evaluations=result['evaluations'],
                               seconds=result['seconds'],
                               artifact=f'attempt_{aid:03d}/result.json')
                attempts.append(attempt)
                if result['status'] == 'validated_descents':
                    matched=[]
                    for end in result['endpoints']:
                        x=np.asarray(end['positions_A']);m=geometry_mol(numbers,x,start['charge'])
                        matched.append(abs(node['energy_eV']-end['energy_eV']) <= protocol.energy_tolerance_eV and
                            molecular_rmsd(mols[node_id],np.asarray(node['positions_A']),m,x) <
                            protocol.geometry_tolerance_A)
                    if not any(matched):
                        attempt['status']='detached_connection'
                        attempt['source_is_endpoint']=False
                        save()
                        continue
                    ends = [node_id if matched[i] else register(e,node['depth_discovered']+1)
                            for i,e in enumerate(result['endpoints'])]
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
                            chemistry=classify_event(mols[ends[0]],mols[ends[1]])
                            edges.append(dict(id=len(edges), nodes=ends, attempt=aid,
                                proposed_from=node_id, ts_energy_eV=result['ts']['energy_eV'],
                                ts_positions_A=result['ts_positions_A'],
                                barriers_eV=[e['barrier_eV'] for e in result['endpoints']],
                                source_connected=node_id in ends,
                                endpoint_chemistry=chemistry,
                                kind='conformational' if chemistry['resonance_equivalent'] else 'chemical'))
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
    except BudgetExceeded:
        report['status'] = 'initialization_budget_exhausted'
    finally:
        save()
    return report
