"""Single-ended GSM baseline (pyGSM) with the same MLIP, budget and post-processing.

Arms (strategy names in network.json):
  pygsm_edits: driving coordinates are ADD/BREAK of each symbolic proposal's net
               bond edits, i.e. the same proposals the arrows/bond_edits arms use;
  pygsm_b2f2:  up to two breaks and two formations in a seeded random order with
               only a connectivity-degree ceiling and a 4 A formation-distance
               cut-off (no chemical prior).
Every evaluation counts against protocol.total_evaluations: the MLIP root
initialization, each pyGSM run (read from the budgeted calculator's counter
file, so aborted runs are billed too) and the shared post-processing (Sella
P-RFO from the GSM TS node, index-one check, endpoint descent). The output is a
network.json in the explore() layout, so score_reference_recovery.py applies.

Usage: run_pygsm_baseline.py --starts F --start-id ID --arm edits|b2f2 --seed S
         --protocol-json P --potential aimnet2-rxn --pygsm DIR --outdir DIR
"""
import argparse
from collections import Counter
import dataclasses
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

import numpy as np
from ase import Atoms
from ase.io import read, write
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts/exploration'))
from mechbridge.event_graph import bond_orders, geometry_mol  # noqa: E402
from mechbridge.potentials import load_potential  # noqa: E402
from mechbridge.provenance import runtime_provenance  # noqa: E402
from mechbridge.reaction_network import (BudgetExceeded, CountedCalculator, aligned_rmsd,  # noqa: E402
                                         atomic_json, descend_endpoints, edge_chemistry,
                                         initialize_root, inspect_point, is_recoverable_failure,
                                         molecular_rmsd)
from mechbridge.oracle_library import OracleReferenceLibrary  # noqa: E402
from mechbridge.symbolic_library import ArrowLibrary, ResonanceAwareLibrary  # noqa: E402
from run_network_exploration import protocol_from_json  # noqa: E402

MAX_DEGREE = {1: 1, 6: 4, 7: 3, 8: 2}    # neutral closed-shell connectivity ceiling (CHNO)
MIN_REMAINING = 300                       # do not start a GSM run with less budget left


def proposal_driving_sets(proposals):
    """Unique ADD/BREAK sets from symbolic proposals, in proposal order."""
    seen, sets = set(), []
    for proposal in proposals:
        coords = []
        for edit in proposal['edits']:
            i, j = sorted(int(a) for a in edit['atoms'])
            if edit['before'] == 0 and edit['after'] > 0:
                coords.append(('ADD', i, j))
            elif edit['before'] > 0 and edit['after'] == 0:
                coords.append(('BREAK', i, j))
        key = tuple(sorted(coords))
        if coords and key not in seen:
            seen.add(key)
            sets.append(list(key))
    return sets


def b2f2_driving_sets(numbers, positions, bonds, rng, count, max_distance=4.0, max_draws=200000):
    """Uniform random valid <=2-break/<=2-form sets (rejection sampling, seeded, no repeats).

    Full enumeration is infeasible for 40-atom encounters (~1e7 sets), and a
    search tries at most max_attempts sets, so sets are drawn uniformly from all
    candidates (class chosen in proportion to its size) and kept if every atom
    stays within the neutral closed-shell degree ceiling.
    """
    from math import comb
    n = len(numbers)
    existing = sorted(tuple(b) for b in bonds)
    bonded = set(existing)
    distance = np.linalg.norm(positions[:, None]-positions[None], axis=-1)
    absent = [(i, j) for i in range(n) for j in range(i+1, n) if (i, j) not in bonded
              and not (numbers[i] == 1 and numbers[j] == 1) and distance[i, j] <= max_distance]
    degree = Counter()
    for i, j in existing:
        degree[i] += 1
        degree[j] += 1
    classes = [(nb, nf) for nb in range(3) for nf in range(3) if nb+nf and nb <= len(existing) and nf <= len(absent)]
    weights = np.array([comb(len(existing), nb)*comb(len(absent), nf) for nb, nf in classes], dtype=float)
    weights /= weights.sum()
    sets, seen = [], set()
    for _ in range(max_draws):
        if len(sets) >= count:
            break
        nb, nf = classes[rng.choice(len(classes), p=weights)]
        broken = [existing[k] for k in sorted(rng.choice(len(existing), size=nb, replace=False))]
        formed = [absent[k] for k in sorted(rng.choice(len(absent), size=nf, replace=False))]
        key = (tuple(broken), tuple(formed))
        if key in seen:
            continue
        seen.add(key)
        change = Counter(degree)
        for i, j in broken:
            change[i] -= 1
            change[j] -= 1
        for i, j in formed:
            change[i] += 1
            change[j] += 1
        if all(change[a] <= MAX_DEGREE.get(int(numbers[a]), 4) for a in range(n)):
            sets.append([('BREAK', int(i), int(j)) for i, j in broken] + [('ADD', int(i), int(j)) for i, j in formed])
    return sets


def run_gsm(workdir, atoms, coords, limit, task, args, start):
    """One pyGSM SE-GSM run; returns (billed evaluations, exit code, TS-node path or None)."""
    workdir.mkdir(parents=True, exist_ok=False)
    write(workdir/'reactant.xyz', atoms)
    (workdir/'isomers.txt').write_text('NEW\n' + ''.join(f'{op} {i+1} {j+1}\n' for op, i, j in coords))
    counter = workdir/'evaluations.count'
    kwargs = dict(potential=args.potential, root=str(ROOT), limit=int(limit), counter_file=str(counter),
                  charge=start['charge'], multiplicity=start['multiplicity'], threads=args.threads)
    command = [sys.executable, str(args.pygsm/'bin/gsm'), '-xyzfile', 'reactant.xyz', '-mode', 'SE_GSM',
               '-isomers', 'isomers.txt', '-package', 'ase',
               '--ase-class', 'mechbridge.gsm_bridge.BudgetedPotential', '--ase-kwargs', json.dumps(kwargs),
               '-num_nodes', str(args.num_nodes), '-ID', str(task), '-reactant_geom_fixed',
               '-charge', str(start['charge'])]
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT/'src'), str(args.pygsm),
                                                       os.environ.get('PYTHONPATH', '')]))
    started = time.time()
    with (workdir/'gsm.log').open('w', encoding='utf-8') as log:
        try:
            code = subprocess.run(command, cwd=workdir, env=env, stdout=log, stderr=subprocess.STDOUT,
                                  timeout=args.gsm_timeout).returncode
        except subprocess.TimeoutExpired:
            code = 'timeout'
    evaluations = int(counter.read_text()) if counter.exists() and counter.read_text().strip() else 0
    ts_path = workdir/f'TSnode_{task}.xyz'
    return evaluations, code, (ts_path if ts_path.exists() else None), time.time()-started


def post_process(ts_path, numbers, calculator, protocol, outdir, charge):
    """Same chain as the search after a TS guess: Sella P-RFO, index-one check, endpoint descent."""
    result = dict(status='started', endpoints=[], is_IRC=False, DFT_verified=False,
                  evidence='MLIP_two_sided_mode_displacement_descent')
    calls = calculator.calls
    try:
        guess = read(ts_path)
        if list(guess.numbers) != list(numbers):
            raise ValueError('TS node atom order differs from the start')
        ts = Atoms(numbers=numbers, positions=guess.positions)
        ts.calc = calculator
        from sella import Sella
        with Sella(ts, order=1, internal=False, logfile=str(outdir/'prfo.log'),
                   trajectory=str(outdir/'prfo.traj')) as opt:
            result['prfo_converged'] = bool(opt.run(fmax=protocol.fmax, steps=protocol.sella_steps))
        residual = float(np.linalg.norm(ts.get_forces(), axis=1).max())
        if residual > protocol.fmax:
            result.update(status='ts_force_unconverged', ts=dict(force_max_eV_A=residual))
            return result
        analysis, modes = inspect_point(ts, protocol)
        result['ts'] = analysis
        write(outdir/'ts.xyz', ts, write_results=False)
        if not analysis['force_converged'] or analysis['imaginary_count'] != 1:
            result['status'] = 'not_index_one'
            return result
        descend_endpoints(ts, modes, analysis['energy_eV'], calculator, outdir, protocol, charge, result)
        if result['status'] != 'unresolved_minimum':
            result['ts_positions_A'] = ts.positions.tolist()
    except BudgetExceeded:
        result['status'] = 'budget_exhausted'
    except Exception as exc:
        if not is_recoverable_failure(exc):
            raise
        result.update(status='calculation_failed', error=f'{type(exc).__name__}: {exc}')
        (outdir/'error.log').write_text(traceback.format_exc(), encoding='utf-8')
    finally:
        result['evaluations'] = calculator.calls - calls
        atomic_json(outdir/'post_result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starts', type=Path, required=True)
    parser.add_argument('--start-id', required=True)
    parser.add_argument('--arm', choices=['edits', 'b2f2', 'oracle'], required=True)
    parser.add_argument('--oracle-references', type=Path,
                        help='References file for the oracle arm (reference net edits; knows the answer)')
    parser.add_argument('--seed', type=int, default=17)
    parser.add_argument('--protocol-json', type=Path, required=True)
    parser.add_argument('--potential', default='aimnet2-rxn')
    parser.add_argument('--pygsm', type=Path, required=True)
    parser.add_argument('--num-nodes', type=int, default=20)
    parser.add_argument('--gsm-attempt-cap', type=int, default=6000,
                        help='Largest share of the budget a single GSM run may use')
    parser.add_argument('--gsm-timeout', type=float, default=7200.)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    protocol = protocol_from_json(args.protocol_json, args.seed)
    start = next(r for r in map(json.loads, args.starts.read_text(encoding='utf-8').splitlines())
                 if r['id'] == args.start_id)
    strategy = f'pygsm_{args.arm}'
    outdir = args.outdir/start['id']/strategy
    outdir.mkdir(parents=True, exist_ok=False)
    backend, potential = load_potential(args.potential, ROOT, threads=args.threads)
    backend.validate_system(start['atomic_numbers'], start['charge'], start['multiplicity'])
    total = protocol.total_evaluations
    calculator = CountedCalculator(backend, total)
    numbers = start['atomic_numbers']
    nodes, mols, edges, attempts = [], [], [], []
    gsm_spent = 0
    report = dict(start=start, strategy=strategy, protocol=dataclasses.asdict(protocol), nodes=nodes,
                  edges=edges, attempts=attempts, status='running', evidence='MLIP_descents_not_DFT_IRC',
                  baseline=dict(engine='pyGSM SE_GSM', pygsm=str(args.pygsm), arm=args.arm,
                                num_nodes=args.num_nodes, gsm_attempt_cap=args.gsm_attempt_cap,
                                pygsm_commit=subprocess.run(['git', '-C', str(args.pygsm), 'rev-parse', 'HEAD'],
                                                            capture_output=True, text=True).stdout.strip()),
                  provenance=runtime_provenance(ROOT), potential=potential)
    started = time.perf_counter()

    def spent():
        return calculator.calls + gsm_spent

    def save():
        report.update(evaluations=spent(), gsm_evaluations=gsm_spent, in_process_evaluations=calculator.calls,
                      elapsed_seconds=time.perf_counter()-started)
        atomic_json(outdir/'network.json', report)

    def register(end, depth):
        x = np.asarray(end['positions_A'])
        mol = geometry_mol(numbers, x, start['charge'])
        for node, old in zip(nodes, mols):
            if (abs(node['energy_eV']-end['energy_eV']) <= protocol.energy_tolerance_eV and
                    molecular_rmsd(old, np.asarray(node['positions_A']), mol, x) < protocol.geometry_tolerance_A):
                return node['id'], mol
        nodes.append(dict(id=len(nodes), depth_discovered=depth, **end))
        mols.append(mol)
        return len(nodes)-1, mol

    try:
        root, checks, init = initialize_root(start, calculator, outdir, protocol)
        report.update(initial_relaxation_checks=checks, initial_point=dict(root), initialization_evaluations=init)
        if not root['full_system_minimum']:
            report['status'] = 'initial_minimum_unresolved'
            return
        register(root, 0)
        root_atoms = Atoms(numbers=numbers, positions=root['positions_A'])
        root_mol = mols[0]
        if args.arm == 'oracle':
            references = {r['id']: r for r in map(json.loads, args.oracle_references.read_text(encoding='utf-8').splitlines())}
            oracle = OracleReferenceLibrary(references[start['id']], start)
            sets = proposal_driving_sets(oracle.propose(root_mol))
            report['baseline']['oracle_policy'] = oracle.policy
        elif args.arm == 'edits':
            library = ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
            if protocol.proposal_resonance_forms:      # same proposals as the symbolic arms
                library = ResonanceAwareLibrary(library, protocol.proposal_resonance_forms)
            sets = proposal_driving_sets(library.propose(root_mol))
        else:
            rng = np.random.default_rng([protocol.random_seed, 11])
            sets = b2f2_driving_sets(np.asarray(numbers), np.asarray(root['positions_A']),
                                     list(bond_orders(root_mol)), rng, protocol.max_attempts)
        report['driving_sets_available'] = len(sets)
        for task, coords in enumerate(sets):
            if len(attempts) >= protocol.max_attempts or total - spent() < MIN_REMAINING:
                break
            attempt_dir = outdir/f'attempt_{task:03d}'
            limit = min(args.gsm_attempt_cap, total - spent())
            gsm_evals, code, ts_path, seconds = run_gsm(attempt_dir/'gsm', root_atoms, coords, limit,
                                                        task, args, start)
            gsm_spent += gsm_evals
            attempt = dict(id=task, source_node=0, driving_coordinates=coords, gsm_exit=code,
                           gsm_evaluations=gsm_evals, gsm_seconds=seconds, ts_node_found=ts_path is not None,
                           artifact=f'attempt_{task:03d}/post/post_result.json')
            post_evals = 0
            if ts_path is not None and total - spent() > 0:
                calculator.attempt_limit = calculator.calls + (total - spent())
                post_dir = attempt_dir/'post'
                post_dir.mkdir()
                result = post_process(ts_path, numbers, calculator, protocol, post_dir, start['charge'])
                post_evals = result['evaluations']
                attempt['status'] = result['status']
                if result['status'] in ('validated_descents', 'validated_core_descents'):
                    registered = [register(e, 1) for e in result['endpoints']]
                    ends = [r[0] for r in registered]
                    attempt['observed_nodes'] = ends
                    if ends[0] == ends[1]:
                        attempt['status'] = 'same_basin_return'
                    elif any(sorted(e['nodes']) == sorted(ends) and
                             abs(e['ts_energy_eV']-result['ts']['energy_eV']) < .03 and
                             aligned_rmsd(np.asarray(e['ts_positions_A']), np.asarray(result['ts_positions_A'])) < .15
                             for e in edges):
                        attempt['status'] = 'duplicate_connection'
                    else:
                        attempt['status'] = 'new_connection'
                        chemistry = edge_chemistry([r[1] for r in registered], [mols[i] for i in ends])
                        edges.append(dict(id=len(edges), nodes=ends, attempt=len(attempts), proposed_from=0,
                                          ts_energy_eV=result['ts']['energy_eV'],
                                          ts_positions_A=result['ts_positions_A'],
                                          barriers_eV=[e['barrier_eV'] for e in result['endpoints']],
                                          endpoint_acceptance=result['status'], **chemistry,
                                          kind=('conformational' if chemistry['endpoint_chemistry']['resonance_equivalent']
                                                else 'chemical')))
            else:
                attempt['status'] = 'gsm_no_ts_node' if ts_path is None else 'budget_exhausted'
            attempt['evaluations'] = gsm_evals + post_evals
            attempts.append(attempt)
            save()
            print(json.dumps(dict(start=start['id'], arm=args.arm, attempt=task, status=attempt['status'],
                                  evaluations=spent())), flush=True)
        report['status'] = 'completed'
        report['stop_reason'] = ('evaluation_budget' if total - spent() < MIN_REMAINING else
                                 'attempt_budget' if len(attempts) >= protocol.max_attempts else 'driving_sets_exhausted')
    except BudgetExceeded:
        report['status'] = 'initialization_budget_exhausted'
    except Exception as exc:
        report.update(status='aborted_error', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
