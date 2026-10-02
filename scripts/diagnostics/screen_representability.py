"""Representability screen: does the MLIP reproduce each reference TS and reaction?

For every admitted reference reaction: Sella P-RFO (order=1) from the reference
TS with a capped budget that is billed to the screen, never to a search method;
index-one check with the finite-difference Hessian; bidirectional MLIP IRC plus
BFGS polish; endpoint connectivities compared with the reference reactant and
product under one reactant automorphism. Only representable reactions enter
the main Fig. 4 denominator; the others are reported separately.

Array use: --task K --tasks N processes every N-th admitted reaction.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.potentials import load_potential  # noqa: E402
from mechbridge.reaction_network import (BudgetExceeded, CountedCalculator, SearchProtocol,  # noqa: E402
                                         aligned_rmsd, inspect_point, integrate_irc,
                                         is_recoverable_failure)
from mechbridge.reference_matching import (automorphisms, bond_set, edge_matches, mol_bonds,  # noqa: E402
                                           permute, stereo_free)


def best_rmsd(x, y, maps):
    """Minimum Kabsch RMSD over reactant automorphisms (atom relabelling)."""
    x, y = np.asarray(x), np.asarray(y)
    values = []
    for mapping in maps:
        order = [mapping[i] for i in range(len(x))]
        values.append(aligned_rmsd(x[order], y))
    return min(values) if values else aligned_rmsd(x, y)


def screen(reference, backend, protocol, outdir, prfo_steps, budget):
    numbers = reference['atomic_numbers']
    maps = automorphisms(numbers, reference['reactant_bonds'])
    react, prod = bond_set(reference['reactant_bonds']), bond_set(reference['product_bonds'])
    calculator = CountedCalculator(backend, budget)
    row = dict(id=reference['id'], atoms=len(numbers), representable=False, stage='started')
    started = time.time()
    outdir.mkdir(parents=True, exist_ok=False)
    try:
        from sella import Sella
        ts = Atoms(numbers=numbers, positions=reference['ts_positions_A'])
        ts.calc = calculator
        with Sella(ts, order=1, internal=False, logfile=str(outdir/'prfo.log'),
                   trajectory=str(outdir/'prfo.traj')) as opt:
            row['prfo_converged'] = bool(opt.run(fmax=protocol.fmax, steps=prfo_steps))
        row['prfo_evaluations'] = calculator.calls
        analysis, modes = inspect_point(ts, protocol)
        frequencies = analysis['frequencies_cm-1']
        row.update(ts_force_max_eV_A=analysis['force_max_eV_A'],
                   imaginary_count=analysis['imaginary_count'],
                   lowest_frequencies_cm=sorted(frequencies)[:2],
                   ts_energy_eV=analysis['energy_eV'],
                   ts_rmsd_to_reference_A=best_rmsd(ts.positions, reference['ts_positions_A'], maps))
        if not analysis['force_converged'] or analysis['imaginary_count'] != 1:
            row['stage'] = 'not_index_one'
            return row
        ends = []
        for direction in ('reverse', 'forward'):
            end = ts.copy()
            end.calc = calculator
            info = integrate_irc(end, outdir, protocol, direction)
            with BFGS(end, maxstep=.1, logfile=str(outdir/f'polish_{direction}.log')) as opt:
                opt.run(fmax=protocol.fmax, steps=protocol.descent_steps)
            mol = geometry_mol(numbers, end.positions, 0)
            ends.append(dict(direction=direction, irc_converged=info['irc_converged'],
                             graph=graph_smiles(mol), key=stereo_free(graph_smiles(mol)),
                             bonds=mol_bonds(mol), energy_eV=float(end.get_potential_energy())))
        mapped = edge_matches(ends[0]['bonds'], ends[1]['bonds'], react, prod, maps)
        unmapped = sorted(e['key'] for e in ends) == sorted([reference['reactant_key'],
                                                             reference['product_key']])
        row.update(stage='screened', irc_converged=[e['irc_converged'] for e in ends],
                   endpoint_graphs=[e['graph'] for e in ends], mapped_match=mapped,
                   unmapped_match=unmapped, representable=bool(mapped),
                   barrier_from_reactant_side_eV=None)
        if mapped:
            # Orient: which IRC end is the reactant? (same automorphism test as above)
            reactant_end = next(e for e in ends if any(permute(e['bonds'], m) == react for m in maps))
            row['barrier_from_reactant_side_eV'] = row['ts_energy_eV'] - reactant_end['energy_eV']
    except BudgetExceeded:
        row['stage'] = 'budget_exhausted'
    except Exception as exc:   # perception, nonfinite output, optimizer failure: not representable
        if not is_recoverable_failure(exc):
            raise
        row.update(stage='calculation_failed', error=f'{type(exc).__name__}: {exc}')
    finally:
        row.update(evaluations=calculator.calls, seconds=time.time()-started)
        (outdir/'screen.json').write_text(json.dumps(row, indent=2, default=list), encoding='utf-8')
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--potential', default='aimnet2-rxn')
    parser.add_argument('--prfo-steps', type=int, default=200)
    parser.add_argument('--budget', type=int, default=4000,
                        help='Screen-only evaluation cap per reaction (P-RFO + Hessian + IRC)')
    parser.add_argument('--task', type=int, default=0)
    parser.add_argument('--tasks', type=int, default=1)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    references = [r for r in map(json.loads, args.references.read_text(encoding='utf-8').splitlines())
                  if r.get('admission') == 'accepted']
    backend, provenance = load_potential(args.potential, ROOT, threads=args.threads)
    validator = getattr(backend, 'validate_system', None)
    protocol = SearchProtocol(hessian_batch_size=32)
    rows = []
    for reference in references[args.task::args.tasks]:
        if validator is not None:
            validator(reference['atomic_numbers'], 0, 1)
        row = screen(reference, backend, protocol, args.out/'reactions'/reference['id'],
                     args.prfo_steps, args.budget)
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ('endpoint_graphs',)}, default=list),
              flush=True)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out/f'task_{args.task:03d}.json').write_text(json.dumps(dict(
        potential=provenance, protocol=dict(fmax=protocol.fmax, hessian_step=protocol.hessian_step,
        irc_dx=protocol.irc_dx, irc_steps=protocol.irc_steps, prfo_steps=args.prfo_steps,
        budget=args.budget), rows=rows), indent=2, default=list), encoding='utf-8')


if __name__ == '__main__':
    main()
