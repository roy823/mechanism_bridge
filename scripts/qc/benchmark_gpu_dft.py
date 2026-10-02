"""Pdft-0: CPU PySCF versus GPU4PySCF for the DFT verification level.

Same functional, basis, grid level, SCF convergence and (no) density fitting as
backends.PySCFCalculator. For each geometry: wall time and results of the SCF
energy, analytic gradient and (optionally) analytic Hessian on both devices,
and the largest differences. The decision rule (E doc, Pdft-0): use the GPU for
Fig. 6 only if |dE| and gradients agree (barrier error <= 0.1 kcal/mol) and the
GPU is >= 5x faster than the multithreaded CPU on the same node.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HARTREE_KCAL = 627.509474


def geometries(limit_atoms):
    """Two DFT-verified TS, one T1x TS and the largest Coley start (<= limit_atoms)."""
    out = []
    for name in ('reports/bimolecular_v5/qc/formaldehyde_dimer_refined', 'reports/local_transfer_v3/qc/acetone'):
        source = json.loads((ROOT/name/'verification.json').read_text())['source']
        out.append((Path(name).name + '_ts', source['atomic_numbers'], source['positions_A']['ts']))
    t1x = json.loads((ROOT/'data/processed/t1x_test_references.jsonl').read_text().splitlines()[0])
    out.append((t1x['id'] + '_ts', t1x['atomic_numbers'], t1x['ts_positions_A']))
    coley = [json.loads(l) for l in (ROOT/'data/processed/coley_benchmark_starts.jsonl').read_text().splitlines()]
    coley = max((s for s in coley if len(s['atomic_numbers']) <= limit_atoms), key=lambda s: len(s['atomic_numbers']))
    out.append((coley['id'] + '_start', coley['atomic_numbers'], coley['positions_A']))
    return out


def run(device, numbers, positions, args):
    from pyscf import gto, lib
    lib.num_threads(args.threads)
    symbols = [gto.mole._symbol(int(z)) for z in numbers]
    mol = gto.M(atom=list(zip(symbols, np.asarray(positions).tolist())), unit='Angstrom', charge=0, spin=0,
                basis=args.basis, verbose=0)
    if device == 'gpu':
        from gpu4pyscf import dft
    else:
        from pyscf import dft
    mf = dft.RKS(mol, xc=args.xc)
    mf.grids.level = args.grid_level
    mf.conv_tol, mf.max_cycle = 1e-10, 150
    timings, out = {}, {}
    started = time.perf_counter()
    energy = float(mf.kernel())
    timings['scf_s'] = time.perf_counter() - started
    if not mf.converged:
        raise RuntimeError(f'{device} SCF did not converge')
    started = time.perf_counter()
    gradient = mf.nuc_grad_method().kernel()
    gradient = np.asarray(gradient.get() if hasattr(gradient, 'get') else gradient)
    timings['gradient_s'] = time.perf_counter() - started
    out.update(energy_Ha=energy, gradient=gradient.tolist())
    if args.hessian:
        started = time.perf_counter()
        hessian = mf.Hessian().kernel()
        hessian = np.asarray(hessian.get() if hasattr(hessian, 'get') else hessian)
        timings['hessian_s'] = time.perf_counter() - started
        out['hessian_norm'] = float(np.linalg.norm(hessian))
        out['hessian'] = hessian
    out['timings'] = timings
    out['nao'] = int(mol.nao_nr())
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xc', default='wb97x')
    parser.add_argument('--basis', default='6-31g(d)')
    parser.add_argument('--grid-level', type=int, default=3)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--limit-atoms', type=int, default=40)
    parser.add_argument('--hessian', action='store_true')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for name, numbers, positions in geometries(args.limit_atoms):
        row = dict(system=name, atoms=len(numbers))
        results = {}
        for device in ('cpu', 'gpu'):
            try:
                results[device] = run(device, numbers, positions, args)
                row[device] = dict(timings=results[device]['timings'], nao=results[device]['nao'],
                                   energy_Ha=results[device]['energy_Ha'])
            except Exception as exc:          # record and continue: this is a benchmark
                row[device] = dict(error=f'{type(exc).__name__}: {exc}')
        if all('energy_Ha' in row.get(d, {}) for d in ('cpu', 'gpu')):
            c, g = results['cpu'], results['gpu']
            row['abs_energy_difference_kcal_mol'] = abs(c['energy_Ha'] - g['energy_Ha']) * HARTREE_KCAL
            row['max_gradient_difference_Ha_bohr'] = float(np.abs(np.asarray(c['gradient']) -
                                                              np.asarray(g['gradient'])).max())
            if args.hessian:
                row['max_hessian_difference'] = float(np.abs(c['hessian'] - g['hessian']).max())
            total = lambda r: sum(r['timings'].values())
            row['speedup_total'] = total(c) / total(g)
        rows.append(row)
        print(json.dumps(row), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(dict(xc=args.xc, basis=args.basis, grid_level=args.grid_level,
                                        cpu_threads=args.threads, hessian=args.hessian, rows=rows), indent=2))


if __name__ == '__main__':
    main()
