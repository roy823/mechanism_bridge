"""Pdft-0: CPU PySCF versus GPU4PySCF for the DFT verification level.

Same functional, basis, grid level, SCF convergence and (no) density fitting as
backends.PySCFCalculator. For each geometry: wall time and results of the SCF
energy, analytic gradient and (optionally) analytic Hessian on each device.
The CPU reference runs on the CPU partition that the CPU route of Fig. 6 would
use, the GPU side on one GPU; --merge compares the two result files. The
decision rule (E doc, Pdft-0): use the GPU for Fig. 6 only if |dE| and
gradients agree (barrier error <= 0.1 kcal/mol) and the GPU is >= 5x faster
than the multithreaded CPU.
Usage: benchmark_gpu_dft.py --devices cpu|gpu [...] --out F.json
       benchmark_gpu_dft.py --merge CPU.json GPU.json --out F.json
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


def run(device, numbers, positions, args, hessian):
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
    if hessian:
        started = time.perf_counter()
        matrix = mf.Hessian().kernel()
        matrix = np.asarray(matrix.get() if hasattr(matrix, 'get') else matrix)
        timings['hessian_s'] = time.perf_counter() - started
        out['hessian'] = matrix.tolist()
    out['timings'] = timings
    out['nao'] = int(mol.nao_nr())
    return out


def device_info(device):
    if device == 'cpu':
        cpuinfo = Path('/proc/cpuinfo')
        names = [l.split(':', 1)[1].strip() for l in cpuinfo.read_text().splitlines()
                 if l.startswith('model name')] if cpuinfo.exists() else []
        return dict(model=names[0] if names else None)
    import cupy
    props = cupy.cuda.runtime.getDeviceProperties(0)
    return dict(name=props['name'].decode() if isinstance(props['name'], bytes) else props['name'],
                compute_capability=cupy.cuda.Device(0).compute_capability, cupy=cupy.__version__)


def merge(cpu_file, gpu_file):
    """Differences and speedups for systems that both devices finished."""
    cpu = {r['system']: r for r in json.loads(cpu_file.read_text())['rows']}
    gpu = {r['system']: r for r in json.loads(gpu_file.read_text())['rows']}
    rows = []
    for system in [s for s in cpu if s in gpu]:
        c, g = cpu[system]['cpu'], gpu[system]['gpu']
        row = dict(system=system, atoms=cpu[system]['atoms'], cpu_timings=c.get('timings'),
                   gpu_timings=g.get('timings'), cpu_error=c.get('error'), gpu_error=g.get('error'))
        if 'energy_Ha' in c and 'energy_Ha' in g:
            row['abs_energy_difference_kcal_mol'] = abs(c['energy_Ha'] - g['energy_Ha'])*HARTREE_KCAL
            row['max_gradient_difference_Ha_bohr'] = float(np.abs(np.subtract(c['gradient'], g['gradient'])).max())
            if 'hessian' in c and 'hessian' in g:
                row['max_hessian_difference'] = float(np.abs(np.subtract(c['hessian'], g['hessian'])).max())
            shared = [k for k in c['timings'] if k in g['timings']]
            row['speedup'] = {k: c['timings'][k]/g['timings'][k] for k in shared}
            row['speedup_total_shared_steps'] = sum(c['timings'][k] for k in shared)/sum(g['timings'][k] for k in shared)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--devices', nargs='+', choices=('cpu', 'gpu'), default=['cpu', 'gpu'])
    parser.add_argument('--xc', default='wb97x')
    parser.add_argument('--basis', default='6-31g(d)')
    parser.add_argument('--grid-level', type=int, default=3)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--limit-atoms', type=int, default=40)
    parser.add_argument('--hessian', action='store_true')
    parser.add_argument('--hessian-max-atoms', type=int, help='Skip the Hessian above this size')
    parser.add_argument('--fail-fast', action='store_true', help='Stop at the first error (frees a GPU early)')
    parser.add_argument('--merge', nargs=2, type=Path, metavar=('CPU_JSON', 'GPU_JSON'))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.merge:
        rows = merge(*args.merge)
        args.out.write_text(json.dumps(dict(cpu=str(args.merge[0]), gpu=str(args.merge[1]), rows=rows), indent=2))
        print(json.dumps(rows, indent=2))
        return
    info = {d: device_info(d) for d in args.devices}
    print(json.dumps(info), flush=True)
    rows = []
    def save():
        args.out.write_text(json.dumps(dict(xc=args.xc, basis=args.basis, grid_level=args.grid_level,
                                            threads=args.threads, devices=info, rows=rows), indent=2))
    for name, numbers, positions in geometries(args.limit_atoms):
        row = dict(system=name, atoms=len(numbers))
        hessian = args.hessian and (args.hessian_max_atoms is None or len(numbers) <= args.hessian_max_atoms)
        for device in args.devices:
            try:
                row[device] = run(device, numbers, positions, args, hessian)
            except Exception as exc:          # record and continue: this is a benchmark
                row[device] = dict(error=f'{type(exc).__name__}: {exc}')
                if args.fail_fast:
                    rows.append(row)
                    save()
                    raise
        rows.append(row)
        print(json.dumps({k: (v if not isinstance(v, dict) else {x: y for x, y in v.items()
                                                                  if x not in ('gradient', 'hessian')})
                          for k, v in row.items()}), flush=True)
        save()


if __name__ == '__main__':
    main()
