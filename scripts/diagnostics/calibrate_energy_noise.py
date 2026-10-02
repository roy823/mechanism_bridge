"""Calibrate MLIP energy noise under rigid motions to set min_barrier_eV.

Energies must be invariant to rotation and translation, so the spread of
E(R x + t) - E(x) over random rigid motions measures the numerical resolution
of the model (float precision, neighbour-list and summation order). TS and
endpoint geometries are sampled from saved network.json files; no new search.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.potentials import load_potential  # noqa: E402

SUPPORTED = {'aimnet2-rxn': {1, 6, 7, 8}}


def structures(paths, elements, rng, limit):
    """(label, numbers, positions, charge) for TSs and endpoints of saved edges."""
    pool = []
    for path in paths:
        network = json.loads(path.read_text(encoding='utf-8'))
        start = network['start']
        numbers = np.asarray(start['atomic_numbers'])
        if start['charge'] != 0 or start['multiplicity'] != 1 or not set(numbers) <= elements:
            continue
        for edge in network['edges']:
            rel = path.relative_to(ROOT).as_posix()
            pool.append((f"{rel}#edge{edge['id']}:ts", numbers, edge['ts_positions_A']))
            for node in edge['nodes']:
                pool.append((f"{rel}#node{node}", numbers, network['nodes'][node]['positions_A']))
    unique = {label: (label, numbers, positions) for label, numbers, positions in pool}
    labels = sorted(unique)
    chosen = rng.choice(len(labels), size=min(limit, len(labels)), replace=False)
    return [unique[labels[k]] for k in sorted(chosen)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--potential', default='aimnet2-rxn',
                        choices=['aimnet2-rxn', 'aimnet2', 'aimnet2-2025'])
    parser.add_argument('--networks', default='reports/**/network.json',
                        help='Glob (relative to the repository) of saved networks to sample')
    parser.add_argument('--structures', type=int, default=60)
    parser.add_argument('--motions', type=int, default=20)
    parser.add_argument('--seed', type=int, default=20261002)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    backend, provenance = load_potential(args.potential, ROOT, threads=args.threads)
    elements = SUPPORTED.get(args.potential, set(backend.supported_species)
                             if hasattr(backend, 'supported_species') else {1, 6, 7, 8})
    rng = np.random.default_rng(args.seed)
    paths = sorted(ROOT.glob(args.networks))
    rows = []
    for label, numbers, positions in structures(paths, elements, rng, args.structures):
        x = np.asarray(positions, dtype=float)
        motions = [x]
        for _ in range(args.motions):
            rotation = Rotation.random(random_state=rng).as_matrix()
            shift = rng.uniform(-5., 5., size=3)
            motions.append((x - x.mean(0)) @ rotation.T + x.mean(0) + shift)
        # The barrier gate compares single-point energies, so measure those; the
        # batched path (finite-difference Hessians) is reported separately.
        energies = []
        for y in motions:
            atoms = Atoms(numbers=numbers, positions=y)
            atoms.calc = backend
            energies.append(float(atoms.get_potential_energy()))
        energies = np.asarray(energies)
        batch = np.asarray(backend.evaluate_many(numbers, np.asarray(motions))['energy'], dtype=float)
        deltas = energies[1:] - energies[0]
        rows.append(dict(structure=label, atoms=int(len(numbers)), energy_eV=float(energies[0]),
                         std_eV=float(deltas.std()), max_abs_eV=float(np.abs(deltas).max()),
                         batch_minus_single_max_abs_eV=float(np.abs(batch - energies).max())))
    if not rows:
        raise SystemExit('No supported structures were found')
    pooled = float(np.sqrt(np.mean([r['std_eV']**2 for r in rows])))
    worst = float(max(r['max_abs_eV'] for r in rows))
    summary = dict(potential=args.potential, model=provenance, networks_glob=args.networks,
                   networks=len(paths), structures=len(rows), motions=args.motions, seed=args.seed,
                   pooled_std_eV=pooled, max_abs_eV=worst,
                   p99_max_abs_eV=float(np.percentile([r['max_abs_eV'] for r in rows], 99)),
                   batch_minus_single_max_abs_eV=float(max(r['batch_minus_single_max_abs_eV'] for r in rows)),
                   suggested_min_barrier_eV=max(1e-3, 5*pooled),
                   rule='min_barrier_eV = max(1e-3 eV, 5 x pooled rigid-motion std)', rows=rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('rows', 'model')}, indent=2))


if __name__ == '__main__':
    main()
