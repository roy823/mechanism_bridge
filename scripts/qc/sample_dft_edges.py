"""Pre-registered Fig. 6 sample: stratified random MLIP edges for DFT verification.

Frame: every new edge whose attempt result is 'validated_descents' in the given
campaign directories. Duplicates across runs (same unordered graph pair, TS
aligned RMSD <= 0.15 A and |dE| <= 0.03 eV) keep the first occurrence in
manifest order. Strata are strategies; each gets --per-stratum edges; a short
stratum is taken whole and the remainder is spread over the others in
proportion to their size. The seed is fixed before any DFT calculation.
Usage: sample_dft_edges.py CAMPAIGN_DIR [...] --out FILE [--per-stratum 16 --seed 20261002]
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reaction_network import aligned_rmsd  # noqa: E402


def frame(campaigns):
    edges, kept = [], []
    for campaign in campaigns:
        manifest = list(csv.DictReader((campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'))
        for task in manifest:
            run = campaign/'runs'/task['outdir']
            for path in sorted(run.glob('*/*/network.json')):
                network = json.loads(path.read_text(encoding='utf-8'))
                for edge in network['edges']:
                    attempt = network['attempts'][edge['attempt']]
                    result = json.loads((path.parent/attempt['artifact']).read_text(encoding='utf-8'))
                    if result.get('status') != 'validated_descents':
                        continue
                    pair = tuple(sorted(network['nodes'][i]['graph_smiles'] for i in edge['nodes']))
                    edges.append(dict(campaign=str(campaign), run=str(run), start=task['start_id'],
                                      strategy=task['strategy'], seed=int(task['seed']), edge=edge['id'],
                                      kind=edge['kind'], pair=pair, ts_energy_eV=edge['ts_energy_eV'],
                                      ts=np.asarray(edge['ts_positions_A'])))
    for edge in edges:
        if any(k['start'] == edge['start'] and k['pair'] == edge['pair'] and
               abs(k['ts_energy_eV'] - edge['ts_energy_eV']) <= .03 and aligned_rmsd(k['ts'], edge['ts']) <= .15
               for k in kept):
            continue
        kept.append(edge)
    return edges, kept


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaigns', type=Path, nargs='+')
    parser.add_argument('--per-stratum', type=int, default=16)
    parser.add_argument('--seed', type=int, default=20261002)
    parser.add_argument('--chemical-only', action='store_true', help='Exclude conformational edges')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    all_edges, unique = frame(args.campaigns)
    if args.chemical_only:
        unique = [e for e in unique if e['kind'] == 'chemical']
    strata = {}
    for edge in unique:
        strata.setdefault(edge['strategy'], []).append(edge)
    target = args.per_stratum*len(strata)
    take = {s: min(args.per_stratum, len(v)) for s, v in strata.items()}
    spare = target - sum(take.values())
    while spare > 0:
        room = {s: len(strata[s]) - take[s] for s in strata if len(strata[s]) > take[s]}
        if not room:
            break
        for s in sorted(room, key=lambda s: -room[s]):
            if spare and take[s] < len(strata[s]):
                take[s] += 1
                spare -= 1
    rng = np.random.default_rng(args.seed)
    sample = []
    for s in sorted(strata):
        picks = rng.choice(len(strata[s]), size=take[s], replace=False)
        sample += [strata[s][k] for k in sorted(picks)]
    rows = [{k: v for k, v in e.items() if k != 'ts'} for e in sample]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(dict(seed=args.seed, per_stratum=args.per_stratum,
                                        chemical_only=args.chemical_only, frame_edges=len(all_edges),
                                        unique_edges=len(unique), strata={s: len(v) for s, v in strata.items()},
                                        taken=take, sample=rows), indent=2), encoding='utf-8')
    print(json.dumps(dict(frame=len(all_edges), unique=len(unique), taken=take), indent=2))


if __name__ == '__main__':
    main()
