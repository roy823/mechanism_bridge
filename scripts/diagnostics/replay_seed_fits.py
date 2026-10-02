"""Replay symbolic seed-geometry fits under larger iteration caps (seed_fit_max_nfev).

Failed fits: every seed_generation_failed attempt of a campaign is rebuilt from
its source node (geometry, graph), proposal, strategy and protocol. The record
does not store the geometric variant, so variants 0-2 are tried and only those
that fail under the run's own cap are kept; each is refit with every --caps
value. Converged fits (controls): up to --controls successful symbolic attempts
per run are rebuilt from their recorded seed metadata. The replay must give the
recorded RNG seed, nfev and cost under the run's cap, and the seed under each
larger cap must be bit-identical, since a converged fit stops before the cap.
Array use: --task K --tasks N takes every N-th run of the manifest.
Usage: replay_seed_fits.py CAMPAIGN_DIR --out FILE [--caps 1000 3000] [--controls 3]
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from ase import Atoms
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol  # noqa: E402
from mechbridge.search_seeds import make_seed  # noqa: E402

FRACTIONS = (0.35, 0.55, 0.75)


def seed_rng(node, positions, proposal, random_seed, variant):
    """The attempt RNG seed exactly as reaction_network.explore derives it."""
    identity = json.dumps(dict(graph=node['graph_smiles'], positions=np.round(positions, 6).tolist(),
                               edits=proposal['edits'] if proposal else None,
                               seed=random_seed, variant=variant), sort_keys=True)
    return int.from_bytes(hashlib.sha256(identity.encode()).digest()[:4], 'little')


def refit(network, node_id, strategy, proposal, variant, cap):
    protocol = network['protocol']
    numbers = network['start']['atomic_numbers']
    node = network['nodes'][node_id]
    positions = np.asarray(node['positions_A'], dtype=float)
    state = Atoms(numbers=numbers, positions=positions)
    mol = geometry_mol(numbers, positions, network['start']['charge'])
    rng = seed_rng(node, state.positions, proposal, protocol['random_seed'], variant)
    started = time.perf_counter()
    try:
        x, _, meta = make_seed(state, mol, strategy, proposal, variant, rng, protocol['symbolic_seed_scale'],
                               protocol['encounter_policy'], protocol['seed_features'], cap,
                               protocol.get('seed_feature_ablation', 'none'))
    except (ValueError, RuntimeError, ArithmeticError) as exc:     # recoverable in explore()
        return dict(ok=False, error=str(exc), seconds=time.perf_counter() - started, rng=rng)
    return dict(ok=True, nfev=int(meta['geometric_fit_nfev']), cost=float(meta['geometric_fit_cost']),
                seconds=time.perf_counter() - started, rng=rng, x=x)


def original_proposal(meta):
    """Proposal fields that make_seed reads, from the recorded seed metadata."""
    return dict(edits=meta['net_edits'], arrows=meta.get('arrows'), predicted_graph=meta['predicted_graph'],
                template_id=meta['template_id'], template_source_graph=meta.get('template_source_graph'),
                transferred=meta.get('transferred'), origin=meta.get('origin'),
                pattern_smarts=meta.get('pattern_smarts'), intermolecular=meta.get('intermolecular', False))


def strip(result):
    return {k: v for k, v in result.items() if k != 'x'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--caps', type=int, nargs='+', default=[1000, 3000])
    parser.add_argument('--controls', type=int, default=3, help='Successful symbolic attempts per run')
    parser.add_argument('--task', type=int, default=0)
    parser.add_argument('--tasks', type=int, default=1)
    args = parser.parse_args()
    if not 0 <= args.task < args.tasks:
        parser.error('--task must be in [0, --tasks)')
    RDLogger.DisableLog('rdApp.*')
    manifest = list(csv.DictReader((args.campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'))
    failures, controls = [], []
    for task in manifest[args.task::args.tasks]:
        if task['strategy'] not in ('bond_edits', 'arrows'):
            continue
        paths = list((args.campaign/'runs'/task['outdir']).glob('*/*/network.json'))
        if not paths:
            continue
        network = json.loads(paths[0].read_text(encoding='utf-8'))
        run_cap = network['protocol'].get('seed_fit_max_nfev', 200)
        strategy, done = task['strategy'], 0
        for attempt in network['attempts']:
            record = dict(run=task['outdir'], attempt=attempt['id'], node=attempt['source_node'], strategy=strategy)
            if attempt['status'] == 'seed_generation_failed':
                proposal = attempt['proposal']
                record['template_id'] = proposal.get('template_id')
                for variant in range(3):
                    base = refit(network, attempt['source_node'], strategy, proposal, variant, run_cap)
                    if base['ok']:
                        continue                      # not the variant that failed in the run
                    row = dict(record, variant=variant, run_cap=run_cap, at_run_cap=strip(base))
                    row['caps'] = {str(c): strip(refit(network, attempt['source_node'], strategy, proposal,
                                                       variant, c)) for c in args.caps}
                    failures.append(row)
                    print(json.dumps(row), flush=True)
            elif done < args.controls and 'net_edits' in (attempt.get('proposal') or {}):
                meta = attempt['proposal']
                variant = FRACTIONS.index(meta['fraction'])
                proposal = original_proposal(meta)
                base = refit(network, attempt['source_node'], strategy, proposal, variant, run_cap)
                row = dict(record, template_id=meta['template_id'], variant=variant,
                           reproduces=bool(base['ok'] and base['rng'] == meta['random_seed'] and
                                           base['nfev'] == meta['geometric_fit_nfev'] and
                                           base['cost'] == meta['geometric_fit_cost']))
                row['identical_under_caps'] = {}
                for cap in args.caps:
                    other = refit(network, attempt['source_node'], strategy, proposal, variant, cap)
                    row['identical_under_caps'][str(cap)] = bool(base['ok'] and other['ok'] and
                                                                 np.array_equal(base['x'], other['x']))
                controls.append(row)
                done += 1
    summary = dict(failed_fits=len(failures), controls=len(controls),
                   controls_reproduced=sum(c['reproduces'] for c in controls),
                   controls_identical={str(c): sum(r['identical_under_caps'][str(c)] for r in controls)
                                       for c in args.caps},
                   rescued={str(c): sum(f['caps'][str(c)]['ok'] for f in failures) for c in args.caps})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(dict(campaign=str(args.campaign), caps=args.caps, task=args.task,
                                        tasks=args.tasks, summary=summary, failures=failures,
                                        controls=controls), indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
