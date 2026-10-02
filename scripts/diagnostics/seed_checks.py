"""P0(b) seed effectiveness and P4 seed-count check on a pilot campaign.

P0(b) (E doc 1.3): for every strategy, two seeds of the same start must give
different attempt sequences on at least one start that has attempts. An
attempt is identified by its proposal (template, active atoms, fraction,
displacement), status and evaluation count; a strategy whose seeds never
differ produces pseudo-replicates.

P4: from the per-system seed spread of pairs_at_budget (the Fig. 3 metric),
the seeds S per (system, strategy) needed for a paired test across N systems
to detect a 20% relative difference at alpha=0.05 and power 0.8. Model of a
per-system difference of seed means: d_i = delta + u_i + e_i with
Var(u) = tau^2 (true between-system heterogeneity) and
Var(e) = (sA^2 + sB^2)/S (pooled within-system seed variances). Requirement:
delta^2 N/(z_{1-a/2} + z_{power})^2 >= tau^2 + (sA^2 + sB^2)/S. tau^2 is
estimated from the pilot as Var(d) - (sA^2 + sB^2)/S_pilot (floored at 0).
'seeds_noise_only' sets tau^2 = 0: the seeds needed if seed noise were the
only obstacle. 'heterogeneity_limited' marks pairs where no S suffices.
Usage: seed_checks.py CAMPAIGN_DIR --out FILE [--ids F.json] [--systems N] [--budget 16000]
"""
import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]


def load_summary():
    path = ROOT/'scripts/exploration/summarize_campaign.py'
    spec = importlib.util.spec_from_file_location('summarize_campaign', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def signature(network):
    """Seed-sensitive identity of the attempt sequence (timings and RNG seeds excluded)."""
    rows = []
    for attempt in network['attempts']:
        proposal = attempt.get('proposal') or {}
        rows.append((proposal.get('template_id'), tuple(proposal.get('active_atoms') or ()),
                     proposal.get('fraction'), proposal.get('displacement_norm_A'),
                     attempt['status'], attempt['evaluations']))
    return rows


def required_seeds(per_seed, a, b, systems, relative=.2, alpha=.05, power=.8):
    common = sorted(set(per_seed.get(a, {})) & set(per_seed.get(b, {})))
    common = [s for s in common if len(per_seed[a][s]) > 1 and len(per_seed[b][s]) > 1]
    if len(common) < 2:
        return dict(systems_in_pilot=len(common), note='fewer than two systems with two or more seeds')
    pilot_seeds = min(min(len(per_seed[a][s]), len(per_seed[b][s])) for s in common)
    var_a = float(np.mean([np.var(per_seed[a][s], ddof=1) for s in common]))
    var_b = float(np.mean([np.var(per_seed[b][s], ddof=1) for s in common]))
    means_a = np.array([np.mean(per_seed[a][s]) for s in common])
    means_b = np.array([np.mean(per_seed[b][s]) for s in common])
    level = float(np.mean((means_a + means_b)/2))
    out = dict(systems_in_pilot=len(common), pilot_seeds=pilot_seeds, planned_systems=systems,
               seed_variance=dict(zip((a, b), (var_a, var_b))), mean_level=level,
               mean_difference=float(np.mean(means_a - means_b)))
    if level <= 0:
        return dict(out, note='both strategies find nothing on these systems')
    delta = relative*level
    tau2 = max(0., float(np.var(means_a - means_b, ddof=1)) - (var_a + var_b)/pilot_seeds)
    capacity = float(delta**2*systems/(norm.ppf(1 - alpha/2) + norm.ppf(power))**2)
    def seeds(room):
        if room <= 0:
            return None                       # no seed count reaches the power
        return max(1, math.ceil((var_a + var_b)/room))
    out.update(delta=delta, tau2=tau2, seeds_required=seeds(capacity - tau2),
               seeds_noise_only=seeds(capacity), heterogeneity_limited=bool(capacity <= tau2))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--ids', type=Path, help='JSON list or {"ids": [...]} restricting the P4 systems')
    parser.add_argument('--systems', type=int, help='Planned number of systems N (default: pilot systems)')
    parser.add_argument('--budget', type=int, default=16000)
    args = parser.parse_args()
    summary = load_summary()
    keep = None
    if args.ids:
        data = json.loads(args.ids.read_text(encoding='utf-8'))
        keep = set(data if isinstance(data, list) else data['ids'])
    manifest = list(csv.DictReader((args.campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'))
    runs, missing = {}, []
    for task in manifest:
        paths = list((args.campaign/'runs'/task['outdir']).glob('*/*/network.json'))
        if not paths:
            missing.append(task['task_id'])
            continue
        network = json.loads(paths[0].read_text(encoding='utf-8'))
        if network['status'] in ('running', 'aborted_error'):
            missing.append(task['task_id'])
            continue
        runs[task['start_id'], task['strategy'], int(task['seed'])] = network
    p0b = {}
    for strategy in sorted({k[1] for k in runs}):
        starts = sorted({k[0] for k in runs if k[1] == strategy})
        compared = differing = 0
        for start in starts:
            seeds = sorted(k[2] for k in runs if k[:2] == (start, strategy))
            if len(seeds) < 2:
                continue
            sequences = [signature(runs[start, strategy, s]) for s in seeds]
            if not any(sequences):
                continue                      # no attempts under any seed: nothing to compare
            compared += 1
            differing += any(seq != sequences[0] for seq in sequences[1:])
        p0b[strategy] = dict(starts_with_attempts=compared, starts_where_seeds_differ=differing,
                             passes=compared > 0 and differing > 0)
    per_seed = {}
    for (start, strategy, seed), network in sorted(runs.items()):
        if keep is not None and start not in keep:
            continue
        value = summary.value_at(summary.anytime(network), args.budget)
        per_seed.setdefault(strategy, {}).setdefault(start, []).append(value)
    systems = args.systems or len({s for d in per_seed.values() for s in d})
    pairs = [(s, 'geometry') for s in sorted(per_seed) if s != 'geometry']
    if 'arrows' in per_seed and 'bond_edits' in per_seed:
        pairs.append(('arrows', 'bond_edits'))
    p4 = {f'{a}_vs_{b}': required_seeds(per_seed, a, b, systems) for a, b in pairs}
    result = dict(campaign=str(args.campaign), budget=args.budget, finished_runs=len(runs),
                  missing_or_unfinished_tasks=missing, p0b=p0b, p4=p4, pairs_at_budget=per_seed)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'pairs_at_budget'}, indent=2))


if __name__ == '__main__':
    main()
