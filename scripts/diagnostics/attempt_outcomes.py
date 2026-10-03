"""Attempt outcomes per strategy for a campaign (explains Fig. 3 efficiency differences).

Per run (the network.json files of a campaign laid out as for summarize_campaign.py):
attempts, evaluations per attempt and attempt statuses. For attempts whose IRC joined
two minima (new_connection, duplicate_connection): whether the source minimum is one
end, and whether the two end graphs are the source graph and the proposal's predicted
graph (stereo-free; 'intended'). Redundancy: the share of attempts whose (source graph,
predicted graph) pair was already proposed earlier in the same run.
Aggregates per strategy (pooled over attempts, plus per-system means of the per-run
values) and per proposal origin within each strategy (e.g. transferred published
templates versus grammar arrows). Runs still running or aborted are skipped.
Usage: attempt_outcomes.py CAMPAIGN_DIR [...] --out FILE [--ids F.json]
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reference_matching import stereo_free  # noqa: E402

CONNECTED = ('new_connection', 'duplicate_connection')


def plain(smiles):
    try:
        return stereo_free(smiles)
    except (ValueError, TypeError):
        return smiles


def run_outcomes(network):
    graphs = [plain(n['graph_smiles']) for n in network['nodes']]
    seen, rows = set(), []
    for attempt in network['attempts']:
        proposal = attempt.get('proposal') or {}
        source = graphs[attempt['source_node']] if attempt.get('source_node') is not None else None
        predicted = plain(proposal['predicted_graph']) if proposal.get('predicted_graph') else None
        ends = [graphs[i] for i in attempt.get('observed_nodes') or []]
        connected = attempt['status'] in CONNECTED and len(ends) == 2
        rows.append(dict(status=attempt['status'], evaluations=attempt['evaluations'],
                         origin=proposal.get('origin') or 'none', connected=connected,
                         involves_source=connected and bool(attempt.get('source_is_endpoint')),
                         intended=(connected and None not in (source, predicted)
                                   and sorted(ends) == sorted([source, predicted])),
                         repeated=(source, predicted) in seen))
        seen.add((source, predicted))
    return rows


def block(attempts):
    n = len(attempts)
    connected = [a for a in attempts if a['connected']]
    statuses = Counter(a['status'] for a in attempts)
    return dict(attempts=n,
                evaluations_per_attempt=float(np.mean([a['evaluations'] for a in attempts])) if n else None,
                status_fractions={s: c/n for s, c in statuses.most_common()} if n else {},
                connected=len(connected),
                connected_fraction=len(connected)/n if n else None,
                involves_source_of_connected=(sum(a['involves_source'] for a in connected)/len(connected)
                                              if connected else None),
                intended_of_connected=sum(a['intended'] for a in connected)/len(connected) if connected else None,
                repeated_fraction=sum(a['repeated'] for a in attempts)/n if n else None,
                mean_evaluations_by_status={s: float(np.mean([a['evaluations'] for a in attempts if a['status'] == s]))
                                            for s in statuses})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaigns', type=Path, nargs='+')
    parser.add_argument('--ids', type=Path, help='JSON list or {"ids": [...]} restricting the starts')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    keep = None
    if args.ids:
        data = json.loads(args.ids.read_text(encoding='utf-8'))
        keep = set(data if isinstance(data, list) else data['ids'])
    pooled, by_origin = defaultdict(list), defaultdict(lambda: defaultdict(list))
    per_run = defaultdict(lambda: defaultdict(list))   # strategy -> start -> per-run dicts
    skipped = []
    for campaign in args.campaigns:
        for task in csv.DictReader((campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'):
            if keep is not None and task['start_id'] not in keep:
                continue
            paths = list((campaign/'runs'/task['outdir']).glob('*/*/network.json'))
            network = json.loads(paths[0].read_text(encoding='utf-8')) if paths else None
            if network is None or network['status'] in ('running', 'aborted_error'):
                skipped.append(f"{campaign.name}:{task['task_id']}")
                continue
            rows = run_outcomes(network)
            strategy = task['strategy']
            pooled[strategy] += rows
            for row in rows:
                by_origin[strategy][row['origin']].append(row)
            n = len(rows)
            per_run[strategy][task['start_id']].append(dict(
                attempts=n, evaluations_per_attempt=np.mean([r['evaluations'] for r in rows]) if n else np.nan,
                connected_fraction=np.mean([r['connected'] for r in rows]) if n else np.nan,
                new_connection_fraction=np.mean([r['status'] == 'new_connection' for r in rows]) if n else np.nan,
                repeated_fraction=np.mean([r['repeated'] for r in rows]) if n else np.nan))
    per_system = {}
    for strategy, systems in per_run.items():
        per_system[strategy] = {start: {k: float(np.nanmean([r[k] for r in runs])) for k in runs[0]}
                                for start, runs in systems.items()}
    result = dict(campaigns=[str(c) for c in args.campaigns], ids=None if args.ids is None else str(args.ids),
                  skipped=skipped,
                  per_strategy={s: block(rows) for s, rows in sorted(pooled.items())},
                  per_origin={s: {o: block(rows) for o, rows in sorted(origins.items())}
                              for s, origins in sorted(by_origin.items())},
                  per_system=per_system)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    for s, b in result['per_strategy'].items():
        top = ', '.join(f'{k} {v:.2f}' for k, v in list(b['status_fractions'].items())[:5])
        print(f"{s:14s} attempts {b['attempts']:5d}  evals/attempt {b['evaluations_per_attempt']:7.1f}  "
              f"connected {b['connected_fraction']:.2f}  intended|conn {b['intended_of_connected'] or 0:.2f}  "
              f"source|conn {b['involves_source_of_connected'] or 0:.2f}  repeated {b['repeated_fraction']:.2f}  [{top}]")
    for s, origins in result['per_origin'].items():
        for o, b in origins.items():
            print(f"  {s:12s} {o:50s} n {b['attempts']:5d}  connected {b['connected_fraction']:.2f}  "
                  f"intended|conn {b['intended_of_connected'] or 0:.2f}  new {b['status_fractions'].get('new_connection', 0):.2f}")


if __name__ == '__main__':
    main()
