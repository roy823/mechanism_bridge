"""Freeze the representable-reaction list from a representability screen.

Reads every task_*.json under the screen directory, reports the funnel
(P-RFO -> index-one -> IRC endpoints matching the reference) and writes the
sorted IDs of representable reactions; these define the main Fig. 4
denominator. The list is written once and never edited by hand.
Usage: summarize_screen.py SCREEN_DIR --out FILE
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('screen', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    tasks = sorted(args.screen.glob('task_*.json'))
    rows = [row for path in tasks for row in json.loads(path.read_text(encoding='utf-8'))['rows']]
    meta = json.loads(tasks[0].read_text(encoding='utf-8'))
    ids = sorted(r['id'] for r in rows if r['representable'])
    stages = {}
    for row in rows:
        stages[row['stage']] = stages.get(row['stage'], 0) + 1
    screened = [r for r in rows if r['stage'] == 'screened']
    rmsd = [r['ts_rmsd_to_reference_A'] for r in rows if r.get('ts_rmsd_to_reference_A') is not None]
    summary = dict(screen=str(args.screen), tasks=len(tasks), reactions=len(rows), stages=stages,
                   index_one=len(screened), representable=len(ids),
                   unmapped_match=sum(bool(r.get('unmapped_match')) for r in screened),
                   mapped_match=sum(bool(r.get('mapped_match')) for r in screened),
                   both_irc_branches_converged=sum(all(r.get('irc_converged', [])) for r in screened),
                   ts_rmsd_to_reference_A=dict(median=float(np.median(rmsd)) if rmsd else None,
                                               p90=float(np.percentile(rmsd, 90)) if rmsd else None),
                   evaluations=dict(median=float(np.median([r['evaluations'] for r in rows])),
                                    total=int(sum(r['evaluations'] for r in rows))),
                   potential=meta['potential'], protocol=meta['protocol'],
                   representable_ids=ids,
                   representable_ids_sha256=hashlib.sha256('\n'.join(ids).encode()).hexdigest())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('representable_ids', 'potential')},
                     indent=2))


if __name__ == '__main__':
    main()
