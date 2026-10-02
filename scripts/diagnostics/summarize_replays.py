"""Summarize P1/P2 replays against the pre-registered decision rules (E doc, section 1.3).

P1 adopts Dimer->Sella if either
  (1) index-one TS rate >= A - 2 points and median evaluations per found TS <= 0.85 x A, or
  (2) index-one TS rate >= A + 5 points and median evaluations per found TS <= 1.10 x A.
P2 adopts MLIP IRC in search if its accepted-connection endpoint graphs agree with
the displacement protocol at least as often as they are accepted, and its evaluations
per attempt are <= 2 x displacement; otherwise IRC is used only on reported edges.
Same TS: aligned RMSD <= 0.15 A and |dE| <= 0.03 eV (registry tolerances).
Usage: summarize_replays.py RESULTS_DIR --out FILE
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reaction_network import aligned_rmsd  # noqa: E402
from mechbridge.reference_matching import wilson  # noqa: E402

INDEX_ONE = {'validated_descents', 'validated_core_descents', 'unresolved_minimum'}
ACCEPTED = {'validated_descents', 'validated_core_descents'}


def load(results_dirs):
    rows = []
    for results in results_dirs:
        for path in sorted(results.glob('task_*.json')):
            for row in json.loads(path.read_text(encoding='utf-8')):
                rows.append(dict(row, results_dir=str(results)))
    return rows


def ts_geometry(row):
    path = (Path(row['results_dir'])/'attempts'/row['group']/row['attempt'].replace('/', '__')/
            row['variant']/'ts.xyz')
    return read(path).positions if path.exists() else None


def rate(rows, statuses):
    hits = sum(r['status'] in statuses for r in rows)
    return dict(hits=hits, n=len(rows), rate=hits/len(rows) if rows else None,
                wilson95=wilson(hits, len(rows)))


def median_evaluations(rows):
    found = [r['evaluations'] for r in rows if r['status'] in INDEX_ONE]
    return float(np.median(found)) if found else None


def p1_summary(rows):
    by = {}
    for row in rows:
        by.setdefault(row['attempt'], {})[row['variant']] = row
    pairs = [v for v in by.values() if {'legacy', 'sella'} <= set(v)]
    legacy, sella = [p['legacy'] for p in pairs], [p['sella'] for p in pairs]
    same = both = 0
    for pair in pairs:
        a, b = pair['legacy'], pair['sella']
        if a['status'] in INDEX_ONE and b['status'] in INDEX_ONE:
            both += 1
            xa, xb = ts_geometry(a), ts_geometry(b)
            if (xa is not None and xb is not None and aligned_rmsd(xa, xb) <= .15
                    and abs(a['ts_energy_eV'] - b['ts_energy_eV']) <= .03):
                same += 1
    rate_a, rate_b = rate(legacy, INDEX_ONE), rate(sella, INDEX_ONE)
    med_a, med_b = median_evaluations(legacy), median_evaluations(sella)
    delta = None if None in (rate_a['rate'], rate_b['rate']) else 100*(rate_b['rate'] - rate_a['rate'])
    ratio = None if None in (med_a, med_b) or not med_a else med_b/med_a
    adopt = bool(delta is not None and ratio is not None and
                 ((delta >= -2 and ratio <= .85) or (delta >= 5 and ratio <= 1.10)))
    statuses = {v: {} for v in ('legacy', 'sella')}
    for variant, group in (('legacy', legacy), ('sella', sella)):
        for r in group:
            statuses[variant][r['status']] = statuses[variant].get(r['status'], 0) + 1
    rescued = sum(p['legacy']['status'] == 'ts_force_unconverged' and p['sella']['status'] in INDEX_ONE
                  for p in pairs)
    return dict(attempts=len(pairs), index_one_rate=dict(legacy=rate_a, sella=rate_b),
                accepted_rate=dict(legacy=rate(legacy, ACCEPTED), sella=rate(sella, ACCEPTED)),
                median_evaluations_per_found_TS=dict(legacy=med_a, sella=med_b),
                total_evaluations=dict(legacy=sum(r['evaluations'] for r in legacy),
                                       sella=sum(r['evaluations'] for r in sella)),
                rate_difference_points=delta, evaluation_ratio=ratio,
                both_found=both, same_TS_when_both_found=same,
                rescued_force_unconverged=rescued, statuses=statuses,
                legacy_historical_match=sum(bool(r.get('historical_match')) for r in legacy),
                legacy_historical_status_match=sum(r['status'] == r['historical_status'] for r in legacy),
                adopt_dimer_sella=adopt,
                sella_irc=combined_summary(by))


def combined_summary(by):
    """Dimer->Sella with and without MLIP IRC endpoints on the same seeds (P1b)."""
    pairs = [v for v in by.values() if {'sella', 'sella_irc'} <= set(v)]
    if not pairs:
        return None
    sella, combined = [p['sella'] for p in pairs], [p['sella_irc'] for p in pairs]
    return dict(attempts=len(pairs), accepted_rate=dict(sella=rate(sella, ACCEPTED),
                                                        sella_irc=rate(combined, ACCEPTED)),
                index_one_rate=dict(sella=rate(sella, INDEX_ONE), sella_irc=rate(combined, INDEX_ONE)),
                unresolved_rescued=sum(p['sella']['status'] == 'unresolved_minimum' and
                                       p['sella_irc']['status'] in ACCEPTED for p in pairs),
                accepted_lost=sum(p['sella']['status'] in ACCEPTED and
                                  p['sella_irc']['status'] not in ACCEPTED for p in pairs),
                median_evaluations=dict(sella=float(np.median([r['evaluations'] for r in sella])),
                                        sella_irc=float(np.median([r['evaluations'] for r in combined]))))


def p2_summary(rows):
    by = {}
    for row in rows:
        by.setdefault(row['attempt'], {})[row['variant']] = row
    pairs = [v for v in by.values() if {'legacy', 'irc'} <= set(v)]
    legacy, irc = [p['legacy'] for p in pairs], [p['irc'] for p in pairs]
    both = [p for p in pairs if p['legacy']['status'] in ACCEPTED and p['irc']['status'] in ACCEPTED]
    agree = sum(sorted(p['legacy']['endpoint_graphs']) == sorted(p['irc']['endpoint_graphs']) for p in both)
    eval_ratio = (np.median([p['irc']['evaluations']/p['legacy']['evaluations'] for p in pairs
                             if p['legacy']['evaluations']]) if pairs else None)
    accepted_legacy, accepted_irc = rate(legacy, ACCEPTED), rate(irc, ACCEPTED)
    adopt = bool(eval_ratio is not None and eval_ratio <= 2 and
                 accepted_irc['hits'] >= accepted_legacy['hits'])
    return dict(attempts=len(pairs), accepted_rate=dict(legacy=accepted_legacy, irc=accepted_irc),
                irc_converged_both_branches=sum(bool(r['is_IRC']) for r in irc),
                both_accepted=len(both), endpoint_graphs_agree=agree,
                median_evaluation_ratio_irc_over_legacy=None if eval_ratio is None else float(eval_ratio),
                legacy_historical_match=sum(bool(r.get('historical_match')) for r in legacy),
                adopt_irc_in_search=adopt,
                caveat='No DFT IRC references in this sample; agreement is with the displacement protocol')


def p3_summary(rows):
    """fmax sensitivity on the P2 seeds: TS and endpoint identity versus fmax 0.005."""
    by = {}
    for row in rows:
        by.setdefault(row['attempt'], {})[row['variant']] = row
    out = {}
    for other in ('sella_irc_f003', 'sella_irc_f010'):
        pairs = [v for v in by.values() if {'sella_irc', other} <= set(v)]
        if not pairs:
            continue
        both = [p for p in pairs if p['sella_irc']['status'] in ACCEPTED and p[other]['status'] in ACCEPTED]
        same_ts = 0
        for p in both:
            xa, xb = ts_geometry(p['sella_irc']), ts_geometry(p[other])
            if (xa is not None and xb is not None and aligned_rmsd(xa, xb) <= .15 and
                    abs(p['sella_irc']['ts_energy_eV'] - p[other]['ts_energy_eV']) <= .03):
                same_ts += 1
        same_ends = sum(sorted(p['sella_irc']['endpoint_graphs']) == sorted(p[other]['endpoint_graphs'])
                        for p in both)
        out[other] = dict(pairs=len(pairs), accepted=dict(fmax_005=sum(p['sella_irc']['status'] in ACCEPTED
                                                                         for p in pairs),
                                                          other=sum(p[other]['status'] in ACCEPTED for p in pairs)),
                          both_accepted=len(both), same_TS=same_ts, same_endpoints=same_ends,
                          median_evaluation_ratio=float(np.median([p[other]['evaluations']/p['sella_irc']['evaluations']
                                                                   for p in pairs if p['sella_irc']['evaluations']])))
    keep = out.get('sella_irc_f003')
    out['keep_fmax_005'] = bool(keep and keep['both_accepted'] and
                                min(keep['same_TS'], keep['same_endpoints']) >= .98*keep['both_accepted'])
    return out


def repeat_consistency(rows):
    """Same (group, attempt, variant) replayed in several result directories."""
    seen, repeats, identical = {}, 0, 0
    for row in rows:
        key = (row['group'], row['attempt'], row['variant'])
        if key in seen:
            repeats += 1
            first = seen[key]
            identical += first['status'] == row['status'] and first['evaluations'] == row['evaluations']
        else:
            seen[key] = row
    return dict(repeated_replays=repeats, identical_status_and_evaluations=identical)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path, nargs='+')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    rows = load(args.results)
    summary = dict(results=[str(r) for r in args.results], rows=len(rows),
                   P1=p1_summary([r for r in rows if r['group'] == 'P1']),
                   P2=p2_summary([r for r in rows if r['group'] == 'P2']),
                   P3=p3_summary([r for r in rows if r['group'] == 'P2']),
                   repeats=repeat_consistency(rows))
    args.out.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
