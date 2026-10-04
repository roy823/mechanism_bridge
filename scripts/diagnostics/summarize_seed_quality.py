"""Summarize seed_quality.py outputs: how close each arm's seeds come to the reference TS.

Per start and arm: the median core RMSD and edit error over matching proposals, samples and
random seeds. Over starts: the median of those values, and paired comparisons with bond_edits and
with full arrows (arrows:none): mean difference, wins (arm closer), losses and the Wilcoxon
signed-rank p (zsplit). Starts are split by whether a lone-pair term was active in the full-arrows
seeds (lone_pair_terms > 0), since those terms are what the reversed and shuffled arms switch off.
Usage: summarize_seed_quality.py DIR [...] --out FILE
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

METRICS = ('core_rmsd_A', 'edit_error_A')


def paired(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or np.allclose(a, b):
        return dict(n=len(a), mean_difference=float(np.mean(a - b)) if len(a) else None, closer=0, farther=0, p=None)
    return dict(n=len(a), mean_difference=float(np.mean(a - b)), closer=int(np.sum(a < b)),
                farther=int(np.sum(a > b)), p=float(wilcoxon(a, b, zero_method='zsplit').pvalue))


def block(per_start, starts):
    arms = sorted({arm for s in starts for arm in per_start[s]})
    out = dict(starts=len(starts), arms={})
    for arm in arms:
        have = [s for s in starts if arm in per_start[s]]
        entry = dict(starts=len(have))
        for m in METRICS:
            entry[f'median_{m}'] = float(np.median([per_start[s][arm][m] for s in have])) if have else None
            for ref in ('bond_edits', 'arrows:none'):
                if arm == ref:
                    continue
                common = [s for s in have if ref in per_start[s]]
                entry[f'{m}_vs_{ref}'] = paired([per_start[s][arm][m] for s in common],
                                                [per_start[s][ref][m] for s in common])
        out['arms'][arm] = entry
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dirs', type=Path, nargs='+')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    per_start, lone_pair, skipped = {}, {}, []
    for directory in args.dirs:
        for path in sorted(directory.glob('*.json')):
            data = json.loads(path.read_text(encoding='utf-8'))
            if data.get('skipped') or not data.get('rows'):
                skipped.append(dict(start=data['start'], reason=data.get('skipped') or 'no matching proposal'))
                continue
            values = defaultdict(lambda: defaultdict(list))
            for row in data['rows']:
                if 'core_rmsd_A' in row:
                    for m in METRICS:
                        values[row['arm']][m].append(row[m])
                    if row['arm'] == 'arrows:none':
                        values['_lp']['n'].append(row.get('lone_pair_terms') or 0)
            per_start[data['start']] = {arm: {m: float(np.median(v[m])) for m in METRICS}
                                        for arm, v in values.items() if arm != '_lp'}
            lone_pair[data['start']] = bool(values['_lp']['n']) and max(values['_lp']['n']) > 0
    starts = sorted(per_start)
    result = dict(dirs=[str(d) for d in args.dirs], skipped=skipped,
                  all=block(per_start, starts),
                  lone_pair_active=block(per_start, [s for s in starts if lone_pair[s]]),
                  lone_pair_inactive=block(per_start, [s for s in starts if not lone_pair[s]]),
                  per_start=per_start, lone_pair=lone_pair)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    for name in ('all', 'lone_pair_active', 'lone_pair_inactive'):
        b = result[name]
        print(f"== {name}: {b['starts']} starts")
        for arm, e in b['arms'].items():
            vb, va = e.get('core_rmsd_A_vs_bond_edits', {}), e.get('core_rmsd_A_vs_arrows:none', {})
            print(f"  {arm:32s} core RMSD {e['median_core_rmsd_A']:.3f}  edit err {e['median_edit_error_A']:.3f}"
                  f"  vs bond_edits {vb.get('mean_difference') or 0:+.3f} ({vb.get('closer', '-')}/{vb.get('farther', '-')}, p={vb.get('p')})"
                  f"  vs arrows {va.get('mean_difference') or 0:+.3f} ({va.get('closer', '-')}/{va.get('farther', '-')}, p={va.get('p')})")


if __name__ == '__main__':
    main()
