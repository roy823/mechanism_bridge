"""Fig. 4: recovery of reference reactions from score_reference_recovery.py outputs.

One column per panel (dataset and track). Top: S1 (mapped) recovery rate per
strategy on the pre-registered denominator, which is the main overlap layers
(unseen + seen_formula) when labels were given, otherwise the representable
reactions, otherwise all rows; error bars are the 95% cluster-bootstrap
intervals over starts from the scorer. Bottom: fraction of rows whose first
mapped hit came within a budget of B evaluations.
Writes FIG.png, FIG.pdf and FIG.csv.
Usage: fig4_recovery.py --panel NAME=SCORES.json [--panel ...] --out FIG [--budget 16000]
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ORDER = ('geometry', 'center_random', 'bond_edits', 'arrows', 'pygsm_edits', 'pygsm_b2f2', 'pygsm_oracle')
COLORS = dict(geometry='#7f7f7f', center_random='#bcbd22', bond_edits='#1f77b4', arrows='#d62728',
              pygsm_edits='#9467bd', pygsm_b2f2='#8c564b', pygsm_oracle='#e377c2')


def denominator(data):
    """(name, rows, scorer summary) of the pre-registered main denominator."""
    rows, subsets = data['rows'], data.get('subsets', {})
    if 'main_layers' in subsets:
        layers = set(data['main_layers'])
        return ('main_layers', [r for r in rows if r['representable'] and r['overlap_label'] in layers],
                subsets['main_layers'])
    if 'representable' in subsets:
        return 'representable', [r for r in rows if r['representable']], subsets['representable']
    return 'all', rows, data['summary']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', action='append', required=True, help='NAME=SCORES.json')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--budget', type=int, default=16000)
    args = parser.parse_args()
    panels = [p.split('=', 1) for p in args.panel]
    fig, axes = plt.subplots(2, len(panels), figsize=(4.4*len(panels), 7), layout='constrained', squeeze=False)
    grid = np.linspace(0, args.budget, 65)
    table = []
    for column, (name, path) in enumerate(panels):
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        layer, rows, summary = denominator(data)
        strategies = [s for s in ORDER if any(r['strategy'] == s for r in rows)]
        top, bottom = axes[0][column], axes[1][column]
        for k, strategy in enumerate(strategies):
            entry = summary[strategy]
            mine = [r for r in rows if r['strategy'] == strategy]
            rate = entry['S1_mapped_rate']
            low, high = entry['S1_mapped_rate_cluster_bootstrap95']
            top.bar(k, rate, color=COLORS[strategy])
            top.errorbar(k, rate, yerr=[[rate - low], [high - rate]], fmt='none', color='black', capsize=4)
            hits = [r['first_mapped_hit']['evaluations_to_hit'] if r.get('first_mapped_hit') else np.inf
                    for r in mine]
            bottom.step(grid, [np.mean([h <= b for h in hits]) for b in grid], where='post',
                        color=COLORS[strategy], label=strategy)
            table.append(dict(panel=name, denominator=layer, strategy=strategy, rows=len(mine),
                              reactions=len({r['reference'] for r in mine}), S1_mapped_rate=float(rate),
                              ci_low=float(low), ci_high=float(high)))
        top.set(xticks=range(len(strategies)), xticklabels=strategies, ylim=(0, 1), ylabel='S1 (mapped) recovery',
                title=f'{name}\n({layer}, {len({r["reference"] for r in rows})} reactions)')
        top.tick_params(axis='x', rotation=30)
        bottom.set(xlabel='MLIP evaluations', ylabel='recovered within budget', ylim=(0, 1))
        bottom.legend(frameon=False, fontsize=8)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix('.png'), dpi=180)
    fig.savefig(args.out.with_suffix('.pdf'))
    plt.close(fig)
    with args.out.with_suffix('.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(table)
    print(json.dumps(table, indent=2))


if __name__ == '__main__':
    main()
