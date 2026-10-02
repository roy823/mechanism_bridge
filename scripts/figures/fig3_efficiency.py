"""Fig. 3: exploration efficiency per strategy from a summarize_campaign.py output.

(a) Anytime curve: distinct root-connected chemical graph pairs (species
connectivity) against cumulative MLIP evaluations. Each run's step curve is
read on a fixed budget grid, averaged over seeds within a system, then over
systems; the band is a 95% percentile bootstrap over systems.
(b) Pairs within the budget: per-system seed means (dots) with the mean and its
bootstrap interval. (c) New species and largest chemical-step depth.
Writes FIG.png, FIG.pdf and FIG.csv (the plotted numbers).
Usage: fig3_efficiency.py SUMMARY.json --out FIG [--grid 500] [--draws 10000]
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ORDER = ('geometry', 'center_random', 'bond_edits', 'arrows')
COLORS = dict(geometry='#7f7f7f', center_random='#bcbd22', bond_edits='#1f77b4', arrows='#d62728')


def value_at(curve, budget):
    found = 0
    for spent, pairs in curve:
        if spent > budget:
            break
        found = pairs
    return found


def per_system(rows, strategy, metric):
    systems = {}
    for row in rows:
        if row['strategy'] == strategy:
            systems.setdefault(row['start'], []).append(metric(row))
    return {s: np.mean(v, axis=0) for s, v in systems.items()}


def bootstrap(values, draws, seed=20261002):
    values = np.asarray(values, dtype=float)
    picks = np.random.default_rng(seed).integers(0, len(values), size=(draws, len(values)))
    means = values[picks].mean(axis=1)
    return np.percentile(means, [2.5, 97.5], axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--grid', type=int, default=500)
    parser.add_argument('--draws', type=int, default=10000)
    args = parser.parse_args()
    data = json.loads(args.summary.read_text(encoding='utf-8'))
    rows, budget = data['rows'], data['budget']
    grid = np.arange(0, budget + 1, args.grid)
    strategies = [s for s in ORDER if any(r['strategy'] == s for r in rows)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout='constrained')
    table = []
    for k, strategy in enumerate(strategies):
        curves = per_system(rows, strategy, lambda r: [value_at(r['anytime'], b) for b in grid])
        stack = np.array(list(curves.values()))
        mean, (low, high) = stack.mean(axis=0), bootstrap(stack, args.draws)
        axes[0].plot(grid, mean, color=COLORS[strategy], label=f'{strategy} (n={len(stack)})')
        axes[0].fill_between(grid, low, high, color=COLORS[strategy], alpha=.15, linewidth=0)
        final = per_system(rows, strategy, lambda r: r['pairs_at_budget'])
        values = np.array(list(final.values()))
        jitter = np.random.default_rng(k).uniform(-.12, .12, len(values))
        axes[1].scatter(np.full(len(values), k) + jitter, values, s=12, color=COLORS[strategy], alpha=.6)
        lo, hi = bootstrap(values, args.draws)
        axes[1].errorbar(k, values.mean(), yerr=[[values.mean() - lo], [hi - values.mean()]], fmt='o',
                         color='black', capsize=4)
        species = np.array(list(per_system(rows, strategy, lambda r: r['new_species']).values()))
        depth = np.array(list(per_system(rows, strategy, lambda r: r['max_species_depth']).values()))
        axes[2].bar(k - .18, species.mean(), width=.36, color=COLORS[strategy])
        axes[2].bar(k + .18, depth.mean(), width=.36, color=COLORS[strategy], alpha=.45, hatch='//')
        table.append(dict(strategy=strategy, systems=len(values), mean_pairs_at_budget=float(values.mean()),
                          ci_low=float(lo), ci_high=float(hi), mean_new_species=float(species.mean()),
                          mean_max_species_depth=float(depth.mean()),
                          anytime=';'.join(f'{b}:{v:.3f}' for b, v in zip(grid, mean))))
    axes[0].set(xlabel='MLIP evaluations', ylabel='root-connected chemical pairs', title='(a) anytime')
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set(xticks=range(len(strategies)), xticklabels=strategies, ylabel=f'pairs within {budget} evaluations',
                title='(b) per-system seed means')
    axes[2].set(xticks=range(len(strategies)), xticklabels=strategies, ylabel='mean over systems',
                title='(c) new species (solid) and depth (hatched)')
    for ax in axes[1:]:
        ax.tick_params(axis='x', rotation=20)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix('.png'), dpi=180)
    fig.savefig(args.out.with_suffix('.pdf'))
    plt.close(fig)
    with args.out.with_suffix('.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(table)
    print(json.dumps([{k: v for k, v in r.items() if k != 'anytime'} for r in table], indent=2))


if __name__ == '__main__':
    main()
