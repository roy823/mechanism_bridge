"""Fig. 5b: same net bond edits, different arrows (from summarize_same_edits.py and the runs).

Top row: intended success (source and predicted graph) per arm for each group,
with Wilson 95% intervals; bottom row: TS energies of accepted attempts per arm
(relative to the lowest TS of the group), which shows whether the arrow sets
reach the same or different saddle points.
--overview draws one compact row for many groups instead (no run files needed): accepted
and intended fractions per arm and group; resonance-derived groups are hatched, and the
number of formally charged atoms in the predicted graph is printed when the summary has it.
Usage: fig5b_same_edits.py SUMMARY.json RUN.json [...] --out FIG [--groups START:GROUP ...]
       fig5b_same_edits.py SUMMARY.json --overview --out FIG
"""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reference_matching import wilson  # noqa: E402

ACCEPTED = ('validated_descents', 'validated_core_descents')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('runs', type=Path, nargs='*')
    parser.add_argument('--groups', nargs='+', help='START:GROUP panels to draw (default: all)')
    parser.add_argument('--overview', action='store_true')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding='utf-8'))
    if args.overview:
        return overview(summary, args.out)
    if not args.runs:
        parser.error('RUN.json files are needed unless --overview is given')
    runs = {r['start']: r for r in (json.loads(p.read_text(encoding='utf-8')) for p in args.runs)}
    entries = [e for e in summary if args.groups is None or f"{e['start']}:{e['group']}" in args.groups]
    fig, axes = plt.subplots(2, len(entries), figsize=(3.2*len(entries), 6.4), layout='constrained', squeeze=False)
    for column, entry in enumerate(entries):
        rows = [r for r in runs[entry['start']]['rows'] if r['group'] == entry['group']]
        arms = sorted(entry['per_arm'], key=lambda a: (a == 'bond_edits', a))
        labels = [entry['arrow_sets'][int(a[6:])] if a.startswith('arrows') else 'bond_edits' for a in arms]
        top, bottom = axes[0][column], axes[1][column]
        for k, arm in enumerate(arms):
            stats = entry['per_arm'][arm]
            rate = stats['intended']/stats['attempts']
            low, high = wilson(stats['intended'], stats['attempts'])
            top.bar(k, rate, color='#1f77b4' if arm == 'bond_edits' else '#d62728')
            top.errorbar(k, rate, yerr=[[rate - low], [high - rate]], fmt='none', color='black', capsize=3)
        energies = [r['ts_energy_eV'] for r in rows if r['status'] in ACCEPTED and r.get('ts_energy_eV') is not None]
        base = min(energies) if energies else 0.
        for k, arm in enumerate(arms):
            values = [r['ts_energy_eV'] - base for r in rows
                      if r['arm'] == arm and r['status'] in ACCEPTED and r.get('ts_energy_eV') is not None]
            jitter = np.random.default_rng(k).uniform(-.15, .15, len(values))
            bottom.scatter(np.full(len(values), k) + jitter, values, s=12,
                           color='#1f77b4' if arm == 'bond_edits' else '#d62728', alpha=.7)
        for ax in (top, bottom):
            ax.set_xticks(range(len(arms)), labels, rotation=35, ha='right', fontsize=7)
        top.set(ylim=(0, 1), title=f"{entry['start']}\ngroup {entry['group']}: {entry['predicted_graph']}")
        bottom.set(ylabel='TS energy - lowest (eV)' if column == 0 else None)
        if column == 0:
            top.set_ylabel('intended success')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix('.png'), dpi=180)
    fig.savefig(args.out.with_suffix('.pdf'))
    plt.close(fig)


def overview(summary, out):
    fig, axes = plt.subplots(2, 1, figsize=(max(6., .62*len(summary)), 5.6), layout='constrained', sharex=True)
    width = .27
    for k, entry in enumerate(summary):
        arms = sorted(entry['per_arm'], key=lambda a: (a == 'bond_edits', a))
        for j, arm in enumerate(arms):
            stats = entry['per_arm'][arm]
            x = k + (j - (len(arms) - 1)/2)*width
            color = '#1f77b4' if arm == 'bond_edits' else ('#d62728', '#ff9896')[j % 2]
            hatch = '//' if arm != 'bond_edits' and entry['from_resonance'][int(arm[6:])] else None
            for ax, key in zip(axes, ('accepted', 'intended')):
                ax.bar(x, stats[key]/stats['attempts'], width, color=color, hatch=hatch, edgecolor='white')
        charged = (entry.get('charged_atoms') or {}).get('predicted')
        if charged is not None:
            axes[0].annotate(f'{charged}±', (k, 1.0), ha='center', va='bottom', fontsize=6, color='grey')
    axes[0].set(ylim=(0, 1.08), ylabel='accepted fraction')
    axes[1].set(ylim=(0, 1), ylabel='intended success')
    axes[1].set_xticks(range(len(summary)), [f"{e['start']}:{e['group']}" for e in summary], rotation=60,
                       ha='right', fontsize=7)
    axes[0].set_title('arrow sets (red, hatched = via resonance form) and bond_edits (blue); '
                      'grey: charged atoms in the predicted graph', fontsize=8)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix('.png'), dpi=180)
    fig.savefig(out.with_suffix('.pdf'))
    plt.close(fig)


if __name__ == '__main__':
    main()
