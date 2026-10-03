"""Fig. 5c: seed-feature ablations of the arrows strategy (from summarize_campaign.py).

Input: one summarize_campaign.py output that pools the main campaigns (arrows and
bond_edits runs) with the ablation campaigns (labels arrows@FP-JCTC-1-ablation-<arm>),
restricted with --ids to the Fig. 5c systems. Only systems with runs in every arm count.
Per arm: root-connected species pairs within the budget as per-system seed means, the
mean over systems and a 95% percentile bootstrap interval over systems (seed 20261002).
Tests (pre-registered in docs/10-02 JCTC主战役执行计划.md §2.3): each ablation versus the
full arrows arm by the paired Wilcoxon signed-rank test over systems (zero_method zsplit,
as in summarize_campaign.py), Holm-corrected over the six ablations; every arm versus
bond_edits is reported uncorrected. --outcomes adds the attempt-level panel from
attempt_outcomes.py (connected fraction and TS-force non-convergence per arm).
Writes FIG.png, FIG.pdf and FIG.json.
Usage: fig5c_ablation.py SUMMARY.json --out FIG [--outcomes OUTCOMES.json] [--draws 10000]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

ABLATIONS = ('no_lone_pair_direction', 'no_chain_order', 'no_generalized_links', 'no_encounter_alignment',
             'reversed_arrows', 'shuffled_arrows')
ARMS = ('arrows',) + tuple(f'arrows@FP-JCTC-1-ablation-{a}' for a in ABLATIONS) + ('bond_edits',)


def short(arm):
    return arm.split('ablation-')[-1]


def paired(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 2 or np.allclose(x, y):
        return dict(mean_difference=float(np.mean(x - y)) if len(x) else None, wins=int(np.sum(x > y)),
                    losses=int(np.sum(x < y)), p=None)
    return dict(mean_difference=float(np.mean(x - y)), wins=int(np.sum(x > y)), losses=int(np.sum(x < y)),
                p=float(wilcoxon(x, y, zero_method='zsplit').pvalue))


def holm(ps):
    order = sorted((p, k) for k, p in enumerate(ps) if p is not None)
    adjusted, running = [None]*len(ps), 0.
    for rank, (p, k) in enumerate(order):
        running = max(running, min(1., (len(order) - rank)*p))
        adjusted[k] = running
    return adjusted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--outcomes', type=Path)
    parser.add_argument('--draws', type=int, default=10000)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding='utf-8'))
    values = {}
    for row in summary['rows']:
        values.setdefault(row['strategy'], {}).setdefault(row['start'], []).append(row['pairs_at_budget'])
    arms = [a for a in ARMS if a in values]
    systems = sorted(set.intersection(*(set(values[a]) for a in arms)))
    means = {a: np.array([np.mean(values[a][s]) for s in systems]) for a in arms}
    rng = np.random.default_rng(20261002)
    picks = rng.integers(0, len(systems), size=(args.draws, len(systems)))
    per_arm = {a: dict(systems=len(systems), seeds_per_system=sorted({len(values[a][s]) for s in systems}),
                       mean=float(m.mean()), ci95=[float(v) for v in np.percentile(m[picks].mean(1), [2.5, 97.5])])
               for a, m in means.items()}
    ablations = [a for a in arms if '@' in a]
    versus_arrows = {a: paired(means[a], means['arrows']) for a in ablations} if 'arrows' in means else {}
    for a, p in zip(versus_arrows, holm([v['p'] for v in versus_arrows.values()])):
        versus_arrows[a]['p_holm'] = p
    versus_bond_edits = ({a: paired(means[a], means['bond_edits']) for a in arms if a != 'bond_edits'}
                         if 'bond_edits' in means else {})
    outcomes = json.loads(args.outcomes.read_text(encoding='utf-8'))['per_strategy'] if args.outcomes else {}
    result = dict(summary=str(args.summary), budget=summary.get('budget'), systems=systems, per_arm=per_arm,
                  versus_arrows=versus_arrows, versus_bond_edits=versus_bond_edits,
                  per_system={a: dict(zip(systems, m.tolist())) for a, m in means.items()},
                  attempt_outcomes={a: {k: outcomes[a][k] for k in ('attempts', 'evaluations_per_attempt',
                                                                    'connected_fraction')}
                                    | dict(ts_force_unconverged=outcomes[a]['status_fractions'].get(
                                           'ts_force_unconverged', 0.))
                                    for a in arms if a in outcomes})
    panels = 3 if outcomes else 2
    fig, axes = plt.subplots(1, panels, figsize=(5.2*panels, 4.4), layout='constrained')
    x = np.arange(len(arms))
    colors = {a: '#1f77b4' if a == 'bond_edits' else '#ff9896' if '@' in a else '#d62728' for a in arms}
    jitter = np.random.default_rng(0).uniform(-.15, .15, len(systems))
    for k, a in enumerate(arms):
        low, high = per_arm[a]['ci95']
        axes[0].bar(k, per_arm[a]['mean'], color=colors[a])
        axes[0].errorbar(k, per_arm[a]['mean'], yerr=[[per_arm[a]['mean'] - low], [high - per_arm[a]['mean']]],
                         fmt='none', color='black', capsize=3)
        axes[0].scatter(k + jitter, means[a], s=10, color='black', alpha=.5, zorder=3)
    axes[0].set(xticks=x, ylabel=f"root-connected pairs within {summary.get('budget')} evaluations",
                title='(a) per-system seed means')
    axes[0].set_xticklabels([short(a) for a in arms], rotation=40, ha='right', fontsize=8)
    for k, a in enumerate(ablations):
        diff = means[a] - means['arrows']
        axes[1].scatter(np.full(len(diff), k) + jitter, diff, s=12, color='#ff9896')
        axes[1].scatter(k, diff.mean(), s=40, color='black', marker='D')
        p = versus_arrows[a]['p_holm']
        axes[1].annotate('n.d.' if p is None else f'p={p:.2g}', (k, max(diff.max(), 0)), textcoords='offset points',
                         xytext=(0, 6), ha='center', fontsize=7)
    axes[1].axhline(0, color='grey', lw=.8)
    axes[1].set(xticks=range(len(ablations)), ylabel='ablation - full arrows (pairs)',
                title='(b) paired differences (Holm p)')
    axes[1].set_xticklabels([short(a) for a in ablations], rotation=40, ha='right', fontsize=8)
    if outcomes:
        shown = [a for a in arms if a in result['attempt_outcomes']]
        stats = result['attempt_outcomes']
        axes[2].bar(np.arange(len(shown)) - .2, [stats[a]['connected_fraction'] for a in shown], .4, label='connected')
        axes[2].bar(np.arange(len(shown)) + .2, [stats[a]['ts_force_unconverged'] for a in shown], .4,
                    label='TS force unconverged')
        axes[2].set(xticks=range(len(shown)), ylabel='fraction of attempts', title='(c) attempt outcomes')
        axes[2].set_xticklabels([short(a) for a in shown], rotation=40, ha='right', fontsize=8)
        axes[2].legend(fontsize=7)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix('.png'), dpi=180)
    fig.savefig(args.out.with_suffix('.pdf'))
    plt.close(fig)
    args.out.with_suffix('.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('per_system',)}, indent=1))


if __name__ == '__main__':
    main()
