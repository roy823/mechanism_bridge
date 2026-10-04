"""Paired attempts: same proposal and random numbers, arrows versus bond_edits seed features.

For one start and run seed, the arrows and bond_edits runs share the root, the proposal schedule
and the attempt RNG (reaction_network: the seed RNG follows node geometry, edits, run seed and
variant, not the strategy), so their attempts are identical apart from the seed features until the
networks diverge. This pairs attempts k = 0, 1, ... while source graph, edits and proposal random
seed agree, and compares outcomes:
  connected : IRC joined two minima (new_connection or duplicate_connection);
  intended  : connected, and the two ends are the source graph and the predicted graph;
  status    : the full attempt status, as a transition table.
Exact McNemar tests on connected and on intended; split by whether a lone-pair term was active in
the arrows seed (proposal arrow_terms.lone_pair_terms > 0). Reads networks only.
Usage: paired_attempts.py CAMPAIGN_DIR [...] --out FILE [--a arrows --b bond_edits]
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reference_matching import stereo_free  # noqa: E402

CONNECTED = ('new_connection', 'duplicate_connection')


def plain(smiles):
    try:
        return stereo_free(smiles)
    except (ValueError, TypeError):
        return smiles


def network(campaign, outdir):
    paths = list((campaign/'runs'/outdir).glob('*/*/network.json'))
    if not paths:
        return None
    data = json.loads(paths[0].read_text(encoding='utf-8'))
    return None if data['status'] in ('running', 'aborted_error') else data


def outcome(net, attempt):
    graphs = [plain(n['graph_smiles']) for n in net['nodes']]
    proposal = attempt.get('proposal') or {}
    ends = [graphs[i] for i in attempt.get('observed_nodes') or []]
    connected = attempt['status'] in CONNECTED and len(ends) == 2
    source = graphs[attempt['source_node']]
    predicted = plain(proposal['predicted_graph']) if proposal.get('predicted_graph') else None
    return dict(status=attempt['status'], connected=connected, evaluations=attempt['evaluations'],
                intended=connected and predicted is not None and sorted(ends) == sorted([source, predicted]))


def key(net, attempt):
    proposal = attempt.get('proposal') or {}
    return (net['nodes'][attempt['source_node']]['graph_smiles'], json.dumps(proposal.get('edits'), sort_keys=True),
            proposal.get('random_seed'))


def mcnemar(pairs, field):
    a_only = sum(p['a'][field] and not p['b'][field] for p in pairs)
    b_only = sum(p['b'][field] and not p['a'][field] for p in pairs)
    return dict(pairs=len(pairs), a=sum(p['a'][field] for p in pairs), b=sum(p['b'][field] for p in pairs),
                a_only=a_only, b_only=b_only,
                p=binomtest(a_only, a_only + b_only, .5).pvalue if a_only + b_only else None)


def block(pairs):
    return dict(pairs=len(pairs), connected=mcnemar(pairs, 'connected'), intended=mcnemar(pairs, 'intended'),
                evaluations_a=sum(p['a']['evaluations'] for p in pairs),
                evaluations_b=sum(p['b']['evaluations'] for p in pairs),
                transitions=Counter(f"{p['a']['status']} | {p['b']['status']}" for p in pairs).most_common(12))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaigns', type=Path, nargs='+')
    parser.add_argument('--a', default='arrows')
    parser.add_argument('--b', default='bond_edits')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = dict(a=args.a, b=args.b, campaigns={})
    everything = []
    for campaign in args.campaigns:
        tasks = list(csv.DictReader((campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'))
        runs = {(t['start_id'], t['seed'], t['strategy']): t['outdir'] for t in tasks}
        pairs, prefix_lengths = [], []
        for (start, seed, strategy), outdir in sorted(runs.items()):
            if strategy != args.a or (start, seed, args.b) not in runs:
                continue
            net_a, net_b = network(campaign, outdir), network(campaign, runs[start, seed, args.b])
            if net_a is None or net_b is None:
                continue
            k = 0
            for att_a, att_b in zip(net_a['attempts'], net_b['attempts']):
                if key(net_a, att_a) != key(net_b, att_b):
                    break
                lp = ((att_a.get('proposal') or {}).get('arrow_terms') or {}).get('lone_pair_terms') or 0
                pairs.append(dict(start=start, seed=seed, k=k, lone_pair=lp > 0,
                                  a=outcome(net_a, att_a), b=outcome(net_b, att_b)))
                k += 1
            prefix_lengths.append(k)
        out['campaigns'][campaign.name] = dict(
            run_pairs=len(prefix_lengths), mean_matched_prefix=sum(prefix_lengths)/max(1, len(prefix_lengths)),
            all=block(pairs), lone_pair_active=block([p for p in pairs if p['lone_pair']]),
            lone_pair_inactive=block([p for p in pairs if not p['lone_pair']]))
        everything += pairs
    out['pooled'] = dict(all=block(everything), lone_pair_active=block([p for p in everything if p['lone_pair']]),
                         lone_pair_inactive=block([p for p in everything if not p['lone_pair']]))
    args.out.write_text(json.dumps(out, indent=2), encoding='utf-8')
    for name, c in list(out['campaigns'].items()) + [('POOLED', out['pooled'])]:
        print(f"== {name}" + (f": {c['run_pairs']} run pairs, mean matched prefix {c['mean_matched_prefix']:.1f}"
                               if 'run_pairs' in c else ''))
        for part in ('all', 'lone_pair_active', 'lone_pair_inactive'):
            b = c[part]
            cn, it = b['connected'], b['intended']
            print(f"   {part:20s} pairs {b['pairs']:5d}  connected {args.a} {cn['a']} vs {args.b} {cn['b']}"
                  f" ({cn['a_only']}/{cn['b_only']}, p={cn['p']})  intended {it['a']} vs {it['b']}"
                  f" ({it['a_only']}/{it['b_only']}, p={it['p']})")


if __name__ == '__main__':
    main()
