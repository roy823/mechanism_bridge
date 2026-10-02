"""How often each arrow-only seed feature is active (seed_features='arrow_features_v1').

Per strategy, over attempts whose proposal records arrow_terms: lone-pair
direction terms, generalized source/sink links, cyclic push-pull chains,
lone-pair encounter alignment, and the chain progress shift lambda; plus run
statuses, recoverable proposal failures and seed failures for the campaign.
Usage: arrow_feature_usage.py CAMPAIGN_DIR --out FILE
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    manifest = list(csv.DictReader((args.campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'))
    usage, statuses, attempt_statuses, proposal_failures = {}, Counter(), {}, Counter()
    for task in manifest:
        paths = list((args.campaign/'runs'/task['outdir']).glob('*/*/network.json'))
        if not paths:
            statuses['missing'] += 1
            continue
        network = json.loads(paths[0].read_text(encoding='utf-8'))
        strategy = task['strategy']
        statuses[network['status']] += 1
        proposal_failures[strategy] += len(network.get('proposal_failures', []))
        counts = attempt_statuses.setdefault(strategy, Counter())
        row = usage.setdefault(strategy, Counter())
        for attempt in network['attempts']:
            counts[attempt['status']] += 1
            terms = (attempt.get('proposal') or {}).get('arrow_terms')
            row['attempts'] += 1
            if not terms:
                continue
            row['with_arrow_terms'] += 1
            row[f"chain_order={terms.get('chain_order')}"] += 1
            row[f"chain_lambda={terms.get('chain_lambda')}"] += 1
            row['lone_pair_terms>0'] += bool(terms.get('lone_pair_terms'))
            row['generalized_links>0'] += bool(terms.get('generalized_links'))
            row['chain_cyclic'] += bool(terms.get('chain_cyclic'))
            row['encounter_lone_pair_alignment'] += bool(terms.get('encounter_lone_pair_alignment'))
            row['nonzero_chain_shift'] += any(abs(v) > 0 for v in terms.get('chain_shift', []))
    result = dict(campaign=str(args.campaign), run_statuses=dict(statuses),
                  proposal_failures=dict(proposal_failures),
                  attempt_statuses={s: dict(c) for s, c in attempt_statuses.items()},
                  arrow_feature_usage={s: dict(c) for s, c in usage.items()})
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
