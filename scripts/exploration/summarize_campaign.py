"""Summarize a campaign (Fig. 3 efficiency / Fig. 5a information ablation).

Per run: final status, evaluations, root-connected distinct chemical graph
pairs, new species, longest enumerated multistep path, and the anytime curve
(distinct root-connected chemical pairs versus cumulative evaluations, attempt
order; serial runs). Per strategy: means over systems of per-system seed means,
and paired Wilcoxon signed-rank tests across systems against a reference
strategy (e.g. arrows vs bond_edits, every strategy vs geometry).
Usage: summarize_campaign.py CAMPAIGN_DIR --out FILE [--budget 16000]
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.network_metrics import growth_metrics  # noqa: E402


def root_pairs(network, edges):
    """Distinct unordered species pairs of chemical edges in the root component."""
    graph = nx.Graph()
    graph.add_nodes_from(range(len(network['nodes'])))
    graph.add_edges_from(e['nodes'] for e in edges)
    if not network['nodes']:
        return set()
    connected = nx.node_connected_component(graph, 0)
    smiles = [n['graph_smiles'] for n in network['nodes']]
    return {tuple(sorted(smiles[i] for i in e['nodes'])) for e in edges
            if e['kind'] == 'chemical' and set(e['nodes']) <= connected}


def anytime(network):
    """[(cumulative evaluations, distinct root-connected chemical pairs)] after each attempt."""
    spent, curve = network.get('initialization_evaluations', 0), []
    for k, attempt in enumerate(network['attempts']):
        spent += attempt['evaluations']
        curve.append((spent, len(root_pairs(network, [e for e in network['edges'] if e['attempt'] <= k]))))
    return curve


def value_at(curve, budget):
    """Pairs found within a budget (step function of the anytime curve)."""
    found = 0
    for spent, pairs in curve:
        if spent > budget:
            break
        found = pairs
    return found


def paired_test(per_system, strategy, reference):
    common = sorted(set(per_system.get(strategy, {})) & set(per_system.get(reference, {})))
    x = np.array([per_system[strategy][s] for s in common])
    y = np.array([per_system[reference][s] for s in common])
    if len(common) < 2 or np.allclose(x, y):
        return dict(systems=len(common), mean_difference=float(np.mean(x-y)) if common else None, p=None)
    return dict(systems=len(common), mean_difference=float(np.mean(x-y)),
                wins=int(np.sum(x > y)), losses=int(np.sum(x < y)),
                p=float(wilcoxon(x, y, zero_method='zsplit').pvalue))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--budget', type=int, default=16000)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    manifest = list(csv.DictReader((args.campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t'))
    rows, missing = [], []
    for task in manifest:
        paths = list((args.campaign/'runs'/task['outdir']).glob('*/*/network.json'))
        if not paths:
            missing.append(task['task_id'])
            continue
        network = json.loads(paths[0].read_text(encoding='utf-8'))
        metrics = growth_metrics(network)
        curve = anytime(network)
        rows.append(dict(task=int(task['task_id']), start=task['start_id'], strategy=task['strategy'],
                         seed=int(task['seed']), status=network['status'], evaluations=network.get('evaluations'),
                         attempts=len(network['attempts']), edges=len(network['edges']),
                         root_chemical_pairs=len(root_pairs(network, network['edges'])),
                         pairs_at_budget=value_at(curve, args.budget),
                         new_species=len(metrics['root_species']) - 1,
                         longest_enumerated_chemical_path=metrics['longest_enumerated_chemical_path'],
                         anytime=curve))
    per_system = {}
    for row in rows:
        per_system.setdefault(row['strategy'], {}).setdefault(row['start'], []).append(row['pairs_at_budget'])
    per_system = {s: {k: float(np.mean(v)) for k, v in d.items()} for s, d in per_system.items()}
    strategies = sorted(per_system)
    summary = {s: dict(systems=len(per_system[s]), mean_pairs_at_budget=float(np.mean(list(per_system[s].values()))))
               for s in strategies}
    tests = {}
    for s in strategies:
        if s != 'geometry' and 'geometry' in per_system:
            tests[f'{s}_vs_geometry'] = paired_test(per_system, s, 'geometry')
    if 'arrows' in per_system and 'bond_edits' in per_system:
        tests['arrows_vs_bond_edits'] = paired_test(per_system, 'arrows', 'bond_edits')
    statuses = {}
    for row in rows:
        statuses[row['status']] = statuses.get(row['status'], 0) + 1
    result = dict(campaign=str(args.campaign), budget=args.budget, tasks=len(manifest), finished=len(rows),
                  missing_tasks=missing, run_statuses=statuses, per_strategy=summary, paired_tests=tests,
                  per_system_means=per_system, rows=rows)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'per_system_means')}, indent=2))


if __name__ == '__main__':
    main()
