"""Summarize a campaign (Fig. 3 efficiency / Fig. 5a information ablation).

Per run: final status, evaluations, root-connected distinct chemical graph
pairs, new species, longest enumerated multistep path, and the anytime curve
(distinct root-connected chemical pairs versus cumulative evaluations, attempt
order; serial runs). Primary connectivity is by species: minima with one graph
are one species, so a reaction from a re-oriented encounter complex of the root
species counts (pilot 3: every formaldehyde-dimer edge started from such a
complex and none touched the root minimum). The physical (node-level) counts
are reported alongside. Per strategy: means over systems of per-system seed means,
and paired Wilcoxon signed-rank tests across systems against a reference
strategy (e.g. arrows vs bond_edits, every strategy vs geometry).
Several campaigns can be pooled (e.g. Fig. 3 = the first 20 Coley starts of the
Coley campaign plus the growth campaign) and restricted with --ids. A campaign
whose protocol is not FP-JCTC-1 labels its runs strategy@protocol (e.g. the
arrows-legacy arm of Fig. 5a, the seed-feature ablations of Fig. 5c), and each
such arm is also tested against its base strategy. Runs that are still running
or ended in aborted_error count as missing, not as zero.
Usage: summarize_campaign.py CAMPAIGN_DIR [...] --out FILE [--budget 16000] [--ids F.json]
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


def root_pairs(network, edges, species=True):
    """Distinct unordered species pairs of chemical edges connected to the root.

    species=True: connectivity over species (graph SMILES); False: over minima.
    """
    nodes = network['nodes']
    if not nodes:
        return set()
    smiles = [n['graph_smiles'] for n in nodes]
    graph = nx.Graph()
    if species:
        graph.add_nodes_from(smiles)
        graph.add_edges_from((smiles[a], smiles[b]) for a, b in (e['nodes'] for e in edges))
        connected = nx.node_connected_component(graph, smiles[0])
        inside = [smiles[e['nodes'][0]] in connected for e in edges]
    else:
        graph.add_nodes_from(range(len(nodes)))
        graph.add_edges_from(e['nodes'] for e in edges)
        connected = nx.node_connected_component(graph, 0)
        inside = [set(e['nodes']) <= connected for e in edges]
    return {tuple(sorted(smiles[i] for i in e['nodes'])) for e, ok in zip(edges, inside)
            if ok and e['kind'] == 'chemical'}


def species_reach(network):
    """(new species joined to the root species by chemical edges, largest chemical-step depth)."""
    nodes = network['nodes']
    if not nodes:
        return 0, 0
    smiles = [n['graph_smiles'] for n in nodes]
    graph = nx.Graph()
    graph.add_node(smiles[0])
    graph.add_edges_from((smiles[a], smiles[b]) for e in network['edges'] if e['kind'] == 'chemical'
                         for a, b in [e['nodes']] if smiles[a] != smiles[b])
    depth = nx.single_source_shortest_path_length(graph, smiles[0])
    return len(depth) - 1, max(depth.values())


def anytime(network, species=True):
    """[(cumulative evaluations, distinct root-connected chemical pairs)] after each attempt."""
    spent, curve = network.get('initialization_evaluations', 0), []
    for k, attempt in enumerate(network['attempts']):
        spent += attempt['evaluations']
        edges = [e for e in network['edges'] if e['attempt'] <= k]
        curve.append((spent, len(root_pairs(network, edges, species))))
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
    parser.add_argument('campaigns', type=Path, nargs='+')
    parser.add_argument('--budget', type=int, default=16000)
    parser.add_argument('--ids', type=Path, help='JSON list or {"ids": [...]} restricting the starts')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    keep = None
    if args.ids:
        data = json.loads(args.ids.read_text(encoding='utf-8'))
        keep = set(data if isinstance(data, list) else data['ids'])
    manifest = []
    for campaign in args.campaigns:
        protocol = Path(json.loads((campaign/'campaign.json').read_text(encoding='utf-8'))['protocol']).stem
        suffix = '' if protocol == 'FP-JCTC-1' else '@' + protocol
        tasks = csv.DictReader((campaign/'manifest.tsv').open(encoding='utf-8'), delimiter='\t')
        manifest += [dict(t, campaign=campaign, protocol=protocol, strategy=t['strategy'] + suffix)
                     for t in tasks if keep is None or t['start_id'] in keep]
    rows, missing = [], []
    for task in manifest:
        paths = list((task['campaign']/'runs'/task['outdir']).glob('*/*/network.json'))
        if not paths:
            missing.append(f"{task['campaign'].name}:{task['task_id']}")
            continue
        network = json.loads(paths[0].read_text(encoding='utf-8'))
        if network['status'] in ('running', 'aborted_error'):
            missing.append(f"{task['campaign'].name}:{task['task_id']}:{network['status']}")
            continue
        metrics = growth_metrics(network)
        curve, physical = anytime(network), anytime(network, species=False)
        new_species, depth = species_reach(network)
        rows.append(dict(campaign=task['campaign'].name, task=int(task['task_id']), start=task['start_id'],
                         strategy=task['strategy'], protocol=task['protocol'],
                         seed=int(task['seed']), status=network['status'], evaluations=network.get('evaluations'),
                         attempts=len(network['attempts']), edges=len(network['edges']),
                         root_chemical_pairs=len(root_pairs(network, network['edges'])),
                         pairs_at_budget=value_at(curve, args.budget),
                         new_species=new_species, max_species_depth=depth,
                         root_chemical_pairs_physical=len(root_pairs(network, network['edges'], species=False)),
                         pairs_at_budget_physical=value_at(physical, args.budget),
                         new_species_physical=len(metrics['root_species']) - 1,
                         longest_enumerated_chemical_path_physical=metrics['longest_enumerated_chemical_path'],
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
    for s in strategies:
        base = s.split('@')[0]
        if '@' in s and base in per_system:
            tests[f'{s}_vs_{base}'] = paired_test(per_system, s, base)
    statuses = {}
    for row in rows:
        statuses[row['status']] = statuses.get(row['status'], 0) + 1
    result = dict(campaigns=[str(c) for c in args.campaigns], ids=None if args.ids is None else str(args.ids),
                  budget=args.budget, tasks=len(manifest), finished=len(rows),
                  missing_tasks=missing, run_statuses=statuses, per_strategy=summary, paired_tests=tests,
                  per_system_means=per_system, rows=rows)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'per_system_means')}, indent=2))


if __name__ == '__main__':
    main()
