"""Score blind searches against held-out reference reactions (Fig. 4 criteria).

S1 (mapped): an accepted MLIP edge whose two endpoint connectivities equal the
reference reactant and product connectivities under one automorphism of the
reactant graph (atom order is shared because starts keep the reference order).
S1 (unmapped): the unordered pair of stereo-free endpoint SMILES equals the
reference pair. Evaluations-to-hit count initialization plus every attempt up
to and including the first hit, in completion order (serial runs).
With --groups (group_t1x_reactants.py), a run from a group representative is
scored against every reaction of that group. Rows that share a start are not
independent (seeds; reactions of one reactant), so the summary adds a 95%
percentile interval from a bootstrap over starts (clusters); the Wilson
interval treats rows as independent and is kept for reference only.

--representable restricts the main denominator to representable reactions;
--labels adds layers by training-overlap label and the pre-registered main
layer (unseen plus seen_formula) on that denominator.
--root-check (root_preservation.py output) adds, for every subset, a copy
without the starts whose MLIP-relaxed root does not keep the reference
reactant connectivity (sensitivity analysis, plan §5, 2026-10-02 22:05 EDT).
--from-scores reuses the rows of an earlier output instead of scoring the
networks again (the subsets and intervals are recomputed).

Usage: python score_reference_recovery.py --runs RUN_DIR [RUN_DIR ...] --references FILE --out FILE
         [--groups GROUPS.json] [--representable IDS.json] [--labels OVERLAP.json] [--root-check FILE]
       python score_reference_recovery.py --from-scores OLD.json --references FILE --out FILE [...]
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.reference_matching import (MAX_AUTOMORPHISMS, automorphisms, bond_set,  # noqa: E402
                                           edge_matches, mol_bonds, stereo_free, wilson)

def score_network(network, reference, maps):
    numbers = network['start']['atomic_numbers']
    charge = network['start']['charge']
    react, prod = bond_set(reference['reactant_bonds']), bond_set(reference['product_bonds'])
    keys = sorted([reference['reactant_key'], reference['product_key']])
    cache = {}

    def node_info(index):
        if index not in cache:
            mol = geometry_mol(numbers, np.asarray(network['nodes'][index]['positions_A']), charge)
            cache[index] = (mol_bonds(mol), stereo_free(graph_smiles(mol)))
        return cache[index]

    attempts = network['attempts']
    cumulative, spent = [], network.get('initialization_evaluations', 0)
    for attempt in attempts:
        spent += attempt['evaluations']
        cumulative.append(spent)
    hits = dict(mapped=None, unmapped=None)
    for edge in network['edges']:
        (bonds_u, key_u), (bonds_v, key_v) = node_info(edge['nodes'][0]), node_info(edge['nodes'][1])
        mapped = edge_matches(bonds_u, bonds_v, react, prod, maps)
        unmapped = sorted([key_u, key_v]) == keys
        for name, found in (('mapped', mapped), ('unmapped', unmapped)):
            if found and (hits[name] is None or edge['attempt'] < hits[name]['attempt']):
                hits[name] = dict(edge=edge['id'], attempt=edge['attempt'],
                                  evaluations_to_hit=cumulative[edge['attempt']])
    statuses = {}
    for attempt in attempts:
        statuses[attempt['status']] = statuses.get(attempt['status'], 0) + 1
    return dict(status=network['status'], attempts=len(attempts), evaluations=network['evaluations'],
                edges=len(network['edges']), attempt_statuses=statuses,
                S1_mapped=hits['mapped'] is not None, S1_unmapped=hits['unmapped'] is not None,
                first_mapped_hit=hits['mapped'], first_unmapped_hit=hits['unmapped'],
                automorphisms=len(maps), automorphisms_truncated=len(maps) >= MAX_AUTOMORPHISMS)


def cluster_bootstrap(rows, key, draws, seed):
    """95% percentile interval of the hit rate, resampling start clusters with replacement."""
    clusters = {}
    for row in rows:
        counts = clusters.setdefault(row['start'], [0, 0])
        counts[0] += row[key]
        counts[1] += 1
    if not clusters:
        return None
    hits, n = np.array(list(clusters.values()), dtype=float).T
    picks = np.random.default_rng(seed).integers(0, len(hits), size=(draws, len(hits)))
    rates = hits[picks].sum(1)/n[picks].sum(1)
    return [float(np.percentile(rates, 2.5)), float(np.percentile(rates, 97.5))]


MAIN_LAYERS = ('unseen', 'seen_formula')


def summarize(rows, draws, seed):
    """Per strategy: S1 counts, rates, Wilson and cluster-bootstrap intervals."""
    summary = {}
    for strategy in sorted({r['strategy'] for r in rows}):
        mine = [r for r in rows if r['strategy'] == strategy]
        entry = dict(runs=len(mine), reactions=len({r['reference'] for r in mine}),
                     starts=len({r['start'] for r in mine}))
        for key in ('S1_mapped', 'S1_unmapped'):
            hits = sum(bool(r[key]) for r in mine)
            entry[key] = hits
            entry[f'{key}_rate'] = hits/len(mine)
            entry[f'{key}_rate_wilson95'] = wilson(hits, len(mine))
            entry[f'{key}_rate_cluster_bootstrap95'] = cluster_bootstrap(mine, key, draws, seed)
        summary[strategy] = entry
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, nargs='+')
    parser.add_argument('--from-scores', type=Path, help='reuse the rows of an earlier output')
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--groups', type=Path, help='Score each representative run against its group')
    parser.add_argument('--representable', type=Path, help='JSON with representable_ids (main denominator)')
    parser.add_argument('--labels', type=Path, help='audit_benchmark_overlap.py output')
    parser.add_argument('--root-check', type=Path, help='root_preservation.py output (sensitivity subsets)')
    parser.add_argument('--bootstrap', type=int, default=10000)
    parser.add_argument('--bootstrap-seed', type=int, default=20261002)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    if not args.runs and not args.from_scores:
        parser.error('give --runs or --from-scores')
    references = {r['id']: r for r in map(json.loads,
                  args.references.read_text(encoding='utf-8').splitlines())}
    groups = json.loads(args.groups.read_text(encoding='utf-8'))['groups'] if args.groups else None
    rows = json.loads(args.from_scores.read_text(encoding='utf-8'))['rows'] if args.from_scores else []
    for run in args.runs or []:
        for path in sorted(run.rglob('network.json')):
            network = json.loads(path.read_text(encoding='utf-8'))
            start = network['start']['id']
            for target in (groups.get(start, []) if groups is not None else [start]):
                reference = references.get(target)
                if reference is None or reference.get('admission') != 'accepted':
                    continue
                maps = automorphisms(reference['atomic_numbers'], reference['reactant_bonds'])
                row = dict(start=start, reference=target, strategy=network['strategy'],
                           seed=network['protocol']['random_seed'], file=str(path),
                           **score_network(network, reference, maps))
                rows.append(row)
    labels = json.loads(args.labels.read_text(encoding='utf-8'))['labels'] if args.labels else {}
    keep = (set(json.loads(args.representable.read_text(encoding='utf-8'))['representable_ids'])
            if args.representable else None)
    for row in rows:
        row['representable'] = keep is None or row['reference'] in keep
        row['overlap_label'] = labels.get(row['reference'], {}).get('label')
    base = [r for r in rows if r['representable']]
    subsets = dict(representable=base) if keep is not None else {}
    if labels:
        for name in sorted({r['overlap_label'] for r in base} - {None}):
            subsets[f'label={name}'] = [r for r in base if r['overlap_label'] == name]
        subsets['main_layers'] = [r for r in base if r['overlap_label'] in MAIN_LAYERS]
    if args.root_check:
        changed = {r['start'] for r in json.loads(args.root_check.read_text(encoding='utf-8'))['not_preserved']}
        for name, subset in list(subsets.items()) + [('all', rows)]:
            subsets[f'{name}_root_preserved'] = [r for r in subset if r['start'] not in changed]
    summary = summarize(rows, args.bootstrap, args.bootstrap_seed)
    layers = {k: summarize(v, args.bootstrap, args.bootstrap_seed) for k, v in subsets.items() if v}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(dict(references=str(args.references), summary=summary, subsets=layers,
                                        main_layers=list(MAIN_LAYERS), rows=rows), indent=2), encoding='utf-8')
    print(json.dumps(dict(all=summary, **layers), indent=2))


if __name__ == '__main__':
    main()
