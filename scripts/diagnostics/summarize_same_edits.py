"""Summarize Fig. 5b runs (same_edits_arrows.py): TS clusters per arrow set and paired tests.

TS identity (G doc 4.3), for accepted attempts of one group:
  same      : identical endpoint graphs, |dE_TS| <= 0.02 eV and reaction-core RMSD
              <= 0.10 A (core = edited atoms and their first neighbours; Kabsch);
  different : different endpoints, |dE_TS| > 0.04 eV or core RMSD > 0.30 A;
  otherwise uncertain. Clusters are the connected components of 'same'.
  (The asynchrony criterion of G doc 4.3 is not applied.)
Intended success (G doc 4.3): an accepted attempt whose two endpoint graphs are
the source graph and the predicted product graph (stereo-free comparison).
Tests per group:
  arrow sets: permutation test (10000 draws, seed 20261002) of the chi-square
              statistic of the arrow-set x TS-cluster table of accepted attempts;
  each arrow set versus bond_edits: exact McNemar tests of accepted and of
              intended success, paired by (sample, random seed) over the seeds
              both arms ran.
Holm correction is applied within each test family across groups.
Usage: summarize_same_edits.py RUN.json [...] --out FILE
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.reaction_network import aligned_rmsd  # noqa: E402
from mechbridge.reference_matching import stereo_free  # noqa: E402

ACCEPTED = ('validated_descents', 'validated_core_descents')


def core_atoms(mol, edits):
    atoms = {i for e in edits for i in e['atoms']}
    return sorted(atoms | {n.GetIdx() for i in atoms for n in mol.GetAtomWithIdx(int(i)).GetNeighbors()})


def relation(a, b, core):
    if a['endpoint_graphs'] != b['endpoint_graphs']:
        return 'different'
    de = abs(a['ts_energy_eV'] - b['ts_energy_eV'])
    rmsd = aligned_rmsd(np.asarray(a['ts_positions_A'])[core], np.asarray(b['ts_positions_A'])[core])
    if de <= .02 and rmsd <= .10:
        return 'same'
    if de > .04 or rmsd > .30:
        return 'different'
    return 'uncertain'


def chi2(table):
    table = np.asarray(table, dtype=float)
    expected = table.sum(1, keepdims=True)*table.sum(0, keepdims=True)/table.sum()
    mask = expected > 0
    return float((((table - expected)**2)[mask]/expected[mask]).sum())


def permutation_p(arms, clusters, draws=10000, seed=20261002):
    arm_ids = sorted(set(arms)); cluster_ids = sorted(set(clusters))
    if len(arm_ids) < 2 or len(cluster_ids) < 2:
        return None
    def table(labels):
        return [[sum(1 for a, c in zip(labels, clusters) if a == x and c == y) for y in cluster_ids] for x in arm_ids]
    observed = chi2(table(arms))
    rng = np.random.default_rng(seed)
    labels = np.array(arms)
    hits = sum(chi2(table(list(rng.permutation(labels)))) >= observed - 1e-12 for _ in range(draws))
    return (hits + 1)/(draws + 1)


def intended(row, source, predicted):
    if row['status'] not in ACCEPTED or not row.get('endpoint_graphs'):
        return False
    try:
        return sorted(stereo_free(g) for g in row['endpoint_graphs']) == sorted([source, predicted])
    except ValueError:
        return False


def mcnemar(rows_a, rows_b, success=lambda r: r['status'] in ACCEPTED):
    """Exact McNemar on success, paired by (sample, seed)."""
    ok_a = {(r['sample'], r['seed']): success(r) for r in rows_a}
    ok_b = {(r['sample'], r['seed']): success(r) for r in rows_b}
    pairs = sorted(set(ok_a) & set(ok_b))
    a_only = sum(ok_a[k] and not ok_b[k] for k in pairs)
    b_only = sum(ok_b[k] and not ok_a[k] for k in pairs)
    p = binomtest(a_only, a_only + b_only, .5).pvalue if a_only + b_only else None
    return dict(pairs=len(pairs), arrows_only=a_only, bond_edits_only=b_only, p=p)


def holm(ps):
    order = sorted((p, k) for k, p in enumerate(ps) if p is not None)
    adjusted, running = [None]*len(ps), 0.
    for rank, (p, k) in enumerate(order):
        running = max(running, min(1., (len(order) - rank)*p))
        adjusted[k] = running
    return adjusted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', type=Path, nargs='+')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = []
    for path in args.runs:
        report = json.loads(path.read_text(encoding='utf-8'))
        mol = geometry_mol(report['atomic_numbers'], np.asarray(report['root_positions_A']), report['charge'])
        source = stereo_free(graph_smiles(mol))
        for group in report['groups']:
            predicted = stereo_free(group['predicted_graph'])
            hit = lambda r: intended(r, source, predicted)
            g = group['group']
            rows = [r for r in report['rows'] if r['group'] == g]
            core = core_atoms(mol, group['edits'])
            accepted = [r for r in rows if r['status'] in ACCEPTED and r.get('ts_positions_A')]
            graph = nx.Graph()
            graph.add_nodes_from(range(len(accepted)))
            counts = dict(same=0, different=0, uncertain=0)
            for i, j in combinations(range(len(accepted)), 2):
                kind = relation(accepted[i], accepted[j], core)
                counts[kind] += 1
                if kind == 'same':
                    graph.add_edge(i, j)
            cluster = {}
            for c, members in enumerate(sorted(nx.connected_components(graph), key=lambda m: -len(m))):
                for i in members:
                    cluster[i] = c
            arms = sorted({r['arm'] for r in rows})
            per_arm = {}
            for arm in arms:
                mine = [r for r in rows if r['arm'] == arm]
                idx = [i for i, r in enumerate(accepted) if r['arm'] == arm]
                per_arm[arm] = dict(attempts=len(mine), accepted=len(idx), intended=sum(hit(r) for r in mine),
                                    evaluations=sum(r['evaluations'] for r in mine),
                                    clusters=sorted({cluster[i] for i in idx}))
            arrow_idx = [i for i, r in enumerate(accepted) if r['arm'] != 'bond_edits']
            entry = dict(start=report['start'], group=g, predicted_graph=group['predicted_graph'],
                         arrow_sets=[s['template_id'] for s in group['arrow_sets']],
                         from_resonance=[('resonance' in (s.get('origin') or '')) for s in group['arrow_sets']],
                         core_atoms=core, pair_relations=counts, clusters=len(set(cluster.values())),
                         per_arm=per_arm,
                         arrow_set_cluster_p=permutation_p([accepted[i]['arm'] for i in arrow_idx],
                                                           [cluster[i] for i in arrow_idx]),
                         versus_bond_edits={arm: mcnemar([r for r in rows if r['arm'] == arm],
                                                         [r for r in rows if r['arm'] == 'bond_edits'])
                                            for arm in arms if arm != 'bond_edits'},
                         intended_versus_bond_edits={arm: mcnemar([r for r in rows if r['arm'] == arm],
                                                                  [r for r in rows if r['arm'] == 'bond_edits'], hit)
                                                     for arm in arms if arm != 'bond_edits'})
            out.append(entry)
    for family in ('arrow_set_cluster_p',):
        for e, p in zip(out, holm([e[family] for e in out])):
            e[family + '_holm'] = p
    for family in ('versus_bond_edits', 'intended_versus_bond_edits'):
        flat = [(e, arm) for e in out for arm in e[family]]
        for (e, arm), p in zip(flat, holm([e[family][arm]['p'] for e, arm in flat])):
            e[family][arm]['p_holm'] = p
    args.out.write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps([{k: v for k, v in e.items() if k not in ('core_atoms',)} for e in out], indent=1))


if __name__ == '__main__':
    main()
