"""Score blind searches against held-out reference reactions (Fig. 4 criteria).

S1 (mapped): an accepted MLIP edge whose two endpoint connectivities equal the
reference reactant and product connectivities under one automorphism of the
reactant graph (atom order is shared because starts keep the reference order).
S1 (unmapped): the unordered pair of stereo-free endpoint SMILES equals the
reference pair. Evaluations-to-hit count initialization plus every attempt up
to and including the first hit, in completion order (serial runs).

Usage: python score_reference_recovery.py --runs RUN_DIR [RUN_DIR ...] --references FILE --out FILE
"""
import argparse
import json
import math
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402

MAX_AUTOMORPHISMS = 5000


def stereo_free(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def bond_set(bonds):
    return frozenset(tuple(sorted(map(int, b[:2]))) for b in bonds)


def automorphisms(numbers, bonds):
    graph = nx.Graph()
    graph.add_nodes_from((i, dict(z=int(z))) for i, z in enumerate(numbers))
    graph.add_edges_from(tuple(b[:2]) for b in bonds)
    matcher = GraphMatcher(graph, graph, node_match=lambda a, b: a['z'] == b['z'])
    maps = []
    for mapping in matcher.isomorphisms_iter():
        maps.append(mapping)
        if len(maps) >= MAX_AUTOMORPHISMS:
            break
    return maps


def permute(bonds, mapping):
    return frozenset(tuple(sorted((mapping[i], mapping[j]))) for i, j in bonds)


def edge_matches(bonds_u, bonds_v, react, prod, maps):
    """One reactant automorphism must map the endpoints onto (reactant, product)."""
    return any((permute(bonds_u, m) == react and permute(bonds_v, m) == prod) or
               (permute(bonds_v, m) == react and permute(bonds_u, m) == prod) for m in maps)


def wilson(successes, n, z=1.96):
    if not n:
        return None
    p = successes/n
    center = (p + z*z/(2*n))/(1 + z*z/n)
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1 + z*z/n)
    return [center - half, center + half]


def score_network(network, reference, maps):
    numbers = network['start']['atomic_numbers']
    charge = network['start']['charge']
    react, prod = bond_set(reference['reactant_bonds']), bond_set(reference['product_bonds'])
    keys = sorted([reference['reactant_key'], reference['product_key']])
    cache = {}

    def node_info(index):
        if index not in cache:
            mol = geometry_mol(numbers, np.asarray(network['nodes'][index]['positions_A']), charge)
            bonds = frozenset(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
                              for b in mol.GetBonds())
            cache[index] = (bonds, stereo_free(graph_smiles(mol)))
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, nargs='+', required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    references = {r['id']: r for r in map(json.loads,
                  args.references.read_text(encoding='utf-8').splitlines())}
    rows = []
    for run in args.runs:
        for path in sorted(run.rglob('network.json')):
            network = json.loads(path.read_text(encoding='utf-8'))
            reference = references.get(network['start']['id'])
            if reference is None or reference.get('admission') != 'accepted':
                continue
            maps = automorphisms(reference['atomic_numbers'], reference['reactant_bonds'])
            row = dict(start=network['start']['id'], strategy=network['strategy'],
                       seed=network['protocol']['random_seed'], file=str(path),
                       **score_network(network, reference, maps))
            rows.append(row)
    summary = {}
    for row in rows:
        entry = summary.setdefault(row['strategy'], dict(runs=0, S1_mapped=0, S1_unmapped=0))
        entry['runs'] += 1
        entry['S1_mapped'] += row['S1_mapped']
        entry['S1_unmapped'] += row['S1_unmapped']
    for entry in summary.values():
        entry['S1_mapped_rate_wilson95'] = wilson(entry['S1_mapped'], entry['runs'])
        entry['S1_unmapped_rate_wilson95'] = wilson(entry['S1_unmapped'], entry['runs'])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(dict(references=str(args.references), summary=summary, rows=rows),
                                   indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
