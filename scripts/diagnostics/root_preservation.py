"""Does the relaxed root keep the reference reactant connectivity? (Fig. 4 sensitivity)

Every run relaxes the start geometry once (initialize_root: deterministic, no random
seed) and searches from that minimum. If the relaxation changes the connectivity, a
reference reaction can only be found after the network reaches the reference reactant
again. Per start with a finished run: the distinct root graphs over its finished runs
(they should agree), the start-file graph, and for every reference of the start (with
--groups every reaction of the group, as in score_reference_recovery.py) whether the
root bonds equal the reference reactant bonds under one reactant automorphism (mapped)
and whether the stereo-free root SMILES equals the reference reactant key (unmapped).
Reads the networks only; no potential evaluations.
Usage: root_preservation.py CAMPAIGN_DIR [...] --references FILE --out FILE [--groups GROUPS.json]
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.reference_matching import automorphisms, bond_set, mol_bonds, permute, stereo_free  # noqa: E402


def connectivity(start, positions):
    mol = geometry_mol(start['atomic_numbers'], np.asarray(positions), start['charge'])
    return stereo_free(graph_smiles(mol)), mol_bonds(mol)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaigns', type=Path, nargs='+')
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--groups', type=Path, help='group_t1x_reactants.py output')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    references = {r['id']: r for r in map(json.loads, args.references.read_text(encoding='utf-8').splitlines())}
    groups = json.loads(args.groups.read_text(encoding='utf-8'))['groups'] if args.groups else None
    roots, starts = {}, {}
    for campaign in args.campaigns:
        for path in sorted((campaign/'runs').rglob('network.json')):
            network = json.loads(path.read_text(encoding='utf-8'))
            if network['status'] in ('running', 'aborted_error') or not network['nodes']:
                continue
            start = network['start']
            starts[start['id']] = start
            roots.setdefault(start['id'], Counter())[connectivity(start, network['nodes'][0]['positions_A'])] += 1
    rows = []
    for name, counts in sorted(roots.items()):
        (root_key, root_bonds), _ = counts.most_common(1)[0]
        start = starts[name]
        start_key = connectivity(start, start['positions_A'])[0] if start.get('positions_A') else None
        checks = []
        for target in (groups.get(name, []) if groups is not None else [name]):
            reference = references.get(target)
            if reference is None or reference.get('admission', 'accepted') != 'accepted':
                continue
            react = bond_set(reference['reactant_bonds'])
            maps = automorphisms(reference['atomic_numbers'], reference['reactant_bonds'])
            checks.append(dict(reference=target, mapped=any(permute(root_bonds, m) == react for m in maps),
                               unmapped=root_key == reference.get('reactant_key')))
        rows.append(dict(start=name, runs=sum(counts.values()), root_graph=root_key,
                         distinct_root_graphs=len({k for k, _ in counts}), start_graph=start_key,
                         start_graph_equals_root=start_key == root_key, references=checks,
                         preserved=bool(checks) and all(c['mapped'] for c in checks)))
    scored = [r for r in rows if r['references']]
    result = dict(campaigns=[str(c) for c in args.campaigns], references=str(args.references),
                  starts=len(rows), starts_with_references=len(scored),
                  inconsistent_roots=[r['start'] for r in rows if r['distinct_root_graphs'] > 1],
                  preserved=sum(r['preserved'] for r in scored),
                  not_preserved=[dict(start=r['start'], root=r['root_graph'], start_file=r['start_graph'],
                                      references=len(r['references']))
                                 for r in scored if not r['preserved']],
                  reactions=sum(len(r['references']) for r in scored),
                  reactions_root_mapped=sum(c['mapped'] for r in scored for c in r['references']),
                  rows=rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=1))


if __name__ == '__main__':
    main()
