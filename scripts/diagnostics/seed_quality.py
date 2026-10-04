"""Seed quality against reference TSs: does the arrow information move the seed toward the TS?

For one start whose reference product the frozen library proposes (connectivity match under a
reactant automorphism, as in symbolic_coverage.py), seeds are built with make_seed from the
MLIP-relaxed root of a finished campaign run, exactly as the search builds them (protocol JSON
settings), for bond_edits, for arrows with the historical seed (legacy) and for arrows under every
seed_feature_ablation of arrow_features_v1, over samples 0-2 and random seeds 0..N-1 (the same
random seeds for every arm, so arms differ only in the seed features).
Metrics against the reference TS (author TS in the start atom order):
  core_rmsd_A : Kabsch RMSD over the edited atoms and their first neighbours;
  edit_error_A: mean |d_seed - d_TS| over the edited atom pairs.
No PES evaluation. One start per call, so starts can run as an array.
Usage: seed_quality.py --starts F --references F --campaign DIR --protocol-json P --start-id ID
         --out FILE [--seeds 3]
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol  # noqa: E402
from mechbridge.reaction_network import aligned_rmsd  # noqa: E402
from mechbridge.reference_matching import automorphisms, bond_set, mol_bonds, permute  # noqa: E402
from mechbridge.search_seeds import SEED_ABLATIONS, make_seed  # noqa: E402
from mechbridge.symbolic_library import ArrowLibrary, ResonanceAwareLibrary  # noqa: E402

ARMS = [('bond_edits', 'bond_edits', 'arrow_features_v1', 'none'),
        ('arrows_legacy', 'arrows', 'legacy', 'none')] + [
       (f'arrows:{a}', 'arrows', 'arrow_features_v1', a) for a in SEED_ABLATIONS]


def root_positions(campaign, start_id):
    for path in sorted((campaign/'runs').glob(f'{start_id}/*/*/{start_id}/*/network.json')):
        network = json.loads(path.read_text(encoding='utf-8'))
        if network['status'] not in ('running', 'aborted_error') and network['nodes']:
            return np.asarray(network['nodes'][0]['positions_A']), str(path)
    return None, None


def matching_proposals(proposals, mol, reference):
    """Proposals whose edits turn the root bonds into the reference product bonds."""
    node_bonds, target = mol_bonds(mol), bond_set(reference['product_bonds'])
    maps = automorphisms(reference['atomic_numbers'], reference['reactant_bonds'])
    seen, out = set(), []
    for p in proposals:
        product = set(node_bonds)
        for e in p['edits']:
            pair = tuple(sorted(int(a) for a in e['atoms']))
            if e['after'] == 0:
                product.discard(pair)
            elif e['before'] == 0:
                product.add(pair)
        key = json.dumps([p['edits'], p.get('arrows')], sort_keys=True)
        if key not in seen and any(permute(frozenset(product), m) == target for m in maps):
            seen.add(key)
            out.append(p)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starts', type=Path, required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--campaign', type=Path, required=True, help='campaign whose runs give the MLIP root')
    parser.add_argument('--protocol-json', type=Path, required=True)
    parser.add_argument('--start-id', required=True)
    parser.add_argument('--seeds', type=int, default=3)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    protocol = json.loads(args.protocol_json.read_text(encoding='utf-8'))
    start = next(s for s in map(json.loads, args.starts.read_text(encoding='utf-8').splitlines())
                 if s['id'] == args.start_id)
    reference = next(r for r in map(json.loads, args.references.read_text(encoding='utf-8').splitlines())
                     if r['id'] == args.start_id)
    result = dict(start=args.start_id, protocol=str(args.protocol_json), rows=[])
    root, source = root_positions(args.campaign, args.start_id)
    ts = reference.get('ts_positions_A')
    if root is None or ts is None:
        result['skipped'] = 'no finished root' if root is None else 'no reference TS'
        args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
        return
    ts = np.asarray(ts)
    numbers = start['atomic_numbers']
    mol = geometry_mol(numbers, root, start['charge'])
    library = ResonanceAwareLibrary(ArrowLibrary(ROOT/'data/raw/synepd/polar.json'), protocol['proposal_resonance_forms'])
    proposals = matching_proposals(library.propose(mol), mol, reference)
    result.update(root_source=source, matching_proposals=len(proposals))
    for k, proposal in enumerate(proposals):
        edits = proposal['edits']
        edited = sorted({int(a) for e in edits for a in e['atoms']})
        core = sorted(set(edited) | {n.GetIdx() for i in edited for n in mol.GetAtomWithIdx(i).GetNeighbors()})
        pairs = [tuple(int(a) for a in e['atoms']) for e in edits]
        ts_d = np.array([np.linalg.norm(ts[i]-ts[j]) for i, j in pairs])
        for arm, strategy, features, ablation in ARMS:
            for sample in range(3):
                for seed in range(args.seeds):
                    row = dict(proposal=k, template_id=proposal['template_id'], origin=proposal.get('origin'),
                               arm=arm, sample=sample, seed=seed)
                    try:
                        x, _, meta = make_seed(Atoms(numbers=numbers, positions=root), mol, strategy, proposal,
                                               sample, seed, protocol['symbolic_seed_scale'],
                                               protocol['encounter_policy'], features,
                                               protocol['seed_fit_max_nfev'], ablation)
                    except (ValueError, RuntimeError) as exc:
                        row['error'] = str(exc)[:200]
                        result['rows'].append(row)
                        continue
                    d = np.array([np.linalg.norm(x[i]-x[j]) for i, j in pairs])
                    row.update(core_rmsd_A=aligned_rmsd(x[core], ts[core]),
                               edit_error_A=float(np.mean(np.abs(d - ts_d))),
                               lone_pair_terms=(meta.get('arrow_terms') or {}).get('lone_pair_terms'))
                    result['rows'].append(row)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    ok = [r for r in result['rows'] if 'core_rmsd_A' in r]
    print(json.dumps(dict(start=args.start_id, proposals=len(proposals), rows=len(result['rows']), ok=len(ok))))


if __name__ == '__main__':
    main()
