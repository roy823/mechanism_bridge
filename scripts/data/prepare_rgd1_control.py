"""RGD1 seen-domain positive control: frame, stratified sample, starts + references.

AIMNet2-rxn was trained on RGD1, so this set only measures the gain from
training overlap and never enters generalization claims (B doc 4.1).
Frame (CSV level, every row): mapped reactant and product SMILES parse with
all hydrogens explicit, map numbers are 1..N, both sides are net neutral and
closed shell, CHNO only, and the connectivities differ. Rows are deduplicated
by the unordered pair of stereo-free graph keys, keeping the row with the
lowest TS energy in the HDF5 (TS conformers of one reaction share a pair).
Strata: net connectivity change b1f1, b2f2, b3f3+ (both counts >= 3), other.
Allocation: proportional to frame stratum sizes (largest remainder).
Sampling: one fixed permutation per stratum from default_rng(seed); candidates
are examined in that order and accepted until the quota is met when the HDF5
group agrees with the CSV: same element order as the atom maps (map k ->
atom k-1), reactant connectivity perceived from RG equals the CSV reactant,
likewise PG versus the CSV product when the group stores PG, and TSG present.
Starts hold only RG; references hold everything used for scoring.
Usage: prepare_rgd1_control.py --csv RGD1CHNO_AMsmiles.csv --h5 RGD1_allrxns.h5 [--size 100]
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.reference_matching import mol_bonds, stereo_free  # noqa: E402
from mechbridge.symbolic_library import parse_explicit  # noqa: E402

HARTREE_EV = 27.211386245988
STRATA = ('b1f1', 'b2f2', 'b3f3+', 'other')


def mapped_side(smiles):
    """(mol, atomic numbers by map, connectivity by map) of one side; maps k -> atom k-1."""
    mol = parse_explicit(smiles)
    maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    if sorted(maps) != list(range(1, len(maps) + 1)):
        raise ValueError('atom maps are not 1..N over all atoms')
    if sum(a.GetFormalCharge() for a in mol.GetAtoms()) or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError('charged or open-shell side')
    numbers = [0]*len(maps)
    for atom in mol.GetAtoms():
        numbers[atom.GetAtomMapNum() - 1] = atom.GetAtomicNum()
    bonds = frozenset(tuple(sorted((b.GetBeginAtom().GetAtomMapNum() - 1, b.GetEndAtom().GetAtomMapNum() - 1)))
                      for b in mol.GetBonds())
    return mol, numbers, bonds


def stratum(reactant, product):
    broken, formed = len(reactant - product), len(product - reactant)
    if broken == formed == 1:
        return 'b1f1'
    if broken == formed == 2:
        return 'b2f2'
    return 'b3f3+' if min(broken, formed) >= 3 else 'other'


def allocate(sizes, total):
    """Largest-remainder proportional allocation, capped by stratum size."""
    frame = sum(sizes.values())
    exact = {s: total*n/frame for s, n in sizes.items()}
    take = {s: min(int(v), sizes[s]) for s, v in exact.items()}
    for s in sorted(exact, key=lambda s: (-(exact[s] - int(exact[s])), s)):
        if sum(take.values()) >= total:
            break
        if take[s] < sizes[s]:
            take[s] += 1
    return take


def groups_by_id(handle):
    if all(name.startswith('MR_') for name in list(handle)[:10]):
        return {name: name for name in handle}
    return {handle[name]['_id'][()].decode(): name for name in handle}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, default=ROOT/'data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv')
    parser.add_argument('--h5', type=Path, default=ROOT/'data/raw/rgd1_zenodo/RGD1_allrxns.h5')
    parser.add_argument('--size', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20261002)
    parser.add_argument('--starts', type=Path, default=ROOT/'data/processed/rgd1_control_starts.jsonl')
    parser.add_argument('--references', type=Path, default=ROOT/'data/processed/rgd1_control_references.jsonl')
    parser.add_argument('--report', type=Path, default=ROOT/'data/processed/rgd1_control_sampling.json')
    args = parser.parse_args()
    for path in (args.starts, args.references, args.report):
        if path.exists():
            raise FileExistsError(path)
    RDLogger.DisableLog('rdApp.*')
    csv_rejected, frame = {}, {}
    with h5py.File(args.h5, 'r') as handle:
        index = groups_by_id(handle)
        with args.csv.open(encoding='utf-8-sig', newline='') as source:
            for row_number, row in enumerate(csv.DictReader(source)):
                try:
                    if row['reaction'] not in index:
                        raise ValueError('no HDF5 group')
                    mol_r, numbers, bonds_r = mapped_side(row['reactant'])
                    mol_p, numbers_p, bonds_p = mapped_side(row['product'])
                    if numbers != numbers_p or not set(numbers) <= {1, 6, 7, 8}:
                        raise ValueError('element mismatch between sides or outside CHNO')
                    if bonds_r == bonds_p:
                        raise ValueError('reactant and product share the same connectivity')
                    key = tuple(sorted((stereo_free(graph_smiles(mol_r)), stereo_free(graph_smiles(mol_p)))))
                    if 'TS_E' not in handle[index[row['reaction']]]:
                        raise ValueError('no TS energy in the HDF5 group')
                except (ValueError, RuntimeError) as exc:
                    reason = str(exc).split(':')[0]
                    csv_rejected[reason] = csv_rejected.get(reason, 0) + 1
                    continue
                ts_energy = float(handle[index[row['reaction']]]['TS_E'][()])
                best = frame.get(key)
                if best is None or ts_energy < best['ts_energy_Ha']:
                    frame[key] = dict(row=row_number, reaction=row['reaction'], reactant=row['reactant'],
                                      product=row['product'], ts_energy_Ha=ts_energy,
                                      stratum=stratum(bonds_r, bonds_p))
        members = {s: sorted((m for m in frame.values() if m['stratum'] == s), key=lambda m: m['row'])
                   for s in STRATA}
        sizes = {s: len(v) for s, v in members.items() if v}
        quota = allocate(sizes, args.size)
        rng = np.random.default_rng(args.seed)
        starts, references, walk = [], [], {}
        for s in STRATA:
            if s not in quota:
                continue
            order = rng.permutation(len(members[s]))
            walk[s] = dict(frame=sizes[s], quota=quota[s], examined=0, accepted=0, rejected={})
            for k in order:
                if walk[s]['accepted'] == quota[s]:
                    break
                candidate = members[s][k]
                walk[s]['examined'] += 1
                group = handle[index[candidate['reaction']]]
                try:
                    mol_r, numbers, bonds_r = mapped_side(candidate['reactant'])
                    mol_p, _, bonds_p = mapped_side(candidate['product'])
                    elements = np.asarray(group['elements']).astype(int).tolist()
                    if elements != numbers:
                        raise ValueError('HDF5 element order differs from the atom maps')
                    if 'TSG' not in group or 'RG' not in group:
                        raise ValueError('missing RG or TSG')
                    rg = np.asarray(group['RG'], dtype=float)
                    perceived_r = geometry_mol(elements, rg, 0)
                    if mol_bonds(perceived_r) != bonds_r:
                        raise ValueError('RG connectivity differs from the CSV reactant')
                    product_source, product_positions = 'csv_mapped_smiles', None
                    if 'PG' in group:
                        product_positions = np.asarray(group['PG'], dtype=float)
                        perceived_p = geometry_mol(elements, product_positions, 0)
                        if mol_bonds(perceived_p) != bonds_p:
                            raise ValueError('PG connectivity differs from the CSV product')
                        mol_p, product_source = perceived_p, 'hdf5_PG_geometry'
                except (ValueError, RuntimeError) as exc:
                    reason = str(exc).split(':')[0]
                    walk[s]['rejected'][reason] = walk[s]['rejected'].get(reason, 0) + 1
                    continue
                walk[s]['accepted'] += 1
                ident = f"rgd1_{candidate['reaction']}"
                graph_r, graph_p = graph_smiles(perceived_r), graph_smiles(mol_p)
                energies = {k: float(group[k][()]) for k in ('R_E', 'TS_E', 'P_E') if k in group}
                record = dict(id=ident, dataset='RGD1', reaction=candidate['reaction'], csv_row=candidate['row'],
                              stratum=s, atomic_numbers=elements, energies_Ha=energies,
                              forward_barrier_eV=(energies['TS_E'] - energies['R_E'])*HARTREE_EV
                              if {'TS_E', 'R_E'} <= set(energies) else None,
                              ts_positions_A=np.asarray(group['TSG'], dtype=float).tolist(),
                              admission='accepted', reactant_graph=graph_r, product_graph=graph_p,
                              reactant_key=stereo_free(graph_r), product_key=stereo_free(graph_p),
                              reactant_bonds=sorted(map(list, bonds_r)), product_bonds=sorted(map(list, bonds_p)),
                              product_graph_source=product_source, training_overlap='RGD1 is in the AIMNet2-rxn training set')
                if product_positions is not None:
                    record['product_positions_A'] = product_positions.tolist()
                references.append(record)
                starts.append(dict(id=ident, atomic_numbers=elements, positions_A=rg.tolist(), charge=0,
                                   multiplicity=1, provenance=dict(dataset='RGD1', reaction=candidate['reaction'],
                                                                   input_smiles=graph_r,
                                                                   reference_TS_or_product_geometry_used=False)))
    for path, rows in ((args.starts, starts), (args.references, references)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(''.join(json.dumps(r) + '\n' for r in rows), encoding='utf-8')
    report = dict(seed=args.seed, size=args.size, csv_sha256=hashlib.sha256(args.csv.read_bytes()).hexdigest(),
                  h5_name=args.h5.name, csv_rejected=csv_rejected, frame_unique_pairs=len(frame),
                  strata=walk, accepted=len(starts),
                  starts_sha256=hashlib.sha256(args.starts.read_bytes()).hexdigest(),
                  references_sha256=hashlib.sha256(args.references.read_bytes()).hexdigest(),
                  definitions=__doc__)
    args.report.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'definitions'}, indent=2))


if __name__ == '__main__':
    main()
