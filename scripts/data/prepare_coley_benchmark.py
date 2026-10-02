"""Coley/Stuyver [3+2] subset -> reactant-only encounter starts plus held-out references.

Eligible records: CHNO, two neutral closed-shell reactant fragments, at most
--max-atoms atoms with hydrogens. Strata: element triplet of the dipole (two
termini and the shared centre) and the dipolarophile bond; a seeded
proportional sample is drawn. Reactants are embedded from the mapped reactant
SMILES (ETKDG + UFF, as in v5) and assembled into one encounter per reaction;
author TS/product geometries are never used for starts. The mapped product
SMILES gives the reference product connectivity in the start's atom order;
hydrogens must stay on the same heavy atoms (checked, otherwise rejected).
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.encounters import assemble_encounter  # noqa: E402
from mechbridge.event_graph import geometry_mol, graph_smiles, resonance_equivalent  # noqa: E402
from mechbridge.reference_matching import stereo_free  # noqa: E402
from mechbridge.symbolic_library import parse_explicit  # noqa: E402


def mapped_bonds(mol):
    """Bonds between mapped atoms, as unordered map-number pairs."""
    return {tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum())))
            for b in mol.GetBonds() if b.GetBeginAtom().GetAtomMapNum() and b.GetEndAtom().GetAtomMapNum()}


def hydrogens_by_map(mol):
    return {a.GetAtomMapNum(): a.GetTotalNumHs(includeNeighbors=True)
            for a in mol.GetAtoms() if a.GetAtomMapNum()}


def stratum(reactant, product):
    """(dipole element triplet, dipolarophile bond) from the two newly formed bonds."""
    formed = mapped_bonds(product) - mapped_bonds(reactant)
    if len(formed) != 2:
        raise ValueError(f'expected two forming bonds, found {len(formed)}')
    by_map = {a.GetAtomMapNum(): a for a in reactant.GetAtoms() if a.GetAtomMapNum()}
    fragment = {a.GetAtomMapNum(): k for k, f in enumerate(Chem.GetMolFrags(reactant))
                for a in (reactant.GetAtomWithIdx(i) for i in f) if a.GetAtomMapNum()}
    ends = [sorted(bond, key=lambda m: fragment[m]) for bond in formed]
    dipole_side = {fragment[ends[0][0]], fragment[ends[0][1]]}
    if len(dipole_side) != 2:
        raise ValueError('forming bond within one fragment')
    termini_a, termini_b = [], []
    for a, b in ends:
        termini_a.append(a)
        termini_b.append(b)
    # The dipole is the fragment whose two termini share a common neighbour.
    for termini, partner in ((termini_a, termini_b), (termini_b, termini_a)):
        x, y = (by_map[m] for m in termini)
        shared = {n.GetIdx() for n in x.GetNeighbors()} & {n.GetIdx() for n in y.GetNeighbors()}
        if shared:
            centre = reactant.GetAtomWithIdx(min(shared))
            # A dipole read in either direction is the same type: C-N-N == N-N-C.
            triplet = min('-'.join([x.GetSymbol(), centre.GetSymbol(), y.GetSymbol()]),
                          '-'.join([y.GetSymbol(), centre.GetSymbol(), x.GetSymbol()]))
            p, q = (by_map[m] for m in partner)
            bond = reactant.GetBondBetweenAtoms(p.GetIdx(), q.GetIdx())
            dipolarophile = f"{p.GetSymbol()}{'-=#'[int(bond.GetBondTypeAsDouble())-1] if bond else '?'}{q.GetSymbol()}"
            return triplet, dipolarophile
    raise ValueError('no three-atom dipole found')


def allocate(strata, count, rng):
    """Largest-remainder proportional allocation, then a seeded draw per stratum."""
    total = sum(len(v) for v in strata.values())
    quotas = {k: count*len(v)/total for k, v in strata.items()}
    take = {k: int(q) for k, q in quotas.items()}
    for k in sorted(quotas, key=lambda k: (-(quotas[k]-take[k]), k))[:count-sum(take.values())]:
        take[k] += 1
    chosen = []
    for k in sorted(strata):
        members = sorted(strata[k], key=lambda r: int(r['rxn_id']))
        picks = rng.choice(len(members), size=min(take[k], len(members)), replace=False)
        chosen += [members[i] for i in sorted(picks)]
    return chosen, take


def embed(mol, seed):
    x = np.zeros((mol.GetNumAtoms(), 3))
    for ids, part in zip(Chem.GetMolFrags(mol), Chem.GetMolFrags(mol, asMols=True)):
        params = AllChem.ETKDGv3(); params.randomSeed = seed
        conformers = list(AllChem.EmbedMultipleConfs(part, numConfs=4, params=params))
        optimized = AllChem.UFFOptimizeMoleculeConfs(part, maxIters=500, numThreads=1)
        good = [i for i, (status, _) in enumerate(optimized) if status == 0]
        if not good:
            raise ValueError('no converged reactant conformer')
        best = min(good, key=lambda i: optimized[i][1])
        x[list(ids)] = part.GetConformer(conformers[best]).GetPositions()
    return x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, default=ROOT/'data/raw/coley_dipolar/full_dataset.csv')
    parser.add_argument('--count', type=int, default=100)
    parser.add_argument('--max-atoms', type=int, default=40)
    parser.add_argument('--seed', type=int, default=20261002)
    parser.add_argument('--starts', type=Path, default=ROOT/'data/processed/coley_benchmark_starts.jsonl')
    parser.add_argument('--references', type=Path,
                        default=ROOT/'data/processed/coley_benchmark_references.jsonl')
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    for path in (args.starts, args.references):
        if path.exists():
            raise FileExistsError(path)
    rows = list(csv.DictReader(args.csv.open(encoding='utf-8')))
    strata, audit, seen = {}, dict(records=len(rows), rejected={}), set()
    for row in rows:
        reactant_smiles, product_smiles = row['rxn_smiles'].split('>>')
        try:
            reactant = parse_explicit(reactant_smiles)
            product = parse_explicit(product_smiles)
            if not {a.GetAtomicNum() for a in reactant.GetAtoms()} <= {1, 6, 7, 8}:
                raise ValueError('elements outside CHNO')
            if len(Chem.GetMolFrags(reactant)) != 2 or Chem.GetFormalCharge(reactant) != 0:
                raise ValueError('not two fragments with zero net charge')
            if any(a.GetNumRadicalElectrons() for a in reactant.GetAtoms()):
                raise ValueError('radical reactant')
            if reactant.GetNumAtoms() > args.max_atoms:
                raise ValueError('too many atoms')
            maps = sorted(a.GetAtomMapNum() for a in reactant.GetAtoms() if a.GetAtomMapNum())
            if maps != sorted(a.GetAtomMapNum() for a in product.GetAtoms() if a.GetAtomMapNum()):
                raise ValueError('reactant and product atom maps differ')
            if hydrogens_by_map(reactant) != hydrogens_by_map(product):
                raise ValueError('hydrogen moves between heavy atoms')
            key = (stereo_free(graph_smiles(reactant)), stereo_free(graph_smiles(product)))
            if key in seen:
                raise ValueError('duplicate reaction')
            seen.add(key)
            row['_stratum'] = '|'.join(stratum(reactant, product))
        except ValueError as exc:
            reason = str(exc)
            audit['rejected'][reason] = audit['rejected'].get(reason, 0) + 1
            continue
        strata.setdefault(row['_stratum'], []).append(row)
    audit['eligible'] = sum(len(v) for v in strata.values())
    rng = np.random.default_rng(args.seed)
    chosen, take = allocate(strata, args.count, rng)
    starts, references = [], []
    for row in chosen:
        reactant_smiles, product_smiles = row['rxn_smiles'].split('>>')
        reactant, product = parse_explicit(reactant_smiles), parse_explicit(product_smiles)
        index_of = {a.GetAtomMapNum(): a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomMapNum()}
        heavy_product = {tuple(sorted((index_of[i], index_of[j]))) for i, j in mapped_bonds(product)}
        hydrogen = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in reactant.GetBonds()
                    if 1 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())}
        reactant_bonds = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
                                for b in reactant.GetBonds())
        product_bonds = sorted(heavy_product | hydrogen)
        plain = Chem.Mol(reactant)
        for atom in plain.GetAtoms():
            atom.SetAtomMapNum(0)
        ident = f"coley_{row['rxn_id']}"
        seed = 40000 + int(row['rxn_id'])
        record = dict(id=ident, dataset='Coley/Stuyver [3+2]', rxn_id=row['rxn_id'],
                      rxn_smiles=row['rxn_smiles'], stratum=row['_stratum'],
                      atomic_numbers=[a.GetAtomicNum() for a in plain.GetAtoms()],
                      reactant_bonds=[list(b) for b in reactant_bonds],
                      product_bonds=[list(b) for b in product_bonds],
                      reactant_key=stereo_free(graph_smiles(plain)),
                      product_key=stereo_free(graph_smiles(product)),
                      solution_G_act_kcal_mol=float(row['G_act']), solvent=row['solvent'],
                      reference_level='B3LYP-D3(BJ)/def2-TZVP//def2-SVP, SMD water (not comparable to gas-phase MLIP)')
        try:
            x = embed(plain, seed)
            numbers = record['atomic_numbers']
            xyz = assemble_encounter(numbers, x, Chem.GetMolFrags(plain), seed)
            observed = geometry_mol(numbers, xyz, 0)
            # As in v5: compare connectivity/resonance only; the 3D embedding may
            # assign stereo the SMILES leaves unspecified.
            expected, perceived = Chem.Mol(plain), Chem.Mol(observed)
            Chem.RemoveStereochemistry(expected); Chem.RemoveStereochemistry(perceived)
            if not resonance_equivalent(expected, perceived):
                raise ValueError('encounter geometry changes the perceived reactant graph')
            record['admission'] = 'accepted'
            starts.append(dict(id=ident, atomic_numbers=numbers, positions_A=xyz.tolist(), charge=0,
                multiplicity=1, provenance=dict(dataset='Coley/Stuyver [3+2]', rxn_id=row['rxn_id'],
                    input_smiles=graph_smiles(observed), random_seed=seed,
                    assembly='ETKDG+UFF fragments, independent SO(3) rotations, vdW gap 0.2 A',
                    reference_TS_or_product_geometry_used=False)))
        except ValueError as exc:
            record.update(admission='rejected', admission_reason=str(exc))
        references.append(record)
    for path, data in ((args.starts, starts), (args.references, references)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(''.join(json.dumps(r) + '\n' for r in data), encoding='utf-8')
    audit.update(strata={k: len(v) for k, v in sorted(strata.items())}, allocation=take,
                 sampled=len(chosen), accepted=len(starts), seed=args.seed,
                 starts_sha256=hashlib.sha256(args.starts.read_bytes()).hexdigest(),
                 references_sha256=hashlib.sha256(args.references.read_bytes()).hexdigest())
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
