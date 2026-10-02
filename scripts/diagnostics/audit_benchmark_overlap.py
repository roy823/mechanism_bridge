"""Training-overlap labels for the Fig. 4 benchmarks against RGD1 (B doc 2.4).

AIMNet2-rxn was trained on RGD1, so every benchmark reaction is labelled
before any method is scored. Keys are computed from stereo-free molecules:
  K1 reaction key: unordered pair of the two sides, each side the sorted tuple
     of standard InChIKeys of its components;
  K2 reactant key: InChIKey connectivity blocks (first 14 characters) of the
     reactant components; 'all' = every component occurs among the components
     of any RGD1 endpoint (reactant or product side), 'any' = at least one;
  K3 composition key: molecular formula of the whole reaction.
Label (first that applies): seen_reaction (K1), seen_reactant (K2 all),
seen_formula (K3), unseen. Standard InChI treats some heteroatom H as mobile,
so tautomers can share a key; this errs towards labelling a reaction as seen.
'strict_pair' repeats K1 with stereo-free canonical SMILES (no mobile-H
merging) for information. Not covered: the AIMNet2 base training structures,
the unpublished enumerated reactions of AIMNet2-rxn, and geometry keys (K4).
Usage: audit_benchmark_overlap.py --out FILE [--rgd1 CSV]
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.provenance import git_provenance  # noqa: E402

BENCHMARKS = ('data/processed/t1x_test_references.jsonl', 'data/processed/coley_benchmark_references.jsonl')
LABELS = ('seen_reaction', 'seen_reactant', 'seen_formula', 'unseen')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def side(smiles):
    """(InChIKeys, connectivity blocks, canonical SMILES, formula) of one stereo-free side."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'unparsable SMILES {smiles}')
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.RemoveStereochemistry(mol)
    mol = Chem.RemoveHs(mol)
    keys, smiles_parts = [], []
    for fragment in Chem.GetMolFrags(mol, asMols=True):
        key = Chem.MolToInchiKey(fragment)
        if not key:
            raise ValueError(f'no InChIKey for {Chem.MolToSmiles(fragment)}')
        keys.append(key)
        smiles_parts.append(Chem.MolToSmiles(fragment, isomericSmiles=False))
    return (tuple(sorted(keys)), frozenset(k[:14] for k in keys), tuple(sorted(smiles_parts)),
            rdMolDescriptors.CalcMolFormula(mol))


def pair(a, b):
    return tuple(sorted((a, b)))


def rgd1_index(path):
    index = dict(k1=defaultdict(list), strict=defaultdict(list), blocks=set(), formulas=set())
    stats = dict(rows=0, parsed=0, parse_errors=0)
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            stats['rows'] += 1
            try:
                reactant, product = side(row['reactant']), side(row['product'])
            except ValueError:
                stats['parse_errors'] += 1
                continue
            stats['parsed'] += 1
            index['k1'][pair(reactant[0], product[0])].append(row['reaction'])
            index['strict'][pair(reactant[2], product[2])].append(row['reaction'])
            index['blocks'] |= reactant[1] | product[1]
            index['formulas'].add(reactant[3])
    return index, stats


def label(reference, index):
    """Overlap label from the reference's stereo-free endpoint keys (the benchmark identity)."""
    reactant, product = side(reference['reactant_key']), side(reference['product_key'])
    k1 = index['k1'].get(pair(reactant[0], product[0]), [])
    strict = index['strict'].get(pair(reactant[2], product[2]), [])
    k2_all = reactant[1] <= index['blocks']
    k2_any = bool(reactant[1] & index['blocks'])
    k3 = reactant[3] in index['formulas']
    name = ('seen_reaction' if k1 else 'seen_reactant' if k2_all else 'seen_formula' if k3 else 'unseen')
    return dict(label=name, k1_rgd1_ids=k1[:5], k1_hits=len(k1), strict_pair_hits=len(strict),
                k2_all=k2_all, k2_any=k2_any, k3=k3, formula=reactant[3],
                identity_reaction_under_inchi=reactant[0] == product[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rgd1', type=Path, default=ROOT/'data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    RDLogger.DisableLog('rdApp.*')
    index, stats = rgd1_index(args.rgd1)
    labels, counts, failures = {}, {}, {}
    for name in BENCHMARKS:
        dataset = Path(name).stem
        counts[dataset] = dict.fromkeys(LABELS, 0)
        for line in (ROOT/name).read_text(encoding='utf-8').splitlines():
            reference = json.loads(line)
            try:
                labels[reference['id']] = result = label(reference, index)
            except ValueError as exc:
                failures[reference['id']] = str(exc)
                continue
            counts[dataset][result['label']] += 1
    report = dict(inputs=dict(rgd1=dict(path=str(args.rgd1), sha256=sha256(args.rgd1)),
                              **{Path(n).stem: dict(path=n, sha256=sha256(ROOT/n)) for n in BENCHMARKS}),
                  code=git_provenance(ROOT), rdkit=rdBase.rdkitVersion, rgd1_stats=stats,
                  rgd1_distinct_k1_pairs=len(index['k1']), counts=counts, failures=failures,
                  definitions=__doc__, labels=labels)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k not in ('labels', 'definitions')}, indent=2))


if __name__ == '__main__':
    main()
