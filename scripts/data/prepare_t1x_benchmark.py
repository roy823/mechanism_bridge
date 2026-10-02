"""Transition1x official test split -> reactant-only starts plus held-out references.

Starts contain only the DFT reactant geometry (re-relaxed on the MLIP by the
search). Product graphs, bond lists and the reference TS go to a separate
references file that is read only for post-search scoring. Every test reaction
is kept and its admission status recorded; nothing is dropped silently.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402

ENERGY = 'wB97x_6-31G(d).energy'


def stereo_free(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Cannot parse graph SMILES: {smiles}')
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def connectivity(mol):
    return sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()] if b.GetBeginAtomIdx() < b.GetEndAtomIdx()
                  else [b.GetEndAtomIdx(), b.GetBeginAtomIdx()] for b in mol.GetBonds())


def frame(group, numbers):
    positions = np.asarray(group['positions'], dtype=float)
    positions = positions.reshape(-1, len(numbers), 3)[0]
    energy = float(np.asarray(group[ENERGY]).reshape(-1)[0]) if ENERGY in group else None
    return positions, energy


def reactions(handle, split):
    for formula in sorted(handle[split]):
        for rxn in sorted(handle[split][formula]):
            yield formula, rxn, handle[split][formula][rxn]


def describe(handle, split):
    formula, rxn, group = next(reactions(handle, split))
    layout = {name: (dict(shape=list(item.shape), dtype=str(item.dtype)) if isinstance(item, h5py.Dataset)
                     else sorted(item)) for name, item in group.items()}
    return dict(split=split, example=f'{formula}/{rxn}', layout=layout,
                reactions=sum(len(handle[split][f]) for f in handle[split]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--h5', type=Path, default=ROOT/'data/raw/transition1x/Transition1x.h5')
    parser.add_argument('--split', default='test')
    parser.add_argument('--starts', type=Path, default=ROOT/'data/processed/t1x_test_starts.jsonl')
    parser.add_argument('--references', type=Path,
                        default=ROOT/'data/processed/t1x_test_references.jsonl')
    parser.add_argument('--inspect', action='store_true', help='Print the HDF5 layout and stop')
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    with h5py.File(args.h5, 'r') as handle:
        if args.inspect:
            print(json.dumps(describe(handle, args.split), indent=2))
            return
        for path in (args.starts, args.references):
            if path.exists():
                raise FileExistsError(path)
        starts, references = [], []
        for formula, rxn, group in reactions(handle, args.split):
            source = group if 'atomic_numbers' in group else group['reactant']
            numbers = np.asarray(source['atomic_numbers']).astype(int).tolist()
            ident = f't1x_{args.split}_{formula}_{rxn}'
            reactant, e_r = frame(group['reactant'], numbers)
            product, e_p = frame(group['product'], numbers)
            ts, e_ts = frame(group['transition_state'], numbers)
            record = dict(id=ident, dataset='Transition1x', split=args.split, formula=formula,
                          reaction=rxn, atomic_numbers=numbers, energies_eV=dict(
                              reactant=e_r, transition_state=e_ts, product=e_p),
                          forward_barrier_eV=None if None in (e_r, e_ts) else e_ts - e_r,
                          ts_positions_A=ts.tolist(), product_positions_A=product.tolist())
            try:
                if not set(numbers) <= {1, 6, 7, 8}:
                    raise ValueError('elements outside CHNO')
                mols = []
                for side, positions in (('reactant', reactant), ('product', product)):
                    try:
                        mols.append(geometry_mol(numbers, positions, 0))
                    except ValueError as exc:
                        raise ValueError(f'{side} perception: {exc}') from exc
                mol_r, mol_p = mols
                graph_r, graph_p = graph_smiles(mol_r), graph_smiles(mol_p)
                if connectivity(mol_r) == connectivity(mol_p):
                    raise ValueError('reactant and product share the same connectivity')
                record.update(admission='accepted', reactant_graph=graph_r, product_graph=graph_p,
                              reactant_key=stereo_free(graph_r), product_key=stereo_free(graph_p),
                              reactant_bonds=connectivity(mol_r), product_bonds=connectivity(mol_p))
                starts.append(dict(id=ident, atomic_numbers=numbers, positions_A=reactant.tolist(),
                    charge=0, multiplicity=1, provenance=dict(dataset='Transition1x',
                        split=args.split, formula=formula, reaction=rxn, input_smiles=graph_r,
                        reference_TS_or_product_geometry_used=False)))
            except ValueError as exc:
                record.update(admission='rejected', admission_reason=str(exc))
            references.append(record)
    for path, rows in ((args.starts, starts), (args.references, references)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(''.join(json.dumps(r) + '\n' for r in rows), encoding='utf-8')
    summary = dict(split=args.split, reactions=len(references), accepted=len(starts),
                   rejected={}, starts_sha256=hashlib.sha256(args.starts.read_bytes()).hexdigest(),
                   references_sha256=hashlib.sha256(args.references.read_bytes()).hexdigest())
    for record in references:
        if record['admission'] == 'rejected':
            # Group by side and error class, e.g. "product perception: Valence ...".
            reason = record['admission_reason'].split(';')[0]
            reason = reason.split(' of atom ')[0] if 'Valence' in reason else reason
            summary['rejected'][reason] = summary['rejected'].get(reason, 0) + 1
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
