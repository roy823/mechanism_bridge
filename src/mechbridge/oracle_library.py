"""Given-reaction-centre ("oracle") proposals for the engine-only benchmark track.

This library knows the answer: it proposes exactly the reference reaction's
net bond edits, and only at the root reactant. It isolates the 3D engine
(how efficiently a correct bond-edit hypothesis becomes a validated TS) from
proposal coverage, and every run using it is labelled as an oracle run.
Arrows are not available for these references, so only edit-based strategies
(bond_edits, center_random) and the pyGSM edits arm can use it.
"""
import numpy as np
from rdkit import Chem

from .event_graph import bond_orders, geometry_mol, graph_smiles
from .symbolic_library import parse_explicit


def _kekule_orders(mol):
    """Integer bond orders (aromatic bonds kekulized) keyed by sorted atom pairs."""
    copy = Chem.Mol(mol)
    Chem.Kekulize(copy, clearAromaticFlags=True)
    return {k: int(round(v)) for k, v in bond_orders(copy).items()}


def reference_edits(reference, start):
    """(edits, predicted product graph) in the start's atom order."""
    numbers = reference['atomic_numbers']
    if 'product_positions_A' in reference:          # T1x: both sides perceived from geometry
        reactant = geometry_mol(numbers, np.asarray(start['positions_A']), 0)
        product = geometry_mol(numbers, np.asarray(reference['product_positions_A']), 0)
        before, after = _kekule_orders(reactant), _kekule_orders(product)
        predicted = graph_smiles(product)
    else:                                           # Coley: mapped SMILES, hydrogens fixed
        reactant_smiles, product_smiles = reference['rxn_smiles'].split('>>')
        reactant, product = parse_explicit(reactant_smiles), parse_explicit(product_smiles)
        index = {a.GetAtomMapNum(): a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomMapNum()}
        before = _kekule_orders(reactant)
        after = {k: v for k, v in before.items() if 1 in (numbers[k[0]], numbers[k[1]])}
        kekule = Chem.Mol(product)
        Chem.Kekulize(kekule, clearAromaticFlags=True)
        for bond in kekule.GetBonds():
            i, j = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if i and j:
                after[tuple(sorted((index[i], index[j])))] = int(round(bond.GetBondTypeAsDouble()))
        predicted = graph_smiles(product)
    edits = [dict(atoms=list(pair), before=int(before.get(pair, 0)), after=int(after.get(pair, 0)))
             for pair in sorted(set(before) | set(after)) if before.get(pair, 0) != after.get(pair, 0)]
    if not edits:
        raise ValueError('reference reaction has no bond edits')
    return edits, predicted


class OracleReferenceLibrary:
    """Proposes the reference reaction once, at the root reactant only."""
    policy = 'oracle_reference_net_edits_v1 (knows the answer; engine-only track)'

    def __init__(self, reference, start):
        self.reference = reference
        self.edits, self.predicted = reference_edits(reference, start)
        self.root_bonds = frozenset(tuple(sorted(b)) for b in reference['reactant_bonds'])
        self.audit = dict(reference=reference['id'], edits=len(self.edits))

    def propose(self, mol, limit=24):
        bonds = frozenset(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds())
        if bonds != self.root_bonds:
            return []
        return [dict(template_id='oracle:'+self.reference['id'], name='reference net edits',
                     edits=[dict(e) for e in self.edits], arrows=[], predicted_graph=self.predicted,
                     origin='oracle_reference_net_edits', intermolecular=len(Chem.GetMolFrags(mol)) > 1,
                     knows_the_answer=True)]
