import pytest
from rdkit import Chem

from mechbridge.event_graph import graph_smiles
from mechbridge.oracle_library import OracleReferenceLibrary
from mechbridge.symbolic_library import parse_explicit


def diels_alder():
    rxn = '[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]>>[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:5][CH2:6]1'
    reactant = parse_explicit(rxn.split('>>')[0])
    numbers = [a.GetAtomicNum() for a in reactant.GetAtoms()]
    bonds = sorted([sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in reactant.GetBonds()])
    return reactant, dict(id='da', rxn_smiles=rxn, atomic_numbers=numbers, reactant_bonds=bonds)


def test_oracle_proposes_reference_edits_at_root_only():
    reactant, reference = diels_alder()
    library = OracleReferenceLibrary(reference, start=None)
    proposals = library.propose(reactant)
    assert len(proposals) == 1 and proposals[0]['knows_the_answer']
    changes = {tuple(e['atoms']): (e['before'], e['after']) for e in proposals[0]['edits']}
    index = {a.GetAtomMapNum(): a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomMapNum()}
    pair = lambda i, j: tuple(sorted((index[i], index[j])))
    assert changes == {pair(1, 2): (2, 1), pair(2, 3): (1, 2), pair(3, 4): (2, 1), pair(5, 6): (2, 1),
                       pair(4, 5): (0, 1), pair(1, 6): (0, 1)}
    assert proposals[0]['predicted_graph'] == 'C1=CCCCC1'
    assert proposals[0]['arrows'] == []
    other = parse_explicit('C=CC=C.CC')
    assert library.propose(other) == []
