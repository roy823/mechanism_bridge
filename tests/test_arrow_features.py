import numpy as np
import pytest
from ase import Atoms
from rdkit.Chem import AllChem

from mechbridge.arrow_features import arrow_chain, chain_progress_shift, lone_pair_terms
from mechbridge.event_graph import graph_smiles
from mechbridge.search_seeds import make_seed
from mechbridge.symbolic_library import parse_explicit, replay


def tautomer():
    """Keto-enol arrows on acetaldehyde (same example as the replay tests)."""
    mol = parse_explicit('[O:1]=[CH:3][CH2:2][H:4]')
    ids = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    arrows = [dict(source=[ids[i] for i in s], sink=[ids[i] for i in t], electrons=2)
              for s, t in [([1], [1, 4]), ([2, 4], [2, 3]), ([1, 3], [1])]]
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    product, edits = replay(mol, arrows)
    proposal = dict(arrows=arrows, edits=edits, template_id='t', predicted_graph=graph_smiles(product))
    atoms = Atoms(numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                  positions=mol.GetConformer().GetPositions())
    return mol, atoms, proposal


def test_chain_rank_linear_and_cyclic():
    linear = [dict(source=[0], sink=[0, 1]), dict(source=[1, 2], sink=[2])]
    rank, succ, cyclic = arrow_chain(linear)
    assert rank == {0: 0, 1: 1} and succ == {0: [1], 1: []} and not cyclic
    ring = [dict(source=[5], sink=[5, 1]), dict(source=[2, 1], sink=[2]),
            dict(source=[2], sink=[2, 14]), dict(source=[5, 14], sink=[5])]
    rank, _, cyclic = arrow_chain(ring)
    assert rank == {0: 0, 1: 1, 2: 2, 3: 3} and cyclic


def test_chain_shift_orders_bonds_and_control_matches_magnitude():
    arrows = [dict(source=[0], sink=[0, 1]), dict(source=[1, 2], sink=[2])]
    edits = [dict(atoms=[0, 1], before=0, after=1), dict(atoms=[1, 2], before=2, after=1)]
    shift, lam = chain_progress_shift(edits, arrows, 0, 11)
    assert lam == -.3 and np.allclose(shift, [.15, -.15])
    control, lam_control = chain_progress_shift(edits, None, 0, 11)
    assert lam_control == lam and sorted(np.abs(control)) == pytest.approx(sorted(np.abs(shift)))
    zero, _ = chain_progress_shift(edits, arrows, 2, 11)       # lambda = 0 on this sample
    assert np.allclose(zero, 0.)


def test_bond_edits_seed_is_independent_of_arrows_under_features():
    mol, atoms, proposal = tautomer()
    scrambled = dict(proposal, arrows=[dict(a, source=a['sink'], sink=a['source'])
                                       for a in reversed(proposal['arrows'])])
    for sample in range(3):
        a = make_seed(atoms, mol, 'bond_edits', proposal, sample, 29, seed_features='arrow_features_v1')
        b = make_seed(atoms, mol, 'bond_edits', scrambled, sample, 29, seed_features='arrow_features_v1')
        assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
        assert a[2]['arrow_terms']['chain_order'] == 'random_bond_permutation_control'


def test_arrow_features_are_active_and_legacy_is_unchanged():
    mol, atoms, proposal = tautomer()
    assert lone_pair_terms(proposal['arrows'], mol, atoms.positions)
    for sample in range(3):
        legacy = make_seed(atoms, mol, 'arrows', proposal, sample, 29)
        explicit = make_seed(atoms, mol, 'arrows', proposal, sample, 29, seed_features='legacy')
        assert np.array_equal(legacy[0], explicit[0]) and 'arrow_terms' not in legacy[2]
        featured = make_seed(atoms, mol, 'arrows', proposal, sample, 29, seed_features='arrow_features_v1')
        terms = featured[2]['arrow_terms']
        assert terms['lone_pair_terms'] >= 1 and terms['chain_order'] == 'arrow_push_pull_rank'
        assert np.linalg.norm(featured[0]-atoms.positions) == pytest.approx(legacy[2]['displacement_norm_A'])
    with pytest.raises(ValueError):
        make_seed(atoms, mol, 'arrows', proposal, 0, 29, seed_features='unknown')


def test_seed_fit_cap_only_changes_fits_that_reach_it():
    mol, atoms, proposal = tautomer()
    for sample in range(3):
        base = make_seed(atoms, mol, 'arrows', proposal, sample, 29, seed_features='arrow_features_v1')
        assert base[2]['geometric_fit_nfev'] < 200
        wider = make_seed(atoms, mol, 'arrows', proposal, sample, 29, seed_features='arrow_features_v1',
                          seed_fit_max_nfev=3000)
        assert np.array_equal(base[0], wider[0]) and np.array_equal(base[1], wider[1])
    with pytest.raises(ValueError, match='did not converge'):
        make_seed(atoms, mol, 'arrows', proposal, 0, 29, seed_features='arrow_features_v1', seed_fit_max_nfev=1)
