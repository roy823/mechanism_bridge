from mechbridge.reference_matching import automorphisms, bond_set, edge_matches, stereo_free


def test_automorphisms_respect_elements_and_bonds():
    # Ethane-like heavy skeleton C0-C1 with three hydrogens on each carbon.
    numbers = [6, 6, 1, 1, 1, 1, 1, 1]
    bonds = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 5], [1, 6], [1, 7]]
    maps = automorphisms(numbers, bonds)
    assert len(maps) == 2 * 6 * 6          # swap carbons, permute each CH3
    react = bond_set(bonds)
    assert all(edge_matches(react, react, react, react, [m]) for m in maps[:5])


def test_stereo_free_key():
    assert stereo_free('C[C@@H](O)C=O') == stereo_free('CC(O)C=O')
