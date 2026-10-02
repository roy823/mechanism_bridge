import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_scorer():
    path = ROOT/'scripts/exploration/score_reference_recovery.py'
    spec = importlib.util.spec_from_file_location('score_reference_recovery', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_mapped_match_allows_only_reactant_automorphisms():
    scorer = load_scorer()
    # Water-like A: O0 bonded to H1 and H2; product moves H2 onto a second O3.
    numbers = [8, 1, 1, 8]
    react = scorer.bond_set([[0, 1], [0, 2]])
    prod = scorer.bond_set([[0, 1], [2, 3]])
    maps = scorer.automorphisms(numbers, [[0, 1], [0, 2]])
    assert len(maps) == 2          # the two hydrogens are equivalent in the reactant
    swapped = scorer.bond_set([[0, 2], [1, 3]])   # the other hydrogen transferred
    assert scorer.edge_matches(react, prod, react, prod, maps)
    assert scorer.edge_matches(swapped, react, react, prod, maps)
    other = scorer.bond_set([[0, 1], [0, 2], [0, 3]])
    assert not scorer.edge_matches(react, other, react, prod, maps)


def test_wilson_interval():
    scorer = load_scorer()
    low, high = scorer.wilson(8, 10)
    assert 0.44 < low < 0.5 and 0.94 < high < 0.98
    assert scorer.wilson(0, 0) is None
