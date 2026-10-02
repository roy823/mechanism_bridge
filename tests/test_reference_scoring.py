import importlib.util
from pathlib import Path

import mechbridge.reference_matching as scorer_module

ROOT = Path(__file__).resolve().parents[1]


def load_scorer():
    return scorer_module


def test_scorer_script_imports():
    path = ROOT/'scripts/exploration/score_reference_recovery.py'
    spec = importlib.util.spec_from_file_location('score_reference_recovery', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.edge_matches is scorer_module.edge_matches


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


def test_cluster_bootstrap_resamples_starts_not_rows():
    path = ROOT/'scripts/exploration/score_reference_recovery.py'
    spec = importlib.util.spec_from_file_location('score_reference_recovery', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # One start with ten hits and one with ten misses: the interval spans both clusters.
    rows = [dict(start='a', S1_mapped=True)]*10 + [dict(start='b', S1_mapped=False)]*10
    low, high = module.cluster_bootstrap(rows, 'S1_mapped', 2000, 1)
    assert low == 0. and high == 1.
    assert module.cluster_bootstrap([], 'S1_mapped', 10, 1) is None


def test_t1x_reactant_groups_share_one_geometry(tmp_path, monkeypatch):
    import json
    path = ROOT/'scripts/data/group_t1x_reactants.py'
    spec = importlib.util.spec_from_file_location('group_t1x_reactants', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    starts = tmp_path/'starts.jsonl'
    rows = [dict(id='r2', atomic_numbers=[1, 1], positions_A=[[0, 0, 0], [0, 0, .74]]),
            dict(id='r1', atomic_numbers=[1, 1], positions_A=[[0, 0, 0], [0, 0, .74]]),
            dict(id='r3', atomic_numbers=[1, 1], positions_A=[[0, 0, 0], [0, 0, .75]]),
            dict(id='r4', atomic_numbers=[1, 1], positions_A=[[0, 0, 0], [0, 0, .76]])]
    starts.write_text(''.join(json.dumps(r) + '\n' for r in rows))
    ids = tmp_path/'ids.json'
    ids.write_text(json.dumps(dict(representable_ids=['r1', 'r2', 'r3'])))
    out = tmp_path/'groups.json'
    monkeypatch.setattr('sys.argv', ['group', '--starts', str(starts), '--ids', str(ids), '--out', str(out)])
    module.main()
    result = json.loads(out.read_text())
    assert result['ids'] == ['r1', 'r3'] and result['groups'] == {'r1': ['r1', 'r2'], 'r3': ['r3']}
