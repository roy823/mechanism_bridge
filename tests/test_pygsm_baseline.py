import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_driver():
    path = ROOT/'scripts/baselines/run_pygsm_baseline.py'
    spec = importlib.util.spec_from_file_location('run_pygsm_baseline', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_proposal_driving_sets_dedupe_and_skip_order_changes():
    driver = load_driver()
    proposals = [dict(edits=[dict(atoms=[3, 1], before=0, after=1), dict(atoms=[1, 2], before=2, after=1)]),
                 dict(edits=[dict(atoms=[1, 3], before=0, after=1)]),
                 dict(edits=[dict(atoms=[0, 2], before=1, after=0)]),
                 dict(edits=[dict(atoms=[1, 2], before=2, after=1)])]
    assert driver.proposal_driving_sets(proposals) == [[('ADD', 1, 3)], [('BREAK', 0, 2)]]


def test_b2f2_sampler_is_seeded_and_respects_degree_ceiling():
    driver = load_driver()
    numbers = np.array([6, 8, 1, 1, 8, 1, 1])          # formaldehyde + water
    positions = np.array([[0., 0, 0], [1.2, 0, 0], [-.5, .9, 0], [-.5, -.9, 0],
                          [0., 0, 2.8], [.8, 0, 3.3], [-.8, 0, 3.3]])
    bonds = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6)]
    first = driver.b2f2_driving_sets(numbers, positions, bonds, np.random.default_rng(5), 10)
    again = driver.b2f2_driving_sets(numbers, positions, bonds, np.random.default_rng(5), 10)
    assert first == again and 0 < len(first) <= 10
    assert len({tuple(s) for s in first}) == len(first)
    for coords in first:
        degree = {i: 0 for i in range(len(numbers))}
        for i, j in bonds:
            degree[i] += 1
            degree[j] += 1
        for op, i, j in coords:
            delta = 1 if op == 'ADD' else -1
            degree[i] += delta
            degree[j] += delta
        assert all(degree[a] <= driver.MAX_DEGREE[int(numbers[a])] for a in degree)
