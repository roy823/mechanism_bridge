import dataclasses
import importlib.util
import json
from pathlib import Path

import pytest

from mechbridge.reaction_network import SearchProtocol

ROOT = Path(__file__).resolve().parents[1]


def load_runner():
    path = ROOT/'scripts/exploration/run_network_exploration.py'
    spec = importlib.util.spec_from_file_location('run_network_exploration', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_protocol_json_lists_every_field_except_seed(tmp_path):
    runner = load_runner()
    full = {k: v for k, v in dataclasses.asdict(SearchProtocol()).items() if k != 'random_seed'}
    path = tmp_path/'protocol.json'
    path.write_text(json.dumps(full))
    assert runner.protocol_from_json(path, 29) == SearchProtocol(random_seed=29)
    for broken in ({k: v for k, v in full.items() if k != 'fmax'},
                   dict(full, random_seed=17), dict(full, unknown_field=1)):
        path.write_text(json.dumps(broken))
        with pytest.raises(ValueError, match='exactly the SearchProtocol fields'):
            runner.protocol_from_json(path, 29)


def test_cli_protocol_flags_map_to_protocol_fields():
    runner = load_runner()
    fields = {f.name for f in dataclasses.fields(SearchProtocol)}
    assert set(runner.PROTOCOL_FLAGS.values()) <= fields


def test_seed_fit_cap_must_be_positive():
    with pytest.raises(ValueError, match='seed_fit_max_nfev'):
        SearchProtocol(seed_fit_max_nfev=0)
