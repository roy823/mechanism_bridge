import importlib.util
import io
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]
WATER = [(8, 0., 0., 0.), (1, .757, .586, 0.), (1, -.757, .586, 0.)]
METHANE = [(6, 0., 0., 5.), (1, .629, .629, 5.629), (1, -.629, -.629, 5.629), (1, -.629, .629, 4.371),
           (1, .629, -.629, 4.371)]


def xyz(atoms):
    symbols = {1: 'H', 6: 'C', 8: 'O'}
    return f'{len(atoms)}\n\n' + ''.join(f'{symbols[z]} {x} {y} {w}\n' for z, x, y, w in atoms)


def test_author_ts_is_mapped_into_the_start_order(tmp_path, monkeypatch):
    path = ROOT/'scripts/data/prepare_coley_ts_references.py'
    spec = importlib.util.spec_from_file_location('prepare_coley_ts_references', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    archive = tmp_path/'profiles.tar.gz'
    files = {'r0_CH4.xyz': xyz(METHANE), 'r1_H2O.xyz': xyz(WATER), 'TS_guess.xyz': xyz(METHANE + WATER),
             'TS_imag_mode.xyz': xyz(METHANE + WATER)}
    with tarfile.open(archive, 'w:gz') as tar:
        for name, text in files.items():
            data = text.encode()
            info = tarfile.TarInfo(f'full_dataset_profiles/7/{name}')
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    references = tmp_path/'refs.jsonl'
    references.write_text(json.dumps(dict(id='coley_7', rxn_id='7', atomic_numbers=[8, 1, 1, 6, 1, 1, 1, 1],
                                          reactant_bonds=[[0, 1], [0, 2], [3, 4], [3, 5], [3, 6], [3, 7]],
                                          product_bonds=[[0, 1], [0, 2], [0, 3], [3, 4], [3, 5], [3, 6], [3, 7]],
                                          admission='accepted')) + '\n')
    out = tmp_path/'out.jsonl'
    monkeypatch.setattr('sys.argv', ['prep', '--profiles', str(archive), '--references', str(references),
                                     '--out', str(out)])
    module.main()
    row = json.loads(out.read_text())
    assert row['admission'] == 'accepted', row['ts_source']
    assert row['ts_positions_A'][0] == [0., 0., 0.]           # oxygen: author atom 5 (r1 after methane)
    assert row['ts_positions_A'][3] == [0., 0., 5.]           # carbon: author atom 0
    assert row['ts_forming_bond_distances_A'] == [5.]         # the O-C pair of the product bonds
