import importlib.util
import io
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]
WATER = [(8, 0., 0., 0.), (1, .757, .586, 0.), (1, -.757, .586, 0.)]
H2 = [(1, 0., 0., 5.), (1, 0., 0., 5.74)]


def xyz(atoms):
    symbols = {1: 'H', 8: 'O'}
    return f'{len(atoms)}\n\n' + ''.join(f'{symbols[z]} {x} {y} {w}\n' for z, x, y, w in atoms)


def test_author_ts_is_mapped_into_the_start_order(tmp_path, monkeypatch):
    path = ROOT/'scripts/data/prepare_coley_ts_references.py'
    spec = importlib.util.spec_from_file_location('prepare_coley_ts_references', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    archive = tmp_path/'profiles.tar.gz'
    files = {'r0_H2.xyz': xyz(H2), 'r1_H2O.xyz': xyz(WATER), 'TS_guess.xyz': xyz(H2 + WATER),
             'TS_imag_mode.xyz': xyz(H2 + WATER)}
    with tarfile.open(archive, 'w:gz') as tar:
        for name, text in files.items():
            data = text.encode()
            info = tarfile.TarInfo(f'full_dataset_profiles/7/{name}')
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    references = tmp_path/'refs.jsonl'
    references.write_text(json.dumps(dict(id='coley_7', rxn_id='7', atomic_numbers=[8, 1, 1, 1, 1],
                                          reactant_bonds=[[0, 1], [0, 2], [3, 4]])) + '\n')
    out = tmp_path/'out.jsonl'
    monkeypatch.setattr('sys.argv', ['prep', '--profiles', str(archive), '--references', str(references),
                                     '--out', str(out)])
    module.main()
    row = json.loads(out.read_text())
    assert row['ts_positions_A'][0] == [0., 0., 0.]                       # the oxygen, author atom 2
    assert sorted(row['ts_positions_A'][3:]) == [[0., 0., 5.], [0., 0., 5.74]]   # H2 hydrogens
