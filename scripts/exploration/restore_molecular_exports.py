"""Restore compact XYZ path exports from a self-contained molecular HTML report."""
import argparse
import json
from pathlib import Path

from ase import Atoms
from ase.io import write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html', type=Path)
    args = parser.parse_args()
    html = args.html.resolve()
    text = html.read_text(encoding='utf-8')
    start = text.index('const DATA=') + len('const DATA=')
    payload, _ = json.JSONDecoder().raw_decode(text[start:])
    output = html.parent
    expected={f"{event['id']}{suffix}" for event in payload['events']
              for suffix in ('_path.xyz','_TS.xyz')}
    for pattern in ('*_path.xyz','*_TS.xyz'):
        for old in output.glob(pattern):
            if old.name not in expected:
                old.unlink()
    for event in payload['events']:
        trajectory = []
        for frame in event['frames']:
            atoms = Atoms(numbers=event['numbers'], positions=frame['positions'])
            atoms.info['relative_energy_eV'] = frame['relative_energy']
            atoms.info['branch'] = frame['branch']
            trajectory.append(atoms)
        write(output / f"{event['id']}_path.xyz", trajectory, format='extxyz')
        write(output / f"{event['id']}_TS.xyz", trajectory[event['ts_index']], write_results=False)
    print(json.dumps(dict(events=len(payload['events']), paths=len(list(output.glob('*_path.xyz'))),
                              transition_states=len(list(output.glob('*_TS.xyz'))))))


if __name__ == '__main__':
    main()
