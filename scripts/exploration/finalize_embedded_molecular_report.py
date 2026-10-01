"""Give embedded aggregate events unique IDs and rebuild sampled XYZ exports."""
import argparse
import json
from pathlib import Path
import re
import sys

from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from mechbridge.report_layout import molecular_document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html', type=Path)
    args = parser.parse_args()
    html = args.html.resolve()
    output = html.parent
    text = html.read_text(encoding='utf-8')
    begin = text.index('const DATA=') + len('const DATA=')
    end = text.index(';\n', begin)
    payload = json.loads(text[begin:end])
    by_edge = {}
    for event in payload['events']:
        by_edge[(event['network_id'], event['aggregate_edge_id'])] = event
    selected = []
    for network in payload['networks']:
        prefix = re.sub(r'[^A-Za-z0-9]+', '_', network['id']).strip('_')
        for edge in network['edges']:
            event = by_edge[(network['id'], edge['id'])]
            event_id = f"{prefix}_edge_{edge['id']}"
            event['id'] = event_id
            event['downloads'] = dict(path=event_id + '_path.xyz', ts=event_id + '_TS.xyz')
            edge['event_id'] = event_id
            selected.append(event)
    if len({event['id'] for event in selected}) != len(selected):
        raise ValueError('Final event IDs are not unique')
    for path in list(output.glob('*_path.xyz')) + list(output.glob('*_TS.xyz')):
        if path.parent.resolve() != output:
            raise ValueError('Refusing to remove export outside report directory')
        path.unlink()
    for event in selected:
        trajectory = []
        for frame in event['frames']:
            atoms = Atoms(numbers=event['numbers'], positions=frame['positions'])
            atoms.info['relative_energy_eV'] = frame['relative_energy']
            atoms.info['branch'] = frame['branch']
            trajectory.append(atoms)
        write(output / event['downloads']['path'], trajectory, format='extxyz')
        write(output / event['downloads']['ts'], trajectory[event['ts_index']], write_results=False)
    payload['events'] = selected
    html.write_text(molecular_document(payload, output), encoding='utf-8')
    result = dict(events=len(selected), networks=len(payload['networks']),
        path_xyz=len(list(output.glob('*_path.xyz'))), ts_xyz=len(list(output.glob('*_TS.xyz'))),
        path_policy='sampled frames embedded by the compact report, with TS frame retained')
    (output / 'finalization.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
