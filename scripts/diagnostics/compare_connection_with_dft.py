"""P2b: MLIP endpoint protocols versus DFT IRC endpoints on DFT-verified events.

For each verified event, the stored TS is refined to an MLIP index-one saddle
with Sella P-RFO, then both connection protocols (+/- mode displacement with
BFGS, and bidirectional MLIP IRC with BFGS polish) are run from that saddle.
The stereo-free endpoint pair of each protocol is compared with the DFT IRC
endpoint pair recorded in verification.json (pre-registered P2 criterion).
"""
import argparse
import dataclasses
import json
from pathlib import Path
import sys
import time

import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.potentials import load_potential  # noqa: E402
from mechbridge.reaction_network import (CountedCalculator, SearchProtocol, inspect_point,  # noqa: E402
                                         integrate_irc)
from mechbridge.reference_matching import stereo_free  # noqa: E402


def endpoints(ts, modes, calculator, protocol, outdir, charge):
    keys, extra = [], []
    for sign, direction in ((-1, 'reverse'), (1, 'forward')):
        end = ts.copy()
        end.calc = calculator
        if protocol.connection_protocol == 'irc':
            extra.append(integrate_irc(end, outdir, protocol, direction)['irc_converged'])
        else:
            end.positions += sign*protocol.mode_displacement*modes[0]
        with BFGS(end, maxstep=.1, logfile=str(outdir/f'polish_{direction}.log')) as opt:
            opt.run(fmax=protocol.fmax, steps=protocol.descent_steps)
        analysis, _ = inspect_point(end, protocol)
        keys.append(dict(key=stereo_free(graph_smiles(geometry_mol(end.numbers, end.positions, charge))),
                         minimum=bool(analysis['force_converged'] and analysis['imaginary_count'] == 0),
                         energy_eV=analysis['energy_eV']))
    return keys, extra


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--potential', default='aimnet2-rxn')
    parser.add_argument('--prfo-steps', type=int, default=200)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    backend, provenance = load_potential(args.potential, ROOT, threads=2)
    base = SearchProtocol(hessian_batch_size=32)
    rows = []
    for path in sorted(ROOT.glob('reports/**/verification.json')):
        record = json.loads(path.read_text(encoding='utf-8'))
        source = record.get('source') or {}
        if not record.get('physical_event_verified') or 'ts' not in (source.get('positions_A') or {}):
            continue
        dft = sorted(stereo_free(e['graph_smiles']) for e in record['endpoints'])
        numbers, charge = source['atomic_numbers'], source['charge']
        backend.validate_system(numbers, charge, source['multiplicity'])
        row = dict(event=path.parent.relative_to(ROOT).as_posix(), dft_endpoints=dft, atoms=len(numbers))
        started = time.time()
        calculator = CountedCalculator(backend, 100000)
        outdir = args.out/path.parent.name
        outdir.mkdir(parents=True, exist_ok=False)
        ts = Atoms(numbers=numbers, positions=source['positions_A']['ts'])
        ts.calc = calculator
        from sella import Sella
        with Sella(ts, order=1, internal=False, logfile=str(outdir/'prfo.log')) as opt:
            row['prfo_converged'] = bool(opt.run(fmax=base.fmax, steps=args.prfo_steps))
        analysis, modes = inspect_point(ts, base)
        row.update(mlip_index_one=bool(analysis['force_converged'] and analysis['imaginary_count'] == 1),
                   ts_evaluations=calculator.calls)
        if row['mlip_index_one']:
            for name in ('mode_displacement', 'irc'):
                protocol = dataclasses.replace(base, connection_protocol=name)
                calls = calculator.calls
                folder = outdir/name
                folder.mkdir()
                keys, irc_converged = endpoints(ts, modes, calculator, protocol, folder, charge)
                pair = sorted(k['key'] for k in keys)
                row[name] = dict(endpoints=pair, matches_dft=pair == dft,
                                 both_minima=all(k['minimum'] for k in keys),
                                 irc_converged=irc_converged or None,
                                 evaluations=calculator.calls - calls)
        row['seconds'] = time.time() - started
        rows.append(row)
        print(json.dumps(row), flush=True)
    summary = dict(potential=provenance['model'], events=len(rows),
                   index_one=sum(r['mlip_index_one'] for r in rows),
                   displacement_matches_dft=sum(r.get('mode_displacement', {}).get('matches_dft', False)
                                                for r in rows),
                   irc_matches_dft=sum(r.get('irc', {}).get('matches_dft', False) for r in rows), rows=rows)
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
