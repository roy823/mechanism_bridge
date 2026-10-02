"""Rerun historical CPU DFT verifications on the GPU (end-to-end check before Fig. 6).

Each record is the 'source' of a historical verification.json: the same MLIP TS
start, Sella TS refinement, Hessians, Sella IRC and endpoint polish, now with
device='gpu'. The GPU run should reach the same status and endpoint graphs;
TS and endpoint energy differences and wall time are reported against the CPU
result. The historical runs used 2 CPU threads.
Usage: gpu_verification_smoke.py --historical V.json [V.json ...] --outdir DIR [--threads 8]
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.verification import verify_event  # noqa: E402

EV_KCAL = 23.060548


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical', type=Path, nargs='+', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    rows = []
    for path in args.historical:
        old = json.loads(path.read_text(encoding='utf-8'))
        protocol = old['protocol']
        new = verify_event(old['source'], args.outdir/path.parent.name, threads=args.threads, frames=0,
                           ts_steps=protocol['ts_max_steps'], irc_steps=protocol['irc_max_steps'],
                           ts_optimizer_fmax=protocol.get('ts_optimizer_fmax_eV_A', .01), device='gpu')
        row = dict(historical=str(path), status=[old['status'], new['status']],
                   verified=[old.get('physical_event_verified'), new.get('physical_event_verified')],
                   gradient_evaluations=[old.get('gradient_evaluations'), new.get('gradient_evaluations')],
                   elapsed_seconds=[old.get('elapsed_seconds'), new.get('elapsed_seconds')],
                   historical_threads=protocol.get('threads'))
        if 'ts' in old and 'ts' in new:
            row['ts_energy_difference_kcal_mol'] = (new['ts']['energy_eV'] - old['ts']['energy_eV'])*EV_KCAL
        ends_old = {e['graph_smiles']: e['energy_eV'] for e in old.get('endpoints', [])}
        ends_new = {e['graph_smiles']: e['energy_eV'] for e in new.get('endpoints', [])}
        row['same_endpoint_graphs'] = sorted(ends_old) == sorted(ends_new)
        if row['same_endpoint_graphs'] and ends_old:
            row['endpoint_energy_differences_kcal_mol'] = {g: (ends_new[g] - ends_old[g])*EV_KCAL for g in ends_old}
        rows.append(row)
        print(json.dumps(row), flush=True)
    (args.outdir/'gpu_verification_smoke.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
