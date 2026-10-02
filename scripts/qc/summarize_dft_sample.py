"""Fig. 6: DFT verification of the pre-registered MLIP edge sample, per stratum.

Reads the sample (sample_dft_edges.py) and the verification directories written
by uf/dft_verify_array.sbatch (edge_NNN/verification.json, index = sample row).
Per stratum and overall:
  confirmed  : physical_event_verified (DFT TS with one imaginary mode, both IRC
               branches converged to distinct minima, non-negative barriers);
  same ends  : the DFT IRC endpoint graphs equal the MLIP endpoint graphs
               (expected_endpoint_match), among confirmed edges;
  not certified : any other final status (DFT non-convergence included), which
               is never evidence that the reaction does not exist;
  barrier error : |barrier(MLIP) - barrier(DFT)| for both directions of confirmed
               edges with matching ends, in kcal/mol.
Rates carry Wilson 95% intervals. Rows whose verification has not finished are
listed as pending and left out of every rate; with --all-finished (every SLURM
task has ended, e.g. some hit the time limit) they count as not certified with
status stopped_<last status>.
Usage: summarize_dft_sample.py SAMPLE.json DFT_DIR --out FILE [--all-finished]
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reference_matching import wilson  # noqa: E402

EV_KCAL = 23.060548
FINAL = ('physical_event_verified', 'physical_event_verified_with_symbolic_hypotheses', 'connection_not_certified',
         'unresolved_ts', 'unresolved_irc', 'calculation_failed')


def mlip_barriers(row):
    """Forward and reverse MLIP barriers (eV) of the sampled edge, keyed by endpoint graph."""
    network = json.loads((Path(row['run'])/row['start']/row['strategy']/'network.json').read_text(encoding='utf-8'))
    edge = next(e for e in network['edges'] if e['id'] == row['edge'])
    result = json.loads((Path(row['run'])/row['start']/row['strategy']/
                         network['attempts'][edge['attempt']]['artifact']).read_text(encoding='utf-8'))
    ts = result['ts']['energy_eV']
    return {e['graph_smiles']: ts - e['energy_eV'] for e in result['endpoints']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample', type=Path)
    parser.add_argument('dft', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--all-finished', action='store_true')
    args = parser.parse_args()
    sample = json.loads(args.sample.read_text(encoding='utf-8'))['sample']
    rows, pending = [], []
    for k, row in enumerate(sample):
        path = args.dft/f'edge_{k:03d}'/'verification.json'
        if not path.exists():
            pending.append(k)
            continue
        result = json.loads(path.read_text(encoding='utf-8'))
        if result.get('status') not in FINAL:
            if not args.all_finished:
                pending.append(k)
                continue
            result = dict(result, status=f"stopped_{result.get('status')}", physical_event_verified=False)
        out = dict(index=k, stratum=row['stratum'], start=row['start'], status=result['status'],
                   confirmed=bool(result.get('physical_event_verified')),
                   same_ends=bool(result.get('expected_endpoint_match')),
                   gradient_evaluations=result.get('gradient_evaluations'), seconds=result.get('elapsed_seconds'))
        if out['confirmed'] and out['same_ends']:
            mlip = mlip_barriers(row)
            dft = {e['graph_smiles']: e['barrier_electronic_eV'] for e in result['endpoints']}
            out['barrier_errors_kcal_mol'] = [abs(mlip[g] - dft[g])*EV_KCAL for g in dft if g in mlip]
        rows.append(out)
    def block(group):
        confirmed = sum(r['confirmed'] for r in group)
        same = sum(r['confirmed'] and r['same_ends'] for r in group)
        errors = [x for r in group for x in r.get('barrier_errors_kcal_mol', [])]
        return dict(edges=len(group), confirmed=confirmed, confirmation_rate_wilson95=wilson(confirmed, len(group)),
                    same_ends_among_confirmed=same, not_certified=len(group) - confirmed,
                    barrier_mae_kcal_mol=float(np.mean(errors)) if errors else None,
                    barrier_errors=len(errors))
    strata = sorted({r['stratum'] for r in rows})
    summary = dict(sample=str(args.sample), dft=str(args.dft), finished=len(rows), pending=pending,
                   overall=block(rows), strata={s: block([r for r in rows if r['stratum'] == s]) for s in strata},
                   rows=rows)
    args.out.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
