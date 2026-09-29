"""Check serial/batched finite differences on minimum, TS and rejected saddle."""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.physics import finite_hessian,vibrational_analysis


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    backend=ReactionPotential(ROOT/'models/aimnet2-rxn')
    files=['reports/local_transfer_v3/qc/acetone/minimum_forward.xyz',
           'reports/local_transfer_v3/qc/acetone/ts.xyz',
           'reports/local_transfer_v3/continuation/acetone_enol_observed/arrows/attempt_000/ts.xyz']
    rows=[]
    for source in files:
        a=read(ROOT/source);a.calc=backend
        serial_force=a.get_forces().copy()
        batch_force=backend.evaluate_many(a.numbers,np.stack([a.positions,a.positions]))['forces'][0]
        timings={1:[],32:[]};matrices={}
        for repeat in range(3):
            for size in ((1,32) if repeat%2==0 else (32,1)):
                start=time.perf_counter()
                h=finite_hessian(a,.005,None if size==1 else
                    lambda z,x:backend.evaluate_many(z,x)['forces'],size)
                timings[size].append(time.perf_counter()-start);matrices[size]=h
        vib={size:vibrational_analysis(a.positions,a.get_masses(),h) for size,h in matrices.items()}
        error=float(np.max(np.abs(matrices[1]-matrices[32])))
        ferror=float(np.max(np.abs(serial_force-batch_force)))
        rows.append(dict(source=source,serial_median_seconds=float(np.median(timings[1])),
            batch_median_seconds=float(np.median(timings[32])),speedup=float(np.median(timings[1])/np.median(timings[32])),
            force_max_difference_eV_A=ferror,hessian_max_difference_eV_A2=error,
            imaginary_counts={size:v['imaginary_count'] for size,v in vib.items()},
            passed=bool(error<.01 and ferror<1e-4 and vib[1]['imaginary_count']==vib[32]['imaginary_count'])))
    result=dict(cases=rows,passed=all(r['passed'] for r in rows),device='cpu',threads=2,
        batch_size=32,step_A=.005,repeats=3,interpretation='Component benchmark, not overall search or DFT speedup; concurrent system load may affect timing')
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
