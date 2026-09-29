"""Post-hoc endpoint curvature diagnosis; never overwrites the primary IRC verdict."""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from ase.io import read,write
from ase.optimize import BFGS
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.backends import PySCFCalculator
from mechbridge.verification import stationary
from mechbridge.event_graph import geometry_mol,graph_smiles


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('verification',type=Path)
    p.add_argument('--branch',choices=['forward','reverse'],required=True)
    args=p.parse_args()
    folder=args.verification.resolve()
    source=json.loads((folder/'verification.json').read_text())
    output=folder/('curvature_recovery_'+args.branch)
    output.mkdir(exist_ok=False)
    start=time.time()
    calc=PySCFCalculator(0,1,source['method'],source['basis'],threads=2)
    atoms=read(folder/('minimum_'+args.branch+'.xyz'))
    report=dict(status='running',selection='Post-hoc diagnosis after primary endpoint failed curvature test',
        original_status=source['status'],original_result_unchanged=True,displacement_A=.08,
        force_tolerance_eV_A=.002,branches=[])
    def save():
        report['seconds']=time.time()-start;report['gradients']=calc.evaluation_count
        (output/'result.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    with threadpool_limits(limits=1,user_api='blas'):
        initial,_,modes=stationary(atoms,calc,fmax=.002)
        report['initial']=initial;save()
        if initial['imaginary_count']==0:
            raise ValueError('This diagnosis requires an endpoint with negative curvature')
        for sign in (-1,1):
            candidate=atoms.copy();candidate.positions+=sign*.08*modes[0];candidate.calc=calc
            with BFGS(candidate,maxstep=.08,logfile=str(output/f'branch_{sign}.log'),
                      trajectory=str(output/f'branch_{sign}.traj')) as opt:
                opt.run(fmax=.002,steps=100)
            result,_,_=stationary(candidate,calc,fmax=.002)
            result['sign']=sign
            result['graph_smiles']=graph_smiles(geometry_mol(candidate.numbers,candidate.positions,0))
            result['minimum_certified']=result['force_converged'] and result['imaginary_count']==0
            write(output/f'minimum_{sign}.xyz',candidate,write_results=False)
            report['branches'].append(result);save()
    report['status']='completed';save()
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
