"""Complete a standalone DFT Hessian for a saved endpoint diagnostic geometry."""
import argparse
import json
from pathlib import Path
import sys
import time
from ase.io import read
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.backends import PySCFCalculator
from mechbridge.verification import stationary
from mechbridge.event_graph import geometry_mol,graph_smiles


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path)
    args=p.parse_args();folder=args.folder.resolve()
    source=json.loads((folder.parent/'verification.json').read_text())
    report=json.loads((folder/'result.json').read_text())
    atoms=read(folder/'last_candidate.xyz')
    calc=PySCFCalculator(0,1,source['method'],source['basis'],threads=2)
    started=time.time()
    with threadpool_limits(limits=1,user_api='blas'):
        result,_,_=stationary(atoms,calc,fmax=.002)
    result.update(sign=-1,graph_smiles=graph_smiles(geometry_mol(atoms.numbers,atoms.positions,0)),
        minimum_certified=result['force_converged'] and result['imaginary_count']==0)
    report.update(status='completed_one_direction',branches=[result],
        stop_reason='Interrupted during Hessian after force convergence; completed saved-geometry Hessian separately',
        resumed_hessian_seconds=time.time()-started,unattempted_signs=[1])
    (folder/'result.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
