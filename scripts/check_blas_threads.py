"""Measure BLAS oversubscription without changing the SCF initial guess or PES."""
import json
from pathlib import Path
import sys
import time
import numpy as np
from ase.io import read
from threadpoolctl import threadpool_info,threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.backends import PySCFCalculator

atoms=read(ROOT/'reports/feasibility/MR_537711_0/ts.xyz')
initial=threadpool_info()
rows=[]; results=[]
for limit in (None,1):
    with threadpool_limits(limits=limit,user_api='blas'):
        atoms.calc=PySCFCalculator(0,1,threads=2,reuse_density=False)
        start=time.time(); energy=atoms.get_potential_energy(); forces=atoms.get_forces()
        rows.append({'blas_threads':limit,'seconds':time.time()-start})
        results.append((energy,forces))
report={'initial_threadpools':initial,'rows':rows,
        'energy_difference_eV':abs(results[0][0]-results[1][0]),
        'force_max_difference_eV_A':float(abs(results[0][1]-results[1][1]).max())}
report['passed']=report['energy_difference_eV']<1e-6 and report['force_max_difference_eV_A']<1e-5
(ROOT/'reports/blas_threads_check.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
