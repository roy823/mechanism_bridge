"""Check SCF density reuse against independent initial guesses at three geometries."""
import json
from pathlib import Path
import sys
import time
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from mechbridge.backends import PySCFCalculator

atoms=read(ROOT/"reports/feasibility/MR_537711_0/ts.xyz")
warm=PySCFCalculator(0,1,threads=4,reuse_density=True)
rows=[]
for displacement in (0.,.001,.002):
    trial=atoms.copy(); trial.positions[0,0]+=displacement
    start=time.time(); trial.calc=PySCFCalculator(0,1,threads=4)
    e=trial.get_potential_energy(); f=trial.get_forces(); cold_time=time.time()-start
    start=time.time(); trial.calc=warm
    ew=trial.get_potential_energy(); fw=trial.get_forces(); warm_time=time.time()-start
    rows.append({'displacement_A':displacement,'energy_difference_eV':float(abs(e-ew)),
                 'force_max_difference_eV_A':float(abs(f-fw).max()),
                 'cold_seconds':cold_time,'warm_seconds':warm_time})
passed=all(r['energy_difference_eV']<1e-6 and r['force_max_difference_eV_A']<1e-5 for r in rows)
report={'passed':passed,'method':'wb97x','basis':'6-31g(d)','rows':rows}
(ROOT/'reports/density_reuse_check.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
if not passed: raise SystemExit(1)
