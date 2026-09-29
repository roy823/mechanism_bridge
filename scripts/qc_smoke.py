#!/usr/bin/env python3
"""Real HF/STO-3G NH3 inversion smoke test; constructed example, not database evidence."""
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from mechbridge.backends import PySCFCalculator
from mechbridge.physics import verify_saddle
from mechbridge.electronic import compute_ibo

out=Path('reports/qc_smoke');out.mkdir(parents=True,exist_ok=True)
positions=[[0.,0.,0.]]+[[np.cos(t),np.sin(t),0.] for t in np.arange(3)*2*np.pi/3]
atoms=Atoms('NH3',positions=positions)
factory=lambda:PySCFCalculator(charge=0,multiplicity=1,method='HF',basis='sto-3g')
atoms.calc=factory()
# Planar symmetry is retained from exactly planar input. Hessian then tests instability.
BFGS(atoms,logfile=str(out/'planar_optimization.log')).run(fmax=.001,steps=100)
write(out/'planar_nh3.xyz',atoms)
result=verify_saddle(atoms,factory,out,fmax=.003,steps=150)
result.update(example_origin='constructed_planar_NH3_not_downloaded_dataset',method='HF/STO-3G',
              purpose='numerical_pipeline_smoke_not_target_level_benchmark')
(out/'verification.json').write_text(json.dumps(result,indent=2))
compute_ibo([atoms],out/'ibo',charge=0,multiplicity=1,method='HF',basis='sto-3g')
print(json.dumps(result,indent=2))
