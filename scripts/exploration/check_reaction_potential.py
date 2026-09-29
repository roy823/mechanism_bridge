"""Check the actual pretrained force interface against finite differences and rotations."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.spatial.transform import Rotation
from ase import Atoms
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.reaction_network import atomic_json


def main():
    calc=ReactionPotential(ROOT/'models/aimnet2-rxn')
    starts=[json.loads(l) for l in (ROOT/'data/processed/network_starts.jsonl').read_text().splitlines()]
    results=[]
    rotation=Rotation.from_euler('xyz',[.2,.4,-.3]).as_matrix()
    for start in starts:
        a=Atoms(numbers=start['atomic_numbers'],positions=start['positions_A'])
        a.calc=calc
        energy=a.get_potential_energy(); forces=a.get_forces()
        x=a.positions.copy()
        errors=[]
        for j in range(3):
            delta=np.zeros_like(x);delta[0,j]=.001
            a.positions=x+delta;ep=a.get_potential_energy()
            a.positions=x-delta;em=a.get_potential_energy()
            errors.append(abs(-(ep-em)/.002-forces[0,j]))
        a.positions=x@rotation+np.array([1.,-2.,.5])
        er=a.get_potential_energy();fr=a.get_forces()
        row=dict(start=start['id'],force_finite_difference_max_error_eV_A=max(errors),
                 rotation_energy_error_eV=abs(er-energy),
                 rotation_force_error_eV_A=float(abs(fr-forces@rotation).max()))
        row['passed']=max(errors)<.01 and abs(er-energy)<.001 and row['rotation_force_error_eV_A']<.001
        results.append(row)
    report=dict(protocol=dict(finite_difference_A=.001, force_tolerance_eV_A=.01,
                rotation_energy_tolerance_eV=.001, rotation_force_tolerance_eV_A=.001),
                results=results,passed=all(r['passed'] for r in results),
                scope='Numerical consistency, not physical accuracy vs DFT')
    atomic_json(ROOT/'reports/network_potential_check.json',report)
    print(json.dumps(report,indent=2))
    if not report['passed']: raise SystemExit(1)


if __name__=='__main__':main()
