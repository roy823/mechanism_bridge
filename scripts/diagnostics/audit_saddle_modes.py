"""Inspect rejected saddle modes; this diagnostic does not refine or accept events."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.reaction_network import CountedCalculator,SearchProtocol,inspect_point,atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('attempts',type=Path,nargs='+')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    calc=CountedCalculator(ReactionPotential(ROOT/'models/aimnet2-rxn'),1000)
    rows=[]
    for folder in a.attempts:
        folder=folder.resolve()
        atoms=read(folder/'ts.xyz');atoms.calc=calc
        point,modes=inspect_point(atoms,SearchProtocol())
        direction=np.load(folder/'seed_direction.npy').ravel()
        negative=[]
        for freq,mode in zip(point['frequencies_cm-1'],modes):
            if freq>=-point['imaginary_threshold_cm-1']:continue
            weights=np.sum(mode**2,axis=1)
            negative.append(dict(frequency_cm=float(freq),
                seed_mode_absolute_overlap=float(abs(np.dot(mode.ravel(),direction))),
                top_displaced_coordinate_indices=np.argsort(-weights)[:4].tolist(),
                squared_cartesian_displacement_fractions=(weights/weights.sum()).tolist()))
        rows.append(dict(source=str(folder.relative_to(ROOT)),point=point,negative_modes=negative))
    result=dict(evaluations=calc.calls,cases=rows,interpretation=
        'Post-hoc mode characterization only. Frequency/overlap alone does not establish a torsional or reaction assignment.',
        changes_to_primary_outcomes=False)
    atomic_json(a.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
