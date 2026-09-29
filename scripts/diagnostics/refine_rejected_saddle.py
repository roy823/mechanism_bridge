"""Post-hoc tighter-tolerance search of one rejected candidate; preserves primary results."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.reaction_network import CountedCalculator,SearchProtocol,inspect_point,search_connection,atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('attempt',type=Path)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    if a.outdir.exists():raise FileExistsError(a.outdir)
    a.outdir.mkdir(parents=True)
    protocol=SearchProtocol(fmax=.005,ts_steps=200,descent_steps=300,
                            total_evaluations=2400,evaluations_per_attempt=2200)
    calc=CountedCalculator(ReactionPotential(ROOT/'models/aimnet2-rxn'),protocol.total_evaluations)
    atoms=read(a.attempt/'ts.xyz');atoms.calc=calc
    initial,modes=inspect_point(atoms,protocol)
    if initial['imaginary_count']<2:raise ValueError('Requires a multiply unstable candidate')
    np.save(a.outdir/'initial_modes.npy',modes)
    manifest=dict(source=str(a.attempt),source_sha256=hashlib.sha256((a.attempt/'result.json').read_bytes()).hexdigest(),
        protocol=asdict(protocol),initial=initial,selection='First failed enol continuation candidate, chosen before this refinement',
        changes_to_primary_outcomes=False,source_code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__),ROOT/'src/mechbridge/reaction_network.py',ROOT/'src/mechbridge/aimnet_backend.py',ROOT/'src/mechbridge/physics.py']})
    atomic_json(a.outdir/'manifest.json',manifest)
    calc.attempt_limit=calc.calls+protocol.evaluations_per_attempt
    result=search_connection(atoms,modes[0],calc,a.outdir/'search',protocol)
    result=dict(result,total_evaluations=calc.calls,primary_counts_unchanged=True,
                interpretation='Post-hoc tighter-convergence diagnostic; no DFT validation of this candidate')
    atomic_json(a.outdir/'result.json',result)
    print(json.dumps({k:result[k] for k in ['status','total_evaluations','primary_counts_unchanged']}))


if __name__=='__main__':main()
