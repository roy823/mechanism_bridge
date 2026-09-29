"""Replay three fixed near-threshold failures using the original-force stopping rule."""
import hashlib
import json
from pathlib import Path
import sys
import zipfile
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.reaction_network import CountedCalculator,SearchProtocol,search_connection,atomic_json

CASES=[('glycolaldehyde',6),('acetone',1),('acetone',10)]


def main():
    out=ROOT/'reports/network_growth_v4/convergence_contract'
    out.mkdir(exist_ok=False)
    sources=list((ROOT/'src/mechbridge').glob('*.py'))+[Path(__file__)]
    with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in sources:z.write(p,p.relative_to(ROOT).as_posix())
    atomic_json(out/'manifest.json',dict(selection='Three fixed projected-converged but physically unconverged primary candidates; no primary outcomes overwritten',
        cases=CASES,source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}))
    backend=ReactionPotential(ROOT/'models/aimnet2-rxn');rows=[]
    for name,attempt in CASES:
        source=ROOT/f'reports/network_growth_v4/search/{name}_s17/{name}/arrows/attempt_{attempt:03d}'
        old=json.loads((source/'result.json').read_text())
        if old['status']!='ts_force_unconverged' or not old['optimizer_converged']:
            raise ValueError('Diagnostic requires the projected/physical stopping mismatch')
        calc=CountedCalculator(backend,1400)
        protocol=SearchProtocol(hessian_batch_size=32,fmax=.005,total_evaluations=1400,
                                evaluations_per_attempt=1400,dimer_extrapolate_forces=False)
        result=search_connection(read(source/'seed.xyz'),np.load(source/'seed_direction.npy'),calc,
                                 out/f'{name}_{attempt}',protocol)
        rows.append(dict(start=name,attempt=attempt,old_status=old['status'],new_status=result['status'],
            old_evaluations=old['evaluations'],new_evaluations=result['evaluations'],
            old_force=old['ts']['force_max_eV_A'],new_force=result.get('ts',{}).get('force_max_eV_A'),
            endpoints=[e['graph_smiles'] for e in result['endpoints']]))
        atomic_json(out/'summary.json',dict(rows=rows,primary_counts_unchanged=True))
        print(json.dumps(rows[-1]),flush=True)


if __name__=='__main__':main()
