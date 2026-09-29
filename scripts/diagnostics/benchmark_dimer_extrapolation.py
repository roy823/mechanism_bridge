"""Paired fixed-seed comparison of ASE's built-in Dimer force extrapolation."""
from dataclasses import asdict
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

CASES={
    'acetone':'reports/local_transfer_v3/search_curvature/acetone_s17/acetone/arrows/attempt_003',
    'enol':'reports/local_transfer_v3/continuation/acetone_enol_observed/arrows/attempt_000',
    'glycolaldehyde':'reports/network_growth_v4/search/glycolaldehyde_s17/glycolaldehyde/arrows/attempt_000',
    'cyclobutanone':'reports/network_growth_v4/search/cyclobutanone_s17/cyclobutanone/arrows/attempt_000'}


def main():
    out=ROOT/'reports/network_growth_v4/dimer_extrapolation'
    out.mkdir(exist_ok=False)
    sources=list((ROOT/'src/mechbridge').glob('*.py'))+[Path(__file__)]
    with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in sources:z.write(p,p.relative_to(ROOT).as_posix())
    atomic_json(out/'manifest.json',dict(cases=CASES,selection='Fixed candidate seeds; engineering microbenchmark, not independent discovery',
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}))
    backend=ReactionPotential(ROOT/'models/aimnet2-rxn');rows=[]
    for index,(name,source) in enumerate(CASES.items()):
        for enabled in ((False,True) if index%2==0 else (True,False)):
            protocol=SearchProtocol(fmax=.005,ts_steps=160,descent_steps=250,hessian_batch_size=32,
                total_evaluations=1600,evaluations_per_attempt=1600,dimer_extrapolate_forces=enabled)
            calc=CountedCalculator(backend,1600)
            seed=read(ROOT/source/'seed.xyz');direction=np.load(ROOT/source/'seed_direction.npy')
            result=search_connection(seed,direction,calc,out/name/str(enabled),protocol)
            row=dict(case=name,extrapolate_forces=enabled,status=result['status'],evaluations=calc.calls,
                model_calls=calc.model_calls,seconds=result['seconds'],
                endpoints=sorted(e['graph_smiles'] for e in result['endpoints']),protocol=asdict(protocol))
            rows.append(row);atomic_json(out/'summary.json',dict(rows=rows))
            print(json.dumps(row),flush=True)
    comparisons=[]
    for name in CASES:
        off=next(r for r in rows if r['case']==name and not r['extrapolate_forces'])
        on=next(r for r in rows if r['case']==name and r['extrapolate_forces'])
        comparisons.append(dict(case=name,evaluations_off=off['evaluations'],evaluations_on=on['evaluations'],
            status_off=off['status'],status_on=on['status'],same_endpoint_graphs=off['endpoints']==on['endpoints'],
            wall_time_ratio=off['seconds']/on['seconds']))
    atomic_json(out/'summary.json',dict(rows=rows,comparisons=comparisons,
        default_enabled=False,interpretation='Validation outcomes take precedence over speed; no automatic promotion from four seeds'))


if __name__=='__main__':main()
