"""Independently refine an observed MLIP edge and verify its DFT IRC endpoints."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.verification import verify_event


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,default=ROOT/'reports/network_exploration_v1')
    p.add_argument('--start',required=True)
    p.add_argument('--strategy',default='arrows')
    p.add_argument('--edge',type=int,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    args=p.parse_args()
    if args.outdir.exists():raise FileExistsError(args.outdir)
    source=args.run/args.start/args.strategy/'network.json'
    network=json.loads(source.read_text(encoding='utf-8'))
    edge=next(e for e in network['edges'] if e['id']==args.edge)
    attempt=network['attempts'][edge['attempt']]
    trial_file=source.parent/attempt['artifact']
    trial=json.loads(trial_file.read_text(encoding='utf-8'))
    if trial['status']!='validated_descents':raise ValueError('Requires an observed MLIP connection')
    record=dict(event_id=f"{args.start}_{args.strategy}_edge{args.edge}_DFT",
        atomic_numbers=network['start']['atomic_numbers'],charge=0,multiplicity=1,
        positions_A={'ts':trial['ts_positions_A']},
        reactant_smiles=trial['endpoints'][0]['graph_smiles'],
        product_smiles=trial['endpoints'][1]['graph_smiles'],
        expected_pair_origin='Actual MLIP descent endpoints, not the symbolic target',
        source_json=str(trial_file.relative_to(ROOT)),
        source_sha256=hashlib.sha256(trial_file.read_bytes()).hexdigest(),
        original_symbolic_proposal=attempt['proposal'],
        source_index_one_frequencies=trial['ts']['frequencies_cm-1'])
    args.outdir.mkdir(parents=True)
    (args.outdir/'input.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    provenance=dict(packages={m:importlib.metadata.version(m) for m in ['pyscf','sella','numpy','ase']},
        blas_threads=1,pyscf_threads=2,frames=0,
        selection='Two predeclared representative root-connected graph pairs from v1, one expected and one unexpected',
        files={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
               for p in [ROOT/'src/mechbridge/verification.py',ROOT/'src/mechbridge/backends.py',Path(__file__)]})
    (args.outdir/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(json.dumps({'event':record['event_id'],'status':'starting_DFT_and_IRC'}),flush=True)
    with threadpool_limits(limits=1,user_api='blas'):
        result=verify_event(record,args.outdir,threads=2,frames=0,ts_steps=100,irc_steps=160)
    print(json.dumps({k:result.get(k) for k in ['event_id','status','physical_event_verified',
        'expected_endpoint_match','gradient_evaluations','elapsed_seconds','error']}),flush=True)


if __name__=='__main__':main()
