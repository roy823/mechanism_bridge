"""Run parallel C2/C3/C4 shared nets covering the canonical formose cycle."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
STARTS=ROOT/'data/processed/formose_cycle_starts.jsonl'
PLAN=[
    ('c2_two_formaldehyde_o0',30),
    ('c3_go_fa_o0',36),
    ('c3_glyceraldehyde_d_o0',36),
    ('c3_dha_o0',36),
    ('c3_enediol_o0',36),
    ('c4_ga_fa_o0',36),
    ('c4_erythrulose_d_o0',30),
    ('c4_two_glycolaldehyde_o0',36),
]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,default=ROOT/'reports/formose_cycle_v12')
    parser.add_argument('--parallel-systems',type=int,default=3)
    parser.add_argument('--workers-per-system',type=int,default=2)
    parser.add_argument('--threads-per-worker',type=int,default=4)
    args=parser.parse_args()
    if args.outdir.exists():raise FileExistsError('Choose a fresh campaign directory')
    args.outdir.mkdir(parents=True)
    starts={json.loads(line)['id'] for line in STARTS.read_text(encoding='utf-8').splitlines()}
    missing={start for start,_ in PLAN}-starts
    if missing:raise ValueError(f'Missing starts: {sorted(missing)}')
    manifest=dict(started_at=time.time(),starts_file=STARTS.relative_to(ROOT).as_posix(),
        plan=[dict(start=start,attempts=attempts) for start,attempts in PLAN],
        parallel_systems=args.parallel_systems,workers_per_system=args.workers_per_system,
        threads_per_worker=args.threads_per_worker,completed=[],failed=[])
    manifest_path=args.outdir/'campaign.json'
    def save():
        temporary=manifest_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(manifest,indent=2),encoding='utf-8');temporary.replace(manifest_path)
    save()
    def execute(job):
        start,attempts=job;destination=args.outdir/'runs'/start
        command=[sys.executable,str(ROOT/'scripts/exploration/run_network_exploration.py'),
            '--outdir',str(destination),'--potential','aimnet2-2025','--starts',str(STARTS),
            '--start-ids',start,'--strategies','hybrid','--workers',str(args.workers_per_system),
            '--threads-per-worker',str(args.threads_per_worker),'--attempts',str(attempts),
            '--evaluations',str(attempts*1200),'--attempt-evaluations','1600','--fmax','.005',
            '--hessian-batch-size','32','--geometry-seeds-per-node','4','--dimer-extrapolate-forces']
        log=args.outdir/f'{start}.log'
        with log.open('w',encoding='utf-8') as stream:
            result=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        network=destination/start/'hybrid/network.json'
        metrics={}
        if network.exists():
            data=json.loads(network.read_text(encoding='utf-8'))
            metrics=dict(status=data['status'],evaluations=data.get('evaluations'),
                         physical_nodes=len(data['nodes']),species=len(data.get('species_nodes',[])),
                         edges=len(data['edges']))
        return dict(start=start,attempts=attempts,exit_code=result.returncode,
                    network=network.relative_to(ROOT).as_posix(),**metrics)
    with ThreadPoolExecutor(max_workers=args.parallel_systems) as executor:
        futures=[executor.submit(execute,job) for job in PLAN]
        for future in as_completed(futures):
            result=future.result();key='completed' if result['exit_code']==0 else 'failed'
            manifest[key].append(result);save();print(json.dumps(result),flush=True)
    manifest['finished_at']=time.time();save()
    if manifest['failed']:raise SystemExit(1)


if __name__=='__main__':main()
