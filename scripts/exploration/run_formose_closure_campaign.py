"""Run targeted neutral tetrose <-> 2GO closure searches across conformers."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
STARTS=ROOT/'data/processed/formose_closure_starts.jsonl'
TEMPLATES=['grammar:formose_retro_aldol_tetrose_to_2go',
           'grammar:formose_inverse_aldol_2go_to_tetrose']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,default=ROOT/'reports/formose_closure_v13')
    parser.add_argument('--parallel-starts',type=int,default=3)
    args=parser.parse_args()
    if args.outdir.exists():raise FileExistsError('Choose a fresh output directory')
    args.outdir.mkdir(parents=True)
    starts=[json.loads(line)['id'] for line in STARTS.read_text(encoding='utf-8').splitlines()]
    manifest=dict(started_at=time.time(),model='aimnet2-2025 member0',starts=starts,
        templates=TEMPLATES,parallel_starts=args.parallel_starts,completed=[],failed=[])
    path=args.outdir/'campaign.json'
    def save():
        temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(manifest,indent=2),encoding='utf-8');temporary.replace(path)
    save()
    def execute(start):
        destination=args.outdir/'runs'/start
        command=[sys.executable,str(ROOT/'scripts/exploration/run_network_exploration.py'),
            '--outdir',str(destination),'--potential','aimnet2-2025','--starts',str(STARTS),
            '--start-ids',start,'--strategies','arrows','--template-ids',*TEMPLATES,
            '--workers','2','--threads-per-worker','4','--attempts','6','--evaluations','12000',
            '--attempt-evaluations','2200','--fmax','.005','--hessian-batch-size','32',
            '--dimer-extrapolate-forces']
        log=args.outdir/f'{start}.log'
        with log.open('w',encoding='utf-8') as stream:
            result=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        network_path=destination/start/'arrows/network.json';metrics={}
        if network_path.exists():
            network=json.loads(network_path.read_text(encoding='utf-8'))
            metrics=dict(status=network['status'],attempts=len(network['attempts']),
                evaluations=network.get('evaluations'),nodes=len(network['nodes']),edges=len(network['edges']))
        return dict(start=start,exit_code=result.returncode,
                    network=network_path.relative_to(ROOT).as_posix(),**metrics)
    with ThreadPoolExecutor(max_workers=args.parallel_starts) as executor:
        futures=[executor.submit(execute,start) for start in starts]
        for future in as_completed(futures):
            result=future.result();manifest['completed' if result['exit_code']==0 else 'failed'].append(result)
            save();print(json.dumps(result),flush=True)
    manifest['finished_at']=time.time();save()
    if manifest['failed']:raise SystemExit(1)


if __name__=='__main__':main()
