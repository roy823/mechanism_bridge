"""Run a fixed, manifest-recorded experiment plan with bounded independent workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=ROOT/'configs/validation_v2.json')
    p.add_argument('--outdir',type=Path,default=ROOT/'reports/validation_v2/search')
    p.add_argument('--potential',choices=['aimnet2-rxn','aimnet2','aimnet2-2025','aimnet2-nse'],default='aimnet2-rxn')
    p.add_argument('--device',choices=['cpu','cuda'],default='cpu')
    p.add_argument('--compile-model',action='store_true')
    a=p.parse_args()
    if a.outdir.exists():raise FileExistsError('Choose a fresh experiment directory')
    a.outdir.mkdir(parents=True)
    config=json.loads(a.config.read_text())
    plan=[]
    for cohort in ('core','extension'):
        for start in config[cohort+'_starts']:
            for seed in config[cohort+'_random_seeds']:
                plan.append(dict(start=start,seed=seed,cohort=cohort,folder=f'{start}_s{seed}'))
    manifest=dict(config=config,config_sha256=hashlib.sha256(a.config.read_bytes()).hexdigest(),
        potential=a.potential,device=a.device,compile_model=a.compile_model,
        started_at=time.time(),plan=plan,completed=[],failed=[])
    path=a.outdir/'campaign.json'
    def save():
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(manifest,indent=2),encoding='utf-8');tmp.replace(path)
    save()
    def execute(job):
        command=[sys.executable,str(ROOT/'scripts/exploration/run_network_exploration.py'),
            '--outdir',str(a.outdir/job['folder']),'--starts',str(ROOT/config['starts_file']),
            '--start-ids',job['start'],'--seed',str(job['seed']),
            '--attempts',str(config['max_attempts']),'--seeds-per-node',str(config['seeds_per_node']),
            '--evaluations',str(config['total_evaluations']),'--attempt-evaluations',str(config['evaluations_per_attempt']),
            '--potential',a.potential,'--device',a.device,'--strategies',*config['strategies']]
        if a.compile_model:command.append('--compile-model')
        with (a.outdir/(job['folder']+'.log')).open('w',encoding='utf-8') as log:
            result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        return dict(**job,exit_code=result.returncode)
    with ThreadPoolExecutor(max_workers=config['workers']) as executor:
        futures=[executor.submit(execute,job) for job in plan]
        for future in as_completed(futures):
            job=future.result();manifest['completed' if job['exit_code']==0 else 'failed'].append(job)
            save();print(json.dumps(job),flush=True)
    manifest['finished_at']=time.time();save()
    if manifest['failed']:raise SystemExit(1)


if __name__=='__main__':main()
