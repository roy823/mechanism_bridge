"""Extend prematurely exhausted geometric baselines under the same PES budget."""
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]


def main():
    root=ROOT/'reports/network_growth_v4'
    summary=json.loads((root/'summary.json').read_text())
    selected=[]
    for row in summary['rows']:
        n=json.loads((ROOT/row['file']).read_text(encoding='utf-8'))
        if row['strategy']=='geometry' and n['stop_reason']=='frontier_exhausted' and row['evaluations']<16000:
            selected.append(row['start'])
    out=root/'geometry_budget_check';out.mkdir(exist_ok=False)
    manifest=dict(selection='All geometry runs whose finite proposal pool ended before 16000 evaluations',
        starts=selected,changes='Raise geometric seeds per node from 9 to 24; same reactants, random seed, solver and evaluation cap',
        started_at=time.time(),completed=[])
    def execute(start):
        cmd=[sys.executable,str(ROOT/'scripts/exploration/run_network_exploration.py'),
            '--outdir',str(out/start),'--starts',str(ROOT/'data/processed/network_growth_starts.jsonl'),
            '--start-ids',start,'--strategies','geometry','--attempts','24','--seeds-per-node','1',
            '--evaluations','16000','--attempt-evaluations','1400','--geometry-seeds-per-node','24','--seed','17']
        with (out/(start+'.log')).open('w',encoding='utf-8') as log:
            result=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        return dict(start=start,exit_code=result.returncode)
    with ThreadPoolExecutor(max_workers=2) as executor:
        for future in as_completed([executor.submit(execute,s) for s in selected]):
            r=future.result();manifest['completed'].append(r);print(json.dumps(r),flush=True)
    manifest['finished_at']=time.time()
    (out/'campaign.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    if any(r['exit_code'] for r in manifest['completed']):raise SystemExit(1)


if __name__=='__main__':main()
