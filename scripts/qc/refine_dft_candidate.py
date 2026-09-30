"""Fresh, separately recorded stricter optimization of an unresolved DFT TS."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile
from ase.io import read
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.verification import verify_event


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verification',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    if a.outdir.exists():raise FileExistsError(a.outdir)
    source=json.loads(a.verification.read_text())
    if source['status']!='unresolved_ts':raise ValueError('Requires unresolved TS evidence')
    record=dict(source['source'])
    record['event_id']+='_stricter_optimizer'
    record['positions_A']={'ts':read(a.verification.parent/'ts.xyz').positions.tolist()}
    record['refinement_parent']=dict(file=a.verification.relative_to(ROOT) if a.verification.is_absolute()
        else a.verification.as_posix(),sha256=hashlib.sha256(a.verification.read_bytes()).hexdigest(),
        reason='Optimizer stopping criterion and final Cartesian force gate differed; physical gate unchanged',
        posthoc=True)
    record['refinement_parent']['file']=str(record['refinement_parent']['file'])
    a.outdir.mkdir(parents=True)
    (a.outdir/'input.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    files=[Path(__file__),ROOT/'src/mechbridge/verification.py',ROOT/'src/mechbridge/backends.py']
    with zipfile.ZipFile(a.outdir/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for f in files:archive.write(f,f.relative_to(ROOT).as_posix())
    with threadpool_limits(limits=1,user_api='blas'):
        result=verify_event(record,a.outdir,method=source['method'],basis=source['basis'],
            threads=2,frames=0,ts_steps=100,irc_steps=160,ts_optimizer_fmax=.01)
    print(json.dumps({k:result.get(k) for k in ['status','physical_event_verified','expected_endpoint_match',
        'gradient_evaluations','elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
