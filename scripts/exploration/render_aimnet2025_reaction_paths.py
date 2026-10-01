"""Render actual saved AIMNet2-2025 TS and descent trajectories without interpolation."""
import json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.molecular_visuals import render_visuals

BASE=ROOT/'reports/aimnet2025_reaction_paths'


def main():
    # Interactive HTML already contains every saved frame.  Remove obsolete
    # static duplicates from prior local renders to keep the artifact compact.
    for path in BASE.rglob('molecules/*'):
        if path.suffix.lower() in ('.png','.pdf'):path.unlink()
    outputs=[]
    for name in ('intramolecular','targeted_bimolecular','epoxide_ammonia'):
        outputs.append(render_visuals(BASE/name,make_static=False))
    run=BASE/'bimolecular'
    paths=sorted((run/'formaldehyde_glycolaldehyde_o0').glob('*/network.json'))
    outputs.append(render_visuals(run,paths,make_static=False))
    run=BASE/'frontier_expansion'
    paths=[path for path in sorted(run.rglob('network.json'))
           if json.loads(path.read_text(encoding='utf-8'))['status']=='completed']
    if paths:outputs.append(render_visuals(run,paths,make_static=False))
    print(json.dumps(outputs,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
