"""Run comparable geometry/net-edit/full-arrow GPA-style searches on real starts."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zipfile
from dataclasses import asdict
from rdkit import RDLogger
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.potentials import load_potential
from mechbridge.symbolic_library import ArrowLibrary
from mechbridge.reaction_network import SearchProtocol, explore, atomic_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--potential',choices=['aimnet2-rxn','aimnet2','aimnet2-2025','aimnet2-nse'],default='aimnet2-rxn')
    p.add_argument('--device',choices=['cpu','cuda'],default='cpu')
    p.add_argument('--compile-model',action='store_true')
    p.add_argument('--strategies', nargs='+', choices=['geometry','center_random','bond_edits','arrows','hybrid','neb_arrows'],
                   default=['geometry','bond_edits','arrows'])
    p.add_argument('--start-ids', nargs='+')
    p.add_argument('--starts',type=Path,default=ROOT/'data/processed/network_starts.jsonl')
    p.add_argument('--attempts', type=int, default=12)
    p.add_argument('--seeds-per-node', type=int, default=1)
    p.add_argument('--hessian-batch-size',type=int,default=32)
    p.add_argument('--fmax',type=float,default=.005)
    p.add_argument('--geometry-seeds-per-node',type=int,default=9)
    p.add_argument('--dimer-extrapolate-forces',action='store_true',
                   help='Enable ASE force extrapolation; optional, validated only on a small fixed-seed benchmark')
    p.add_argument('--evaluations', type=int, default=6000)
    p.add_argument('--attempt-evaluations', type=int, default=700)
    p.add_argument('--seed', type=int, default=17)
    a = p.parse_args()
    if a.outdir.exists():
        raise FileExistsError('Use a new output directory; prior experiments are preserved')
    a.outdir.mkdir(parents=True)
    RDLogger.DisableLog('rdApp.*')
    protocol = SearchProtocol(max_attempts=a.attempts, seeds_per_node=a.seeds_per_node,
        total_evaluations=a.evaluations, evaluations_per_attempt=a.attempt_evaluations,
        random_seed=a.seed,hessian_batch_size=a.hessian_batch_size,fmax=a.fmax,
        geometry_seeds_per_node=a.geometry_seeds_per_node,
        dimer_extrapolate_forces=a.dimer_extrapolate_forces)
    starts = [json.loads(l) for l in a.starts.read_text(encoding='utf-8').splitlines()]
    if a.start_ids:
        starts = [s for s in starts if s['id'] in a.start_ids]
        if {s['id'] for s in starts} != set(a.start_ids):
            raise ValueError('Unknown start ID')
    (a.outdir/'starts.jsonl').write_text(''.join(json.dumps(s)+'\n' for s in starts),encoding='utf-8')
    library = ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    backend,potential=load_potential(a.potential,ROOT,a.device,a.compile_model)
    sources = list((ROOT/'src/mechbridge').glob('*.py')) + [Path(__file__)]
    with zipfile.ZipFile(a.outdir/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for source in sources:archive.write(source,source.relative_to(ROOT).as_posix())
    manifest = dict(protocol=asdict(protocol), starts=[s['id'] for s in starts],
        strategies=a.strategies,**potential,
        source_sha256={s.relative_to(ROOT).as_posix():hashlib.sha256(s.read_bytes()).hexdigest() for s in sources},
        starts_sha256=hashlib.sha256(a.starts.read_bytes()).hexdigest(),
        symbolic_library_sha256=hashlib.sha256((ROOT/'data/raw/synepd/polar.json').read_bytes()).hexdigest(),
        reference_TS_used_in_search=False, reference_product_geometry_used_in_search=False,
        pretrained_overlap=('AIMNet2-rxn includes RGD1; engineering feasibility, not unseen chemistry'
            if a.potential=='aimnet2-rxn' else 'Broad AIMNetCentral training overlap not audited for this run'),
        symbolic_hypotheses=library.policy,
        symbolic_library_audit=dict(library.audit),
        started_at_unix=time.time())
    atomic_json(a.outdir/'manifest.json', manifest)
    for start in starts:
        for strategy in a.strategies:
            report = explore(start, library, backend, strategy, a.outdir/start['id']/strategy, protocol)
            print(json.dumps(dict(start=start['id'], strategy=strategy, status=report['status'],
                                  edges=len(report['edges']), evaluations=report['evaluations'])), flush=True)


if __name__ == '__main__':
    main()
