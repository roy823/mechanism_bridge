"""Run comparable geometry/net-edit/full-arrow GPA-style searches on real starts."""
import argparse
import dataclasses
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
from mechbridge.symbolic_library import ArrowLibrary,FilteredArrowLibrary
from mechbridge.reaction_network import (ENCOUNTER_POLICIES, SearchProtocol, explore,
                                         atomic_json)
from mechbridge.parallel_network import START_METHODS, explore_shared
from mechbridge.provenance import runtime_provenance


# CLI flag destination -> SearchProtocol field.
PROTOCOL_FLAGS = dict(attempts='max_attempts', seeds_per_node='seeds_per_node',
    symbolic_seed_scale='symbolic_seed_scale', evaluations='total_evaluations',
    attempt_evaluations='evaluations_per_attempt', ts_steps='ts_steps',
    descent_steps='descent_steps', hessian_batch_size='hessian_batch_size', fmax='fmax',
    initial_fmax='initial_fmax', initial_steps='initial_steps',
    endpoint_acceptance='endpoint_acceptance', endpoint_core_fmax='endpoint_core_fmax',
    geometry_seeds_per_node='geometry_seeds_per_node',
    dimer_extrapolate_forces='dimer_extrapolate_forces', min_barrier='min_barrier_eV',
    encounter_policy='encounter_policy')


def protocol_from_json(path, seed):
    """Frozen protocol: every SearchProtocol field except random_seed, which is per run."""
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    fields = {f.name for f in dataclasses.fields(SearchProtocol)} - {'random_seed'}
    missing, unknown = sorted(fields - set(data)), sorted(set(data) - fields)
    if missing or unknown:
        raise ValueError(f'Protocol JSON must list exactly the SearchProtocol fields except '
                         f'random_seed; missing={missing}, unknown={unknown}')
    return SearchProtocol(**data, random_seed=seed)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--potential',choices=['aimnet2-rxn','aimnet2','aimnet2-2025','aimnet2-nse'],default='aimnet2-rxn')
    p.add_argument('--device',choices=['cpu','cuda'],default='cpu')
    p.add_argument('--compile-model',action='store_true')
    p.add_argument('--workers',type=int,default=1,
                   help='CPU process workers sharing one centrally registered TransitionNet')
    p.add_argument('--threads-per-worker',type=int,default=2)
    p.add_argument('--start-method',choices=START_METHODS,default='spawn',
                   help='Process start method for --workers > 1 (recorded in network.json)')
    p.add_argument('--strategies', nargs='+', choices=['geometry','center_random','bond_edits','arrows','hybrid'],
                   default=['geometry','bond_edits','arrows'])
    p.add_argument('--template-ids',nargs='+',help='Restrict symbolic proposals to exact reviewed template IDs')
    p.add_argument('--start-ids', nargs='+')
    p.add_argument('--starts',type=Path,default=ROOT/'data/processed/network_starts.jsonl')
    p.add_argument('--attempts', type=int, default=12)
    p.add_argument('--ts-steps',type=int,default=160)
    p.add_argument('--descent-steps',type=int,default=250)
    p.add_argument('--seeds-per-node', type=int, default=1)
    p.add_argument('--symbolic-seed-scale',type=float,default=1.0)
    p.add_argument('--hessian-batch-size',type=int,default=32)
    p.add_argument('--fmax',type=float,default=.005)
    p.add_argument('--initial-fmax',type=float,default=.003)
    p.add_argument('--initial-steps',type=int,default=250)
    p.add_argument('--endpoint-acceptance',choices=['full_system','carbon_skeleton'],default='full_system')
    p.add_argument('--endpoint-core-fmax',type=float,default=.02)
    p.add_argument('--geometry-seeds-per-node',type=int,default=9)
    p.add_argument('--dimer-extrapolate-forces',action='store_true',
                   help='Enable ASE force extrapolation; optional, validated only on a small fixed-seed benchmark')
    p.add_argument('--min-barrier',type=float,default=-1e-4,
                   help='Required E_TS - E_endpoint in eV for both endpoints (legacy -1e-4)')
    p.add_argument('--encounter-policy',choices=ENCOUNTER_POLICIES,default='symbolic_only',
                   help="'matched_controls' gives geometry/center_random the same rigid encounter search")
    p.add_argument('--evaluations', type=int, default=6000)
    p.add_argument('--attempt-evaluations', type=int, default=700)
    p.add_argument('--seed', type=int, default=17)
    p.add_argument('--protocol-json',type=Path,
                   help='Frozen protocol (all SearchProtocol fields except random_seed); '
                        'protocol flags above must then stay at their defaults')
    a = p.parse_args()
    if a.protocol_json:
        changed = sorted('--'+k.replace('_','-') for k in PROTOCOL_FLAGS if getattr(a,k)!=p.get_default(k))
        if changed:
            p.error(f'--protocol-json cannot be combined with protocol flags: {changed}')
    if a.workers<1 or a.threads_per_worker<1 or a.seeds_per_node<1 or a.symbolic_seed_scale<=0:
        p.error('workers, threads-per-worker, seeds-per-node and symbolic-seed-scale must be positive')
    protocol = protocol_from_json(a.protocol_json, a.seed) if a.protocol_json else SearchProtocol(
        max_attempts=a.attempts, seeds_per_node=a.seeds_per_node,
        symbolic_seed_scale=a.symbolic_seed_scale,
        total_evaluations=a.evaluations, evaluations_per_attempt=a.attempt_evaluations,
        ts_steps=a.ts_steps,descent_steps=a.descent_steps,
        random_seed=a.seed,hessian_batch_size=a.hessian_batch_size,fmax=a.fmax,
        initial_fmax=a.initial_fmax,initial_steps=a.initial_steps,
        endpoint_acceptance=a.endpoint_acceptance,endpoint_core_fmax=a.endpoint_core_fmax,
        geometry_seeds_per_node=a.geometry_seeds_per_node,
        dimer_extrapolate_forces=a.dimer_extrapolate_forces,
        min_barrier_eV=a.min_barrier,encounter_policy=a.encounter_policy)
    if a.outdir.exists():
        raise FileExistsError('Use a new output directory; prior experiments are preserved')
    a.outdir.mkdir(parents=True)
    RDLogger.DisableLog('rdApp.*')
    starts = [json.loads(l) for l in a.starts.read_text(encoding='utf-8').splitlines()]
    if a.start_ids:
        starts = [s for s in starts if s['id'] in a.start_ids]
        if {s['id'] for s in starts} != set(a.start_ids):
            raise ValueError('Unknown start ID')
    (a.outdir/'starts.jsonl').write_text(''.join(json.dumps(s)+'\n' for s in starts),encoding='utf-8')
    library = ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    if a.template_ids:library=FilteredArrowLibrary(library,a.template_ids)
    backend,potential=load_potential(a.potential,ROOT,a.device,a.compile_model,
                                     threads=a.threads_per_worker)
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
        template_ids=a.template_ids,
        scheduler=('central_species_registry_process_workers' if a.workers>1 else 'sequential'),
        workers=a.workers,threads_per_worker=a.threads_per_worker,
        start_method=a.start_method if a.workers>1 else None,
        protocol_json=(dict(path=str(a.protocol_json),
            sha256=hashlib.sha256(a.protocol_json.read_bytes()).hexdigest()) if a.protocol_json else None),
        provenance=runtime_provenance(ROOT),
        started_at_unix=time.time())
    atomic_json(a.outdir/'manifest.json', manifest)
    for start in starts:
        for strategy in a.strategies:
            if a.workers>1:
                report=explore_shared(start,library,backend,strategy,a.outdir/start['id']/strategy,
                    a.potential,ROOT,protocol,a.workers,a.threads_per_worker,a.device,a.compile_model,
                    a.start_method)
            else:
                report=explore(start,library,backend,strategy,a.outdir/start['id']/strategy,protocol)
            print(json.dumps(dict(start=start['id'], strategy=strategy, status=report['status'],
                                  edges=len(report['edges']), evaluations=report['evaluations'])), flush=True)


if __name__ == '__main__':
    main()
