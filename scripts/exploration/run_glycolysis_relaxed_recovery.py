"""Retry missing glycolysis segments with carbon-skeleton endpoint acceptance."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
STARTS = ROOT / 'data/processed/glycolysis_starts.jsonl'

TEMPLATES = {
    'hmp_g6p': ['SynEPD:1531:forward', 'SynEPD:1536:reverse'],
    'hmp_f6p': ['SynEPD:1531:reverse', 'SynEPD:1536:forward'],
    'triose_gap': ['SynEPD:1762:forward', 'SynEPD:1763:reverse'],
    'triose_dhap': ['SynEPD:1762:reverse', 'SynEPD:1763:forward'],
    'pg_shift_3pg': ['grammar:vicinal_phosphate_migration'],
    'pg_shift_2pg': ['grammar:vicinal_phosphate_migration'],
    'pep_hydrolysis_pep_water': ['grammar:phosphate_monoester_hydrolysis',
        'SynEPD:17:forward', 'SynEPD:17:reverse', 'SynEPD:1762:forward',
        'SynEPD:1762:reverse', 'SynEPD:1763:forward', 'SynEPD:1763:reverse'],
    'pep_hydrolysis_pyruvate_pi': ['grammar:phosphate_monoester_condensation',
        'SynEPD:17:forward', 'SynEPD:17:reverse', 'SynEPD:1762:forward',
        'SynEPD:1762:reverse', 'SynEPD:1763:forward', 'SynEPD:1763:reverse'],
    'glucose_phosphorylation_glucose_pi': ['grammar:phosphate_monoester_condensation'],
    'glucose_phosphorylation_g6p_water': ['grammar:phosphate_monoester_hydrolysis'],
    'f6p_phosphorylation_f6p_pi': ['grammar:phosphate_monoester_condensation'],
    'f6p_phosphorylation_fbp_water': ['grammar:phosphate_monoester_hydrolysis'],
    'bpg_hydrolysis_bpg_water': ['grammar:phosphate_monoester_hydrolysis'],
    'bpg_hydrolysis_3pg_pi': ['grammar:phosphate_monoester_condensation'],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path,
        default=ROOT / 'reports/glycolysis_relaxed_recovery_v15')
    parser.add_argument('--parallel-systems', type=int, default=3)
    parser.add_argument('--workers-per-system', type=int, default=2)
    parser.add_argument('--threads-per-worker', type=int, default=3)
    args = parser.parse_args()
    args.outdir = args.outdir.resolve()
    if args.outdir.exists():
        raise FileExistsError('Choose a fresh campaign directory')
    args.outdir.mkdir(parents=True)
    starts = [json.loads(line) for line in STARTS.read_text(encoding='utf-8').splitlines()]
    plan = [dict(start=row['id'], system=row['provenance']['system'],
                 templates=TEMPLATES[row['provenance']['system']])
            for row in starts if row['provenance']['system'] in TEMPLATES]
    manifest = dict(started_at=time.time(), model='aimnet2-2025 member0', plan=plan,
        endpoint_acceptance='all carbon-containing fragments; force max <= 0.02 eV/A',
        TS_acceptance='full-system force max <= 0.01 eV/A and exactly one imaginary mode',
        completed=[], failed=[])
    manifest_path = args.outdir / 'campaign.json'

    def save():
        temporary = manifest_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        temporary.replace(manifest_path)

    def execute(job):
        destination = args.outdir / 'runs' / job['start']
        command = [sys.executable, str(ROOT / 'scripts/exploration/run_network_exploration.py'),
            '--outdir', str(destination), '--potential', 'aimnet2-2025', '--starts', str(STARTS),
            '--start-ids', job['start'], '--strategies', 'arrows', '--template-ids', *job['templates'],
            '--workers', str(args.workers_per_system), '--threads-per-worker', str(args.threads_per_worker),
            '--attempts', '6', '--evaluations', '18000', '--attempt-evaluations', '4000',
            '--ts-steps', '500', '--descent-steps', '400', '--initial-steps', '700',
            '--initial-fmax', '.02', '--fmax', '.01', '--endpoint-acceptance', 'carbon_skeleton',
            '--endpoint-core-fmax', '.02', '--hessian-batch-size', '32',
            '--dimer-extrapolate-forces']
        log = args.outdir / f"{job['start']}.log"
        with log.open('w', encoding='utf-8') as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        network = destination / job['start'] / 'arrows/network.json'
        metrics = {}
        if network.exists():
            data = json.loads(network.read_text(encoding='utf-8'))
            metrics = dict(status=data['status'], evaluations=data.get('evaluations'),
                physical_nodes=len(data['nodes']), species=len(data.get('species_nodes', [])),
                edges=len(data['edges']))
        return dict(job, exit_code=result.returncode,
                    network=network.relative_to(ROOT).as_posix(), **metrics)

    save()
    with ThreadPoolExecutor(max_workers=args.parallel_systems) as executor:
        futures = [executor.submit(execute, job) for job in plan]
        for future in as_completed(futures):
            result = future.result()
            manifest['completed' if result['exit_code'] == 0 else 'failed'].append(result)
            save()
            print(json.dumps(result), flush=True)
    manifest['finished_at'] = time.time()
    save()
    if manifest['failed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
