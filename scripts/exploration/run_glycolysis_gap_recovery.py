"""Run targeted continuation searches for four missing glycolysis segments."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
STARTS = ROOT / 'data/processed/glycolysis_gap_recovery_starts.jsonl'
PLANS = {
    'pep_enol_to_pyruvate': dict(actions=['grammar:pyruvate_enol_to_keto',
        'grammar:pyruvate_keto_to_enol'], attempts=12, evaluations=20000),
    'glucose_c6': dict(actions=['grammar:glucose_c6_phosphate_condensation',
        'grammar:phosphate_anion_monoester_hydrolysis'], attempts=14, evaluations=28000),
    'pg_shift_soft': dict(actions=['grammar:phosphoglycerate_3_to_2_site_specific'],
        attempts=14, evaluations=26000),
    'pg_shift_strong': dict(actions=['grammar:phosphoglycerate_3_to_2_site_specific'],
        attempts=14, evaluations=26000),
    'fbp_cleavage_strong': dict(actions=['grammar:beta_hydroxy_carbonyl_retro_aldol',
        'grammar:alpha_hydroxy_carbonyl_aldol_addition'], attempts=14, evaluations=28000),
    'fbp_cleavage_very_strong': dict(actions=['grammar:beta_hydroxy_carbonyl_retro_aldol',
        'grammar:alpha_hydroxy_carbonyl_aldol_addition'], attempts=14, evaluations=28000),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ROOT / 'reports/glycolysis_gap_recovery_v18')
    parser.add_argument('--parallel-systems', type=int, default=3)
    parser.add_argument('--workers-per-system', type=int, default=2)
    parser.add_argument('--threads-per-worker', type=int, default=3)
    args = parser.parse_args()
    args.outdir = args.outdir.resolve()
    if args.outdir.exists():
        raise FileExistsError('Choose a fresh campaign directory')
    args.outdir.mkdir(parents=True)
    starts = [json.loads(line) for line in STARTS.read_text(encoding='utf-8').splitlines()]
    jobs = []
    for start in starts:
        group = start['provenance']['recovery_group']
        jobs.append(dict(start=start['id'], group=group,
            seed_scale=start['provenance']['seed_scale'], **PLANS[group]))
    manifest = dict(started_at=time.time(), model='aimnet2-2025 member0', jobs=jobs,
        purpose='targeted recovery of four missing v17 skeleton segments',
        source_products_or_TS_used=False, completed=[], failed=[])
    manifest_path = args.outdir / 'campaign.json'

    def save():
        temporary = manifest_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        temporary.replace(manifest_path)

    def execute(job):
        destination = args.outdir / 'runs' / job['start']
        command = [sys.executable, str(ROOT / 'scripts/exploration/run_network_exploration.py'),
            '--outdir', str(destination), '--potential', 'aimnet2-2025', '--starts', str(STARTS),
            '--start-ids', job['start'], '--strategies', 'arrows', '--template-ids', *job['actions'],
            '--workers', str(args.workers_per_system), '--threads-per-worker', str(args.threads_per_worker),
            '--attempts', str(job['attempts']), '--evaluations', str(job['evaluations']),
            '--attempt-evaluations', '4000', '--ts-steps', '550', '--descent-steps', '500',
            '--initial-steps', '900', '--initial-fmax', '.02', '--fmax', '.01',
            '--endpoint-acceptance', 'carbon_skeleton', '--endpoint-core-fmax', '.02',
            '--hessian-batch-size', '32', '--dimer-extrapolate-forces',
            '--seeds-per-node', '2', '--symbolic-seed-scale', str(job['seed_scale'])]
        log = args.outdir / f"{job['start']}.log"
        with log.open('w', encoding='utf-8') as stream:
            process = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        network = destination / job['start'] / 'arrows/network.json'
        metrics = {}
        if network.exists():
            data = json.loads(network.read_text(encoding='utf-8'))
            metrics = dict(status=data['status'], evaluations=data.get('evaluations'),
                species=len(data.get('species_nodes', [])), edges=len(data['edges']))
        return dict(job, exit_code=process.returncode,
                    network=network.relative_to(ROOT).as_posix(), **metrics)

    save()
    with ThreadPoolExecutor(max_workers=args.parallel_systems) as executor:
        futures = [executor.submit(execute, job) for job in jobs]
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
