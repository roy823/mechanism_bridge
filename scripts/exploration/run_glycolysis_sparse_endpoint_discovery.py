"""Discover nine glycolysis skeleton segments from one endpoint per inventory."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
STARTS = ROOT / 'data/processed/glycolysis_sparse_endpoint_starts.jsonl'
TAUTOMER = ['grammar:alpha_hydroxy_carbonyl_to_enediol',
    'grammar:enediol_to_left_carbonyl', 'grammar:enediol_to_right_carbonyl']
PLAN = {
    'hexose_isomerization': dict(actions=TAUTOMER, attempts=12, evaluations=24000),
    'fbp_cleavage': dict(actions=['grammar:beta_hydroxy_carbonyl_retro_aldol',
        'grammar:alpha_hydroxy_carbonyl_aldol_addition'], attempts=10, evaluations=22000),
    'triose_isomerization': dict(actions=TAUTOMER, attempts=10, evaluations=20000),
    'phosphoglycerate_shift': dict(actions=['grammar:vicinal_phosphate_anion_migration'],
        attempts=8, evaluations=18000),
    'phosphoglycerate_dehydration': dict(actions=['grammar:phosphoglycerate_dehydration',
        'grammar:phosphoenolpyruvate_hydration'], attempts=8, evaluations=18000),
    'pep_hydrolysis': dict(actions=['grammar:phosphate_anion_monoester_hydrolysis'],
        attempts=10, evaluations=22000),
    'glucose_phosphorylation': dict(actions=['grammar:phosphate_anion_monoester_condensation',
        'grammar:phosphate_anion_monoester_hydrolysis', *TAUTOMER], attempts=12, evaluations=26000),
    'f6p_phosphorylation': dict(actions=['grammar:fbp_c1_phosphate_hydrolysis_to_f6p',
        'grammar:f6p_c1_phosphate_condensation_to_fbp',
        'grammar:metaphosphate_water_hydration_to_h2po4', *TAUTOMER], attempts=14, evaluations=30000),
    'bpg_hydrolysis': dict(actions=['grammar:phosphate_anion_monoester_hydrolysis',
        'grammar:phosphate_anion_monoester_condensation'], attempts=10, evaluations=22000),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path,
        default=ROOT / 'reports/glycolysis_sparse_endpoint_v17')
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
        system = start['provenance']['system']
        jobs.append(dict(start=start['id'], system=system, charge=start['charge'], **PLAN[system]))
    manifest = dict(started_at=time.time(), model='aimnet2-2025 member0', starts_file=STARTS.relative_to(ROOT).as_posix(),
        design='nine fixed inventories; one chemical endpoint root per inventory; two conformer/orientation variants',
        target_products_available_to_search=False, product_or_TS_geometry_used=False,
        physical_scope='charged closed-shell CHNOP gas phase; one stoichiometric Pi or water only where present in root',
        jobs=jobs, completed=[], failed=[])
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
            '--attempt-evaluations', '4000', '--ts-steps', '500', '--descent-steps', '450',
            '--initial-steps', '800', '--initial-fmax', '.02', '--fmax', '.01',
            '--endpoint-acceptance', 'carbon_skeleton', '--endpoint-core-fmax', '.02',
            '--hessian-batch-size', '32', '--dimer-extrapolate-forces']
        log = args.outdir / f"{job['start']}.log"
        with log.open('w', encoding='utf-8') as stream:
            process = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        network = destination / job['start'] / 'arrows/network.json'
        metrics = {}
        if network.exists():
            data = json.loads(network.read_text(encoding='utf-8'))
            metrics = dict(status=data['status'], evaluations=data.get('evaluations'),
                physical_nodes=len(data['nodes']), species=len(data.get('species_nodes', [])),
                edges=len(data['edges']))
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
