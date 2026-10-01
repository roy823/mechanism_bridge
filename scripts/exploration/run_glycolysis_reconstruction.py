"""Run the enzyme-free, all-metabolite glycolysis reconstruction screen."""
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
    'cleavage_fbp': ['grammar:beta_hydroxy_carbonyl_retro_aldol'],
    'cleavage_gap_dhap': ['grammar:alpha_hydroxy_carbonyl_aldol_addition'],
    'triose_gap': ['SynEPD:1762:forward', 'SynEPD:1763:reverse'],
    'triose_dhap': ['SynEPD:1762:reverse', 'SynEPD:1763:forward'],
    'pg_shift_3pg': ['grammar:vicinal_phosphate_migration'],
    'pg_shift_2pg': ['grammar:vicinal_phosphate_migration'],
    'dehydration_2pg': ['grammar:phosphoglycerate_dehydration'],
    'dehydration_pep_water': ['grammar:phosphoenolpyruvate_hydration'],
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
    parser.add_argument('--outdir', type=Path, default=ROOT / 'reports/glycolysis_reconstruction_v14')
    parser.add_argument('--parallel-systems', type=int, default=3)
    parser.add_argument('--workers-per-system', type=int, default=2)
    parser.add_argument('--threads-per-worker', type=int, default=3)
    parser.add_argument('--attempts', type=int, default=6)
    parser.add_argument('--evaluations', type=int, default=18000)
    args = parser.parse_args()
    args.outdir = args.outdir.resolve()
    if args.outdir.exists():
        raise FileExistsError('Choose a fresh campaign directory')
    args.outdir.mkdir(parents=True)
    starts = [json.loads(line) for line in STARTS.read_text(encoding='utf-8').splitlines()]
    plan = []
    for start in starts:
        system = start['provenance']['system']
        plan.append(dict(start=start['id'], system=system, templates=TEMPLATES[system]))
    manifest = dict(started_at=time.time(), model='aimnet2-2025 member0',
        starts_file=STARTS.relative_to(ROOT).as_posix(), plan=plan,
        physical_scope='neutral fully protonated phosphates; no enzyme, ATP/ADP, NAD+/NADH, metal, or solvent',
        parallel_systems=args.parallel_systems, workers_per_system=args.workers_per_system,
        threads_per_worker=args.threads_per_worker, attempts_per_start=args.attempts,
        evaluations_per_start=args.evaluations, completed=[], failed=[])
    manifest_path = args.outdir / 'campaign.json'

    def save():
        temporary = manifest_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        temporary.replace(manifest_path)

    def execute(job):
        destination = args.outdir / 'runs' / job['start']
        command = [sys.executable, str(ROOT / 'scripts/exploration/run_network_exploration.py'),
            '--outdir', str(destination), '--potential', 'aimnet2-2025',
            '--starts', str(STARTS), '--start-ids', job['start'], '--strategies', 'arrows',
            '--template-ids', *job['templates'], '--workers', str(args.workers_per_system),
            '--threads-per-worker', str(args.threads_per_worker), '--attempts', str(args.attempts),
            '--evaluations', str(args.evaluations), '--attempt-evaluations', '4200', '--ts-steps', '400',
            '--descent-steps', '800',
            '--initial-steps', '1200', '--initial-fmax', '.005', '--fmax', '.005',
            '--hessian-batch-size', '32', '--dimer-extrapolate-forces']
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
