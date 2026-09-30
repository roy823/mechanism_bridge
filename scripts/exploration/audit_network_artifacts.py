"""Audit completed search budgets, endpoint provenance, and source-code integrity."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('--source-archive',type=Path,help='Explicit archived source for a historical run')
    args=p.parse_args()
    manifest=json.loads((args.run/'manifest.json').read_text())
    checks=[]
    archive_path=args.source_archive or args.run/'source_snapshot.zip'
    with zipfile.ZipFile(archive_path) as archive:
        for relative,sha in manifest['source_sha256'].items():
            checks.append(dict(check='archived_source_hash',path=relative,
                passed=relative in archive.namelist() and hashlib.sha256(archive.read(relative)).hexdigest()==sha))
    starts=[json.loads((args.run/s/manifest['strategies'][0]/'network.json').read_text(encoding='utf-8'))['start']
            for s in manifest['starts']]
    for start in starts:
        checks.append(dict(check='starting_input_no_reference_geometry',start=start['id'],
            passed=set(start)=={'id','atomic_numbers','positions_A','charge','multiplicity','provenance'}
            and isinstance(start['positions_A'],list)))
    for start in manifest['starts']:
        for strategy in manifest['strategies']:
            folder=args.run/start/strategy
            n=json.loads((folder/'network.json').read_text(encoding='utf-8'))
            checks.append(dict(check='run_finished',start=start,strategy=strategy,
                passed=n['status'] in ('completed','initial_minimum_unresolved','initialization_budget_exhausted')))
            initialization=n.get('initialization_evaluations',n['evaluations'] if not n['attempts'] else -1)
            calls=initialization+sum(a['evaluations'] for a in n['attempts'])
            checks.append(dict(check='evaluation_accounting',start=start,strategy=strategy,
                passed=calls==n['evaluations'] and calls<=n['protocol']['total_evaluations']
                and all(a['evaluations']<=n['protocol']['evaluations_per_attempt'] for a in n['attempts'])))
            for edge in n['edges']:
                trial=n['attempts'][edge['attempt']]
                detail=json.loads((folder/trial['artifact']).read_text())
                observed=sorted(e['graph_smiles'] for e in detail['endpoints'])
                registered=sorted(n['nodes'][i]['graph_smiles'] for i in edge['nodes'])
                checks.append(dict(check='edge_actual_endpoints',start=start,strategy=strategy,edge=edge['id'],
                    passed=observed==registered and edge['nodes'][0]!=edge['nodes'][1]
                    and detail['status']=='validated_descents' and detail['is_IRC'] is False
                    and detail['DFT_verified'] is False))
                checks.append(dict(check='stationary_point_gates',start=start,strategy=strategy,edge=edge['id'],
                    passed=detail['ts']['imaginary_count']==1 and detail['ts']['force_converged']
                    and detail['ts']['force_max_eV_A']<=n['protocol']['fmax']
                    and all(e['imaginary_count']==0 and e['force_converged'] and
                            e['force_max_eV_A']<=n['protocol']['fmax'] and e['barrier_eV']>=-1e-4
                            for e in detail['endpoints'])))
    result=dict(passed=all(c['passed'] for c in checks),checks=checks)
    (args.run/'artifact_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks),failures=[c for c in checks if not c['passed']]),indent=2))
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
