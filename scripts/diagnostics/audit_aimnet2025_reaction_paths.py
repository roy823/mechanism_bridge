"""Audit AIMNet2-2025 formal-algorithm reaction-path artifacts."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'reports/aimnet2025_reaction_paths'


def main():
    summary=json.loads((BASE/'summary.json').read_text(encoding='utf-8'));checks=[]
    unique_pairs={(e['system'],tuple(sorted(e['endpoints']))) for e in summary['events']}
    checks += [dict(check='edge_count',passed=summary['edges']==len(summary['events'])),
               dict(check='unique_graph_pairs',passed=summary['unique_graph_pairs']==len(unique_pairs)),
               dict(check='literature_pairs',passed=summary['unique_literature_pairs']>=5),
               dict(check='formal_registration_contract',passed='even when A is neither endpoint' in summary['registration']),
               dict(check='detached_events_are_retained',passed=any(not e['source_connected'] for e in summary['events']))]
    for campaign in ['intramolecular','targeted_bimolecular','epoxide_ammonia','bimolecular','frontier_expansion']:
        manifest=json.loads((BASE/campaign/'manifest.json').read_text(encoding='utf-8'))
        checks.append(dict(check='aimnet2_2025_manifest',campaign=campaign,passed=manifest['model']=='aimnet2-2025'))
    for event in summary['events']:
        path=ROOT/event['result_file'];result=json.loads(path.read_text(encoding='utf-8'))
        graphs=[endpoint['graph_smiles'] for endpoint in result['endpoints']]
        checks.append(dict(check='event_physics_and_graphs',file=event['result_file'],passed=(
            result['status']=='validated_descents' and result['ts']['imaginary_count']==1 and
            all(e['force_converged'] and e['imaginary_count']==0 and e['barrier_eV']>=0 for e in result['endpoints']) and
            sorted(graphs)==sorted(event['endpoints']))))
    for page in ['intramolecular/molecules/index.html','targeted_bimolecular/molecules/index.html',
                 'epoxide_ammonia/molecules/index.html','bimolecular/molecules/index.html',
                 'frontier_expansion/molecules/index.html','aggregate/molecules/index.html']:
        checks.append(dict(check='molecular_view_exists',page=page,passed=(BASE/page).exists()))
    target_pairs={tuple(sorted(e['endpoints'])) for e in summary['events'] if e['literature_relevant']}
    required=[('C=O.O=CCO','O=C[C@@H](O)CO'),('C=O.O/C=C/O','O=C[C@@H](O)CO'),
              ('C#[N+][N-]C.C=C','CN1CCC=N1'),('C1CO1.N','NCCO'),('C=O.C=O','O=CCO')]
    checks.append(dict(check='five_named_literature_pairs',passed=all(tuple(sorted(pair)) in target_pairs for pair in required)))
    aggregate=summary['aggregate']
    checks += [dict(check='aggregate_raw_edges',passed=aggregate['raw_edges']==summary['edges']),
               dict(check='aggregate_deduplicates_nodes',passed=aggregate['merged_nodes']<aggregate['raw_nodes']),
               dict(check='aggregate_total_edges',passed=aggregate['unique_ts_edges']>=40),
               dict(check='aggregate_inventory_groups',passed=len(aggregate['systems'])==4),
               dict(check='c3_inventory_combines_three_start_systems',passed=any(
                   row['system']=='C3H6O3' and len(row['source_systems'])==3 for row in aggregate['systems']))]
    report=dict(checks=checks,passed=sum(c['passed'] for c in checks),total=len(checks),all_passed=all(c['passed'] for c in checks))
    (BASE/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    if not report['all_passed']:raise SystemExit(1)


if __name__=='__main__':main()
