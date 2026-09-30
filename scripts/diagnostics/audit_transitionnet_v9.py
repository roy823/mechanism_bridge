"""Audit the strict connectivity and aggregate claims of TransitionNet v9."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'reports/aimnet2025_transitionnet_v9'


def main():
    summary=json.loads((BASE/'summary.json').read_text())
    receipt=json.loads((BASE/'data_receipt.json').read_text())
    checks=[]
    archive=ROOT/'data/raw/landscape17/Landscape17.zip'
    checks.append(dict(check='landscape17_size',passed=archive.stat().st_size==receipt['size']))
    checks.append(dict(check='landscape17_md5',passed=hashlib.md5(archive.read_bytes()).hexdigest()==receipt['md5']))
    pes=json.loads((BASE/'landscape17_malonaldehyde_pes_loop/summary.json').read_text())
    checks.append(dict(check='reference_minima_recovered',passed=sum(r['matched'] for r in pes['node_reference_matches'])==2))
    checks.append(dict(check='one_additional_minimum',passed=sum(not r['matched'] for r in pes['node_reference_matches'])==1))
    accepted=[c for c in pes['connections'] if c['accepted']]
    checks.append(dict(check='three_validated_model_edges',passed=len(accepted)==3))
    for connection in accepted:
        selected=[r for r in connection['refinements'] if r['target_pair']]
        checks.append(dict(check='accepted_edge_has_target_pair',job=connection['job'],passed=len(selected)==1 and selected[0]['status']=='validated_descents'))
    for path in [BASE/'landscape17_malonaldehyde_autonomous',BASE/'glyceraldehyde',BASE/'glycolaldehyde']:
        for network_path in path.rglob('network.json'):
            network=json.loads(network_path.read_text())
            checks.append(dict(check='registered_edges_are_source_connected',file=network_path.relative_to(ROOT).as_posix(),
                passed=all(edge.get('source_connected') for edge in network['edges'])))
    landscape=summary['landscape17']
    checks += [dict(check='assisted_count',passed=landscape['assisted_recovered']==2),
               dict(check='exact_reference_edge_count',passed=landscape['exact_reference_edges']==1),
               dict(check='complete_network_claim_is_false',passed=summary['conclusions']['complete_transitionnet_recovered'] is False)]
    report=dict(checks=checks,passed=sum(c['passed'] for c in checks),total=len(checks),all_passed=all(c['passed'] for c in checks))
    (BASE/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    if not report['all_passed']:raise SystemExit(1)


if __name__=='__main__':main()
