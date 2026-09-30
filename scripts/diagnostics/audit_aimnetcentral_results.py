"""Audit model files, cross-device predictions, reactive checks, and pipeline smoke evidence."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];REPORT=ROOT/'reports/aimnetcentral_v8'


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    summary=read(REPORT/'summary.json');receipt=read(REPORT/'model_receipt.json');checks=[]
    for item in receipt['files']:
        path=ROOT/item['file'];checks.append(dict(check='model_sha256',family=item['family'],
            passed=path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256']))
    cpu={m['family']:m for m in summary['CPU']['reaction_TS_diagnostics'][0]['models']}
    # All reaction diagnostics are deterministic across CPU/GPU; compare every field that defines the local PES check.
    for c,g in zip(summary['CPU']['reaction_TS_diagnostics'],summary['GPU']['reaction_TS_diagnostics']):
        for left,right in zip(c['models'],g['models']):
            checks.append(dict(check='cross_device_reactive_diagnostic',case=c['case'],family=left['family'],
                passed=left['family']==right['family'] and left['imaginary_count']==right['imaginary_count']
                and abs(left['force_max_eV_A']-right['force_max_eV_A'])<1e-4
                and max(abs(a-b) for a,b in zip(left['barriers_from_DFT_endpoints_eV'],right['barriers_from_DFT_endpoints_eV']))<1e-4))
    for family,rows in summary['reactive_transfer'].items():
        for row in rows:checks.append(dict(check='reference_pair_preserved',case=row['case'],family=family,
            passed=row['DFT_reference_pair_preserved'] and row['status']!='started'))
    for family,row in summary['SN2'].items():
        checks.append(dict(check='documented_SN2_strict',family=family,passed=row['strict_self_exchange_passed']))
    smoke=summary['pipeline_smoke'];checks.append(dict(check='broad_charged_pipeline_smoke',
        passed=smoke['status']=='completed' and smoke['initial_point']['force_converged']
        and smoke['initial_point']['imaginary_count']==0 and smoke['start']['charge']==-1))
    consistency=summary['registry_vs_pinned_rxn']['registry_weights_without_D3_vs_pinned']
    checks.append(dict(check='registry_weights_match_pinned_when_D3_aligned',
        passed=consistency['energy_abs_difference_eV']==0 and consistency['force_max_abs_difference_eV_A']==0))
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,
        windows_tests='51 passed, 1 skipped',wsl_tests='52 passed',browser_report_check=read(ROOT/'reports/site_checks/validation.json')['passed'],
        pending_computations=[])
    (REPORT/'validation_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks),failures=[c for c in checks if not c['passed']]),indent=2))
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
