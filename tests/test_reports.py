import json
from pathlib import Path
from mechbridge.symbolic_library import parse_explicit,replay
from mechbridge.intermolecular_actions import electron_actions
from mechbridge.event_graph import graph_smiles


def test_epoxide_hypotheses_preserve_atoms_and_distinguish_symbolic_steps():
    mol=parse_explicit('N.C1CO1');products=set()
    for name,arrows in electron_actions(mol):
        product,_=replay(mol,arrows)
        assert product.GetNumAtoms()==mol.GetNumAtoms()
        products.add(graph_smiles(product))
    assert products=={'[NH3+]CC[O-]','NCCO'}


def test_quantum_report_distinguishes_original_alternative_and_unresolved(tmp_path,monkeypatch):
    from mechbridge import report_catalog
    monkeypatch.setattr(report_catalog,'ROOT',tmp_path)
    for i,(passed,matched) in enumerate([(True,True),(True,False),(False,False)]):
        path=tmp_path/'reports'/'stage'/str(i);path.mkdir(parents=True)
        (path/'verification.json').write_text(json.dumps(dict(event_id=str(i),status='test',physical_event_verified=passed,
            expected_endpoint_match=matched,method='wb97x',basis='6-31g(d)',endpoints=[],source={})))
    checks=report_catalog.quantum_checks()
    assert [q['label'] for q in checks]==['原图对通过','替代图对有效','未通过 / 未完成']


def test_molecular_document_uses_shared_assets_and_safe_embedded_json(tmp_path):
    from mechbridge.report_layout import molecular_document
    payload=dict(events=[],networks=[],figures=[],report_url='index.html',example='</script><script>')
    text=molecular_document(payload,tmp_path)
    assert '__MOLECULAR_DATA__' not in text and '__SITE_CSS__' not in text
    assert '<\\/script><script>' in text
    assert 'aria-label="全站导航"' in text
    assert '该连接未做 DFT / IRC 验证' not in text


def test_reference_target_requires_both_actual_endpoints():
    from mechbridge.event_classification import matches_reference_pair
    r=parse_explicit('C=O.C=O');p=parse_explicit('O=CCO');other=parse_explicit('CO.[C-]#[O+]')
    assert matches_reference_pair([r,p],r,[p])
    assert matches_reference_pair([p,r],r,[p])
    assert not matches_reference_pair([other,p],r,[p])
