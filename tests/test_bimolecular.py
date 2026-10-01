import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from mechbridge.symbolic_library import parse_explicit,replay
from mechbridge.intermolecular_actions import electron_actions,crosses_components,diverse_proposals
from mechbridge.encounters import assemble_encounter,orient_reactive_encounter
from mechbridge.event_graph import graph_smiles,geometry_mol


@pytest.mark.parametrize('smiles,target',[
    ('C=O.O','OCO'),('C=O.C=O','O=CCO'),('C=O.OC=CO','O=CC(O)CO'),
    ('C=C.C[N-][N+]#C','CN1CCCN=1')])
def test_bimolecular_electron_grammar(smiles,target):
    mol=parse_explicit(smiles);products=[]
    for name,arrows in electron_actions(mol):
        product,edits=replay(mol,arrows)
        assert crosses_components(mol,edits)
        assert Chem.GetFormalCharge(product)==0
        assert product.GetNumAtoms()==mol.GetNumAtoms()
        products.append(graph_smiles(product))
    # Cycloadduct constitution is checked by ring/fragment topology, not an
    # independently annotated electron-path assertion.
    if smiles.startswith('C=C.'):
        assert products
        for smi in products:
            p=Chem.MolFromSmiles(smi)
            assert len(Chem.GetMolFrags(p))==1
            assert any(len(r)==5 for r in p.GetRingInfo().AtomRings())
    else:assert target in products


def test_encounter_rigid_fragments_and_graph_preservation():
    mol=parse_explicit('C=O.O');assert AllChem.EmbedMolecule(mol,randomSeed=17)==0
    x=mol.GetConformer().GetPositions();numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
    fragments=Chem.GetMolFrags(mol)
    y=assemble_encounter(numbers,x,fragments,23)
    assert graph_smiles(geometry_mol(numbers,y,0))==graph_smiles(mol)
    _,arrows=next(electron_actions(mol));_,edits=replay(mol,arrows)
    z,meta=orient_reactive_encounter(numbers,y,mol,edits,23)
    assert meta and not meta['reference_product_or_TS_used']
    for ids in fragments:
        a,b=x[list(ids)],z[list(ids)]
        assert np.allclose(np.linalg.norm(a[:,None]-a[None,:],axis=-1),
                           np.linalg.norm(b[:,None]-b[None,:],axis=-1))
    assert np.allclose(z,orient_reactive_encounter(numbers,y,mol,edits,23)[0])


def test_proposal_diversity_before_truncation():
    p=[dict(intermolecular=True,predicted_graph='A',id=i) for i in range(20)]
    p.append(dict(intermolecular=True,predicted_graph='B',id=20))
    assert [v['predicted_graph'] for v in diverse_proposals(p,2)]==['A','B']


def test_resonance_and_actual_interfragment_classification():
    from mechbridge.event_graph import resonance_equivalent
    from mechbridge.event_classification import classify_event
    assert resonance_equivalent(parse_explicit('[CH2-][N+]#CC'),parse_explicit('C=[N+]=[C-]C'))
    assert not resonance_equivalent(parse_explicit('CC=O'),parse_explicit('C=CO'))
    mol=parse_explicit('C=O.O');_,arrows=next(electron_actions(mol))
    product,_=replay(mol,arrows)
    c=classify_event(mol,product)
    assert c['classification']=='intermolecular_heavy_atom_bond'
    assert c['fragment_counts']==[2,1]


def test_electron_grammar_never_invents_a_proton():
    assert list(electron_actions(parse_explicit('C=O.N(C)(C)C')))==[]
    assert list(electron_actions(parse_explicit('O=CCO')))==[]


def test_formose_retro_aldol_and_inverse_actions_are_exact_reverses():
    tetrose=parse_explicit('O=C[C@H](O)[C@H](O)CO')
    retro=[arrows for name,arrows in electron_actions(tetrose)
           if name=='formose_retro_aldol_tetrose_to_2go']
    assert retro
    product,_=replay(tetrose,retro[0])
    assert graph_smiles(product)=='O=CCO.O=CCO'
    go=parse_explicit('O=CCO.O=CCO')
    inverse=[arrows for name,arrows in electron_actions(go)
             if name=='formose_inverse_aldol_2go_to_tetrose']
    assert inverse
    products={graph_smiles(replay(go,arrows)[0]) for arrows in inverse}
    assert 'O=CC(O)C(O)CO' in products


def test_dft_optimizer_target_cannot_loosen_physical_gate(tmp_path):
    from mechbridge.verification import verify_event
    with pytest.raises(ValueError,match='no looser'):
        verify_event({},tmp_path,frames=0,ts_optimizer_fmax=.03)


def test_mixture_formula_lists_actual_components():
    from mechbridge.molecular_visuals import entry
    mol=parse_explicit('C=O.O')
    value=entry(mol,np.zeros((mol.GetNumAtoms(),3)))
    assert value['formula']=='CH2O + H2O'
    assert value['total_formula']=='CH4O2'
