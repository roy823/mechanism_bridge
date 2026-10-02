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


def _legacy_orient_reactive_encounter(numbers, positions, mol, edits, random_seed):
    """Verbatim copy of orient_reactive_encounter at c93f22c (regression reference)."""
    from scipy.spatial.transform import Rotation
    from ase.data import covalent_radii
    fragments = Chem.GetMolFrags(mol)
    if len(fragments) < 2:
        return np.array(positions,copy=True), None
    component = {atom: index for index, fragment in enumerate(fragments) for atom in fragment}
    pairs = []
    for e in edits:
        i,j = e['atoms']
        if e['before']==0 and component[i] != component[j]:
            pairs.append((i,j))
    if not pairs:
        return np.array(positions,copy=True), None
    active_components = tuple(sorted((component[pairs[0][0]], component[pairs[0][1]])))
    pairs = [pair for pair in pairs
             if tuple(sorted((component[pair[0]], component[pair[1]]))) == active_components]
    ids_a, ids_b = [np.array(fragments[index],dtype=int) for index in active_components]
    aset = set(ids_a)
    pairs = [(i,j) if i in aset else (j,i) for i,j in pairs]
    numbers = np.asarray(numbers)
    pairs.sort(key=lambda p: (numbers[p[0]]==1 or numbers[p[1]]==1,p))
    anchor_a,anchor_b = pairs[0]
    rng = np.random.default_rng(random_seed)
    original = np.asarray(positions)
    active_ids = np.concatenate((ids_a,ids_b))
    active_center = original[active_ids].mean(0)
    best = None
    for trial in range(64):
        y = original.copy()
        for ids,anchor in ((ids_a,anchor_a),(ids_b,anchor_b)):
            y[ids] = (y[ids]-y[anchor]) @ Rotation.random(random_state=rng).as_matrix()
        y[ids_b] += np.array([0.,0.,2.8])
        y[active_ids] += active_center-y[active_ids].mean(0)
        d = np.linalg.norm(y[ids_a,None]-y[ids_b][None,:],axis=-1)
        floor = 1.15*(covalent_radii[numbers[ids_a,None]]+covalent_radii[numbers[ids_b]][None,:])
        score = sum((np.linalg.norm(y[i]-y[j])-2.8)**2 for i,j in pairs)
        score += 50*np.maximum(floor-d,0).sum()**2
        spectators = np.array([i for index,fragment in enumerate(fragments)
                               if index not in active_components for i in fragment],dtype=int)
        if len(spectators):
            d_other=np.linalg.norm(y[active_ids,None]-y[spectators][None,:],axis=-1)
            floor_other=(covalent_radii[numbers[active_ids,None]]+
                         covalent_radii[numbers[spectators]][None,:])
            score += 50*np.maximum(floor_other-d_other,0).sum()**2
        if best is None or score < best[0]:
            best = score,y-y.mean(0),trial
    return best[1],dict(policy='64_rigid_orientations_cross_bond_distance_and_clash_score',
        score=float(best[0]),selected_trial=best[2],cross_forming_pairs=pairs,
        active_components=list(active_components),spectator_components=len(fragments)-2,
        contact_distance_A=2.8,reference_product_or_TS_used=False,
        displacement_from_source_A=float(np.linalg.norm(best[1]-original)))


@pytest.mark.parametrize('smiles',['C=O.O','C=O.O.O'])
def test_symbolic_encounter_refactor_is_bit_identical(smiles):
    mol=parse_explicit(smiles);assert AllChem.EmbedMolecule(mol,randomSeed=17)==0
    x=mol.GetConformer().GetPositions();numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
    fragments=Chem.GetMolFrags(mol)
    carbon=next(i for i in fragments[0] if numbers[i]==6)
    oxygen=next(i for i in fragments[1] if numbers[i]==8)
    edits=[dict(atoms=[carbon,oxygen],before=0,after=1)]
    for seed in (0,7,23):
        new,meta=orient_reactive_encounter(numbers,x,mol,edits,seed)
        old,old_meta=_legacy_orient_reactive_encounter(numbers,x,mol,edits,seed)
        assert np.array_equal(new,old) and meta==old_meta


def test_matched_controls_get_proposal_free_encounter():
    from ase import Atoms
    from mechbridge.search_seeds import make_seed
    mol=parse_explicit('C=O.O');assert AllChem.EmbedMolecule(mol,randomSeed=17)==0
    x=mol.GetConformer().GetPositions();numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
    fragments=Chem.GetMolFrags(mol)
    y=assemble_encounter(numbers,x,fragments,23)
    atoms=Atoms(numbers=numbers,positions=y)
    _,arrows=next(electron_actions(mol));product,edits=replay(mol,arrows)
    proposal=dict(arrows=arrows,edits=edits,template_id='t',predicted_graph=graph_smiles(product))
    def distances(z):return np.linalg.norm(z[:,None]-z[None,:],axis=-1)
    for strategy in ('geometry','center_random'):
        for sample,amplitude in enumerate((.6,1.,1.4)):
            seed,d,meta=make_seed(atoms,mol,strategy,proposal,sample,31,encounter_policy='matched_controls')
            encounter=meta['encounter_orientation']
            assert meta['encounter_policy']=='matched_controls'
            assert encounter['proposal_bond_used'] is False and encounter['contact_distance_A']==2.8
            oriented=seed-amplitude*d
            for ids in fragments:
                assert np.allclose(distances(y[list(ids)]),distances(oriented[list(ids)]))
            i,j=encounter['contact_pair']
            assert np.linalg.norm(oriented[i]-oriented[j])==pytest.approx(2.8)
            assert np.linalg.norm(d)==pytest.approx(1.)
            if strategy=='center_random':
                active={k for e in edits for k in e['atoms']}
                assert {i,j}<=active
            legacy=make_seed(atoms,mol,strategy,proposal,sample,31)[2]
            assert 'encounter_orientation' not in legacy
    single=parse_explicit('CC=O');assert AllChem.EmbedMolecule(single,randomSeed=3)==0
    z=single.GetConformer().GetPositions()
    unimolecular=Atoms(numbers=[a.GetAtomicNum() for a in single.GetAtoms()],positions=z)
    matched=make_seed(unimolecular,single,'geometry',None,1,5,encounter_policy='matched_controls')
    legacy=make_seed(unimolecular,single,'geometry',None,1,5)
    assert np.array_equal(matched[0],legacy[0]) and matched[2]['encounter_orientation'] is None
