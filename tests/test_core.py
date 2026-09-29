import json
import numpy as np
import pytest
from mechbridge.chemistry import canonical, audit_reaction, graph_pair_key, net_changes
from mechbridge.adapters import flow_lines, arrow_sites, rgd1, transition1x
from mechbridge.pairing import pair_candidates, split_records
from mechbridge.physics import finite_hessian, vibrational_analysis
from mechbridge.io import event

def test_maps_hydrogens_and_stereo():
    assert canonical('[H:8][O:9][H:10]') == canonical('O')
    assert canonical('F[C@H](Cl)Br') != canonical('F[C@@H](Cl)Br')
    assert canonical('[2H]O') != canonical('O')
    assert not audit_reaction('[CH3:1][OH:2]', '[CH3:1][OH:1]')['complete_unique_atom_maps']
    assert not audit_reaction('[CH3:1][OH:2]', '[CH2:1]=[O:2]')['composition_conserved']
    assert not audit_reaction('[Na+:1].[Cl-:2]', '[Na:1].[Cl-:2]')['charge_conserved']

def test_net_changes_are_not_arrows():
    r, p = '[CH2:1]=[CH2:2]', '[CH3:1][CH3:2]'
    assert net_changes(r,p) == [{'atom_maps':[1,2],'before':2.0,'after':1.0}]
    assert graph_pair_key(r,p) == graph_pair_key(p,r)

def test_arrow_flat_sequence_preserved():
    assert len(arrow_sites([(1,2),([2,3],3)])) == 2
    with pytest.raises(ValueError): arrow_sites([(1,2,3)])

def test_no_false_pairs_for_overall_or_unbalanced():
    symbolic = list(flow_lines(['[CH3:1][OH:2]>>[CH3:1][OH:2]|1'], 'test', {}))
    physical = event('physical','1','physical_event_candidate',{})
    physical['symbolic'] = dict(symbolic[0]['symbolic'])
    hits = list(pair_candidates(symbolic,[physical]))
    assert len(hits) == 1 and not hits[0]['verified_pair']
    symbolic[0]['kind'] = 'symbolic_overall_reaction'
    assert list(pair_candidates(symbolic,[physical])) == []

def test_sequence_reverse_and_conformer_splits_stay_together():
    records = list(flow_lines(['CCO>>COC|seq1','COC>>CCO|seq2','COC>>CC=O|seq2'], 'test', {}))
    split = list(split_records(records))
    assert len(set(r['component_id'] for r in split)) == 1
    assert len(set(r['split'] for r in split)) == 1

@pytest.mark.parametrize('xyz,rigid,nmodes', [([[0,0,0],[0,0,1]],5,1), ([[0,0,0],[1,0,0],[0,1,0]],6,3)])
def test_rigid_modes_and_single_negative(xyz,rigid,nmodes):
    x=np.array(xyz,dtype=float); masses=np.ones(len(x))
    baseline=vibrational_analysis(x,masses,np.eye(x.size))
    v=baseline['modes'][0].ravel()
    h=np.eye(x.size)-3*np.outer(v,v)
    result=vibrational_analysis(x,masses,h)
    assert result['rigid_mode_count']==rigid
    assert len(result['frequencies_cm-1'])==nmodes
    assert result['imaginary_count']==1
    shifted=vibrational_analysis(x+123,masses,h)
    assert np.allclose(result['frequencies_cm-1'],shifted['frequencies_cm-1'])

def test_finite_hessian_sign_and_restore():
    from ase import Atoms
    from ase.calculators.calculator import Calculator,all_changes
    class Quadratic(Calculator):
        implemented_properties=['energy','forces']
        def calculate(self,atoms=None,properties=None,system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            x=atoms.positions
            self.results={'energy':float((x*x).sum()),'forces':-2*x}
    atoms=Atoms('H2',positions=[[0,0,0],[0,0,1]])
    atoms.calc=Quadratic(); initial=atoms.positions.copy()
    assert np.allclose(finite_hessian(atoms),2*np.eye(6))
    assert np.array_equal(atoms.positions,initial)

def test_rgd1_adapter_does_not_certify_endpoints(tmp_path):
    import h5py
    path=tmp_path/'rgd.h5'
    with h5py.File(path,'w') as f:
        g=f.create_group('fixture_synthetic')
        g['elements']=[1,1]
        for key in ['RG','PG','TSG']:g[key]=[[0.,0,0],[0,0,.74]]
        g['Rsmiles']=b'[H][H]';g['Psmiles']=b'[H][H]'
        for key in ['R_E','P_E','TS_E']:g[key]=-1.
    record=next(rgd1(path))
    assert record['physical']['endpoint_geometry_status']=='unoptimized_in_this_file'
    assert record['validation']['connectivity']=='not_checked'

def test_transition1x_adapter_never_calls_samples_irc(tmp_path):
    import h5py
    path=tmp_path/'t1x.h5'
    with h5py.File(path,'w') as f:
        g=f.create_group('data/H2/fixture_synthetic');g['atomic_numbers']=[1,1]
        for key in ['reactant','transition_state','product']:
            sg=g.create_group(key);sg['positions']=[[[0.,0,0],[0,0,.74]]]
            sg['wB97x_6-31G(d).energy']=[-30.]
    record=next(transition1x(path))
    assert record['physical']['irc'] is None
    assert record['system']['multiplicity'] is None

def test_rgd_csv_reversed_direction_is_aligned_without_changing_energies(tmp_path):
    import h5py,csv
    path=tmp_path/'rgd.h5';csvpath=tmp_path/'mapped.csv'
    with h5py.File(path,'w') as f:
        g=f.create_group('synthetic_1');g['elements']=[6,6,8]
        for key in ['RG','PG','TSG']:g[key]=np.zeros((3,3))
        g['Rsmiles']=b'COC';g['Psmiles']=b'CCO'
        g['R_E']=-2.;g['P_E']=-3.;g['TS_E']=-1.
    with csvpath.open('w') as f:
        w=csv.DictWriter(f,fieldnames=['reaction','reactant','product']);w.writeheader()
        w.writerow({'reaction':'synthetic_1','reactant':'[CH3:1][CH2:2][OH:3]',
                    'product':'[CH3:1][O:3][CH3:2]'})
    r=next(rgd1(path,mapped_csv=csvpath))
    assert r['provenance']['mapped_csv_direction']=='reverse'
    assert canonical(r['symbolic']['reactant_smiles'])=='COC'
    assert r['physical']['source_energies_hartree']['R_E']==-2.
