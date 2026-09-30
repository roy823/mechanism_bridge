"""Fixed literature cases with explicit distinctions between illustration and evaluation."""
import json,sys,hashlib,zipfile
from pathlib import Path
import numpy as np
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.symbolic_library import parse_explicit,ArrowLibrary
from mechbridge.encounters import assemble_encounter
from mechbridge.event_graph import graph_smiles,geometry_mol,resonance_equivalent


def main():
    RDLogger.DisableLog('rdApp.*')
    out=ROOT/'data/processed/literature_benchmark_starts.jsonl'
    if out.exists():raise FileExistsError(out)
    audit=json.loads((ROOT/'reports/literature_benchmark_v6/flower_test_selection_audit.json').read_text())
    candidate=audit['candidates'][0]
    cases=[dict(id='flower_epoxide_ammonia',reactant='C1CO1.N',targets=['[NH3+]CC[O-]','NCCO'],
        source='https://arxiv.org/pdf/2502.12979',location='Fig. 1b, PDF page 4',
        evidence='paper_illustration_not_individually_DFT_validated',reference_model_success_known=False),
        dict(id='atlas_glycolaldehyde_hydration',reactant='O=CCO.O',targets=['OCC(O)O'],
        source='https://arxiv.org/html/2606.30778v1',location='S2.3.4, glycolaldehyde gem-diol discussion',
        evidence='paper_reports_gas_phase_DFT_barrier_for_same_reaction',reference_model_success_known=True),
        dict(id='atlas_formaldehyde_dimer',reactant='C=O.C=O',targets=['O=CCO'],
        source='https://arxiv.org/html/2606.30778v1',location='S2.3.4, formose initiation',
        evidence='paper_reports_alternative_TS_for_same_endpoint_pair',reference_model_success_known=True),
        dict(id='flower_test_tropylium_cyanide',reactant=candidate['reactant'],targets=[candidate['product']],
        source='https://doi.org/10.6084/m9.figshare.28359407.v3',location=f"flower_dataset/test.txt line {candidate['line']}",
        evidence='exact_public_test_record_individual_prediction_unavailable',reference_model_success_known=False,
        record_sha256=hashlib.sha256(candidate['record'].encode()).hexdigest())]
    library=ArrowLibrary(ROOT/'data/raw/synepd/polar.json');starts=[]
    for k,case in enumerate(cases):
        mol=parse_explicit(case['reactant']);fragments=Chem.GetMolFrags(mol)
        xyz=np.zeros((mol.GetNumAtoms(),3))
        for ids,fragment in zip(fragments,Chem.GetMolFrags(mol,asMols=True)):
            params=AllChem.ETKDGv3();params.randomSeed=32000+k
            confs=list(AllChem.EmbedMultipleConfs(fragment,numConfs=4,params=params))
            energies=AllChem.UFFOptimizeMoleculeConfs(fragment,maxIters=500,numThreads=1)
            good=[i for i,(status,_) in enumerate(energies) if status==0]
            if not good:raise ValueError('No reactant conformer: '+case['id'])
            best=min(good,key=lambda i:energies[i][1]);xyz[list(ids)]=fragment.GetConformer(confs[best]).GetPositions()
        numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
        proposals=library.propose(mol)
        case['symbolic_proposals']=len(proposals);case['atoms']=len(numbers)
        case_starts=[];case['admission']='accepted'
        for orientation in range(2):
            seed=33000+k*10+orientation;x=assemble_encounter(numbers,xyz,fragments,seed)
            try:observed=geometry_mol(numbers,x,0)
            except ValueError as exc:
                case['admission']='blocked_geometry_perception'
                case['admission_reason']=str(exc)
                case_starts=[]
                break
            if not resonance_equivalent(mol,observed):raise ValueError('Input perception mismatch')
            case_starts.append(dict(id=case['id']+f'_o{orientation}',atomic_numbers=numbers,positions_A=x.tolist(),
                charge=0,multiplicity=1,provenance=dict(system=case['id'],input_smiles=graph_smiles(observed),
                source=case['source'],source_evidence=case['evidence'],random_seed=seed,
                reference_TS_or_product_geometry_used=False)))
        starts.extend(case_starts)
    out.write_text(''.join(json.dumps(r)+'\n' for r in starts),encoding='utf-8')
    with zipfile.ZipFile(ROOT/'data/raw/flower/data.zip') as z:beam=z.read('data/flower_dataset/beam.txt').decode().strip()
    report=dict(cases=cases,FlowER_original_flagship=dict(record=beam,executed=False,
        reason='Contains Br, F, Cl and more than two components; outside current MLIP domain; no atom deletion'),
        reference_models_executed=False,model='AIMNet2-rxn ensemble_0',reference_targets_used_only_for_post_search_scoring=True)
    (ROOT/'reports/literature_benchmark_v6/cases.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps([{k:c[k] for k in ['id','atoms','symbolic_proposals','evidence']} for c in cases],indent=2))


if __name__=='__main__':main()
