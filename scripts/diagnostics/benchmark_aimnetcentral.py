"""Benchmark broad AIMNetCentral models on current search operations and reaction evidence."""
import argparse,hashlib,json,platform,statistics,sys,time
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnetcentral_backend import AIMNetCentralPotential,FAMILIES
from mechbridge.physics import finite_hessian,vibrational_analysis
from mechbridge.symbolic_library import parse_explicit


def synchronize(device):
    if str(device).startswith('cuda'):
        import torch;torch.cuda.synchronize()


def fragments_geometry(smiles,seed):
    mol=parse_explicit(smiles);xyz=np.zeros((mol.GetNumAtoms(),3));offset=0.
    for ids,fragment in zip(Chem.GetMolFrags(mol),Chem.GetMolFrags(mol,asMols=True)):
        params=AllChem.ETKDGv3();params.randomSeed=seed
        if AllChem.EmbedMolecule(fragment,params)!=0:raise ValueError('Embedding failed: '+smiles)
        try:AllChem.UFFOptimizeMolecule(fragment,maxIters=500)
        except ValueError:pass
        part=fragment.GetConformer().GetPositions();part-=part.mean(0);part[:,0]+=offset
        xyz[list(ids)]=part;offset+=max(np.ptp(part[:,0])+4.,5.)
    return np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),xyz


def timed(call,repeats,device):
    values=[]
    for _ in range(repeats):
        start=time.perf_counter();call();synchronize(device);values.append(time.perf_counter()-start)
    return dict(median_seconds=statistics.median(values),mean_seconds=statistics.mean(values),
                min_seconds=min(values),repeats=repeats)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--device',choices=['cpu','cuda'],default='cpu');p.add_argument('--compile',action='store_true')
    p.add_argument('--families',nargs='+',choices=list(FAMILIES),default=list(FAMILIES));a=p.parse_args()
    if a.outdir.exists():raise FileExistsError(a.outdir)
    a.outdir.mkdir(parents=True);RDLogger.DisableLog('rdApp.*')
    systems={
        'chno_11':fragments_geometry('O=CCO.O',41),
        'flowER_halogenated_organic_pair':fragments_geometry('CC(=O)CC(=O)C(F)(F)F.NNc1cccc(Br)c1',43),
        'main_group_P_S_Cl':fragments_geometry('CP(=O)(O)S.CCl',47)}
    rng=np.random.default_rng(1701);models=[];families={}
    for family in a.families:
        started=time.perf_counter()
        model=AIMNetCentralPotential(family,ROOT/'models/aimnetcentral',device=a.device,
            compile_model=a.compile,threads=2)
        synchronize(a.device);loaded=time.perf_counter()-started
        row=dict(family=family,description=model.describe(),load_seconds=loaded,systems={})
        if a.device=='cuda':
            import torch;torch.cuda.reset_peak_memory_stats()
        for name,(numbers,xyz) in systems.items():
            try:model.validate_system(numbers,0,1)
            except ValueError as exc:
                row['systems'][name]=dict(supported=False,reason=str(exc));continue
            geometries=xyz[None]+rng.normal(0,.01,(32,*xyz.shape))
            model.evaluate_many(numbers,geometries[:1]);model.evaluate_many(numbers,geometries)
            single=timed(lambda:model.evaluate_many(numbers,geometries[:1]),20,a.device)
            batch=timed(lambda:model.evaluate_many(numbers,geometries),7,a.device)
            row['systems'][name]=dict(supported=True,atoms=len(numbers),elements=sorted(map(int,set(numbers))),
                single=single,batch32=batch,
                batch32_geometries_per_second=32/batch['median_seconds'])
        # This matches the current finite-difference workflow rather than a
        # model-native Hessian that cannot be combined with compilation.
        numbers,xyz=systems['chno_11'];atoms=Atoms(numbers=numbers,positions=xyz)
        started=time.perf_counter();h=finite_hessian(atoms,step=.005,
            batch_forces=lambda n,x:model.evaluate_many(n,x)['forces'],batch_size=32)
        synchronize(a.device);hessian_seconds=time.perf_counter()-started
        vib=vibrational_analysis(xyz,atoms.get_masses(),h)
        row['finite_difference_hessian']=dict(seconds=hessian_seconds,
            imaginary_count=vib['imaginary_count'],geometries=6*len(numbers),batch_size=32)
        import torch
        def native_hessian():
            with torch.compiler.set_stance('force_eager'):
                return model.predict(model._input(numbers,xyz),forces=True,hessian=True)
        native=native_hessian();synchronize(a.device)
        native_time=timed(native_hessian,3,a.device)
        native_h=np.asarray(native['hessian'].detach().cpu()).reshape(h.shape)
        native_vib=vibrational_analysis(xyz,atoms.get_masses(),native_h)
        row['native_hessian']=dict(**native_time,imaginary_count=native_vib['imaginary_count'],
            max_abs_difference_from_finite_eV_A2=float(np.max(abs(native_h-h))))
        if a.device=='cuda':
            import torch;row['peak_cuda_memory_MiB']=torch.cuda.max_memory_allocated()/2**20
        families[family]=model;models.append(row)
    ts_cases=[]
    case_specs=[
        ('acetone_tautomerization',ROOT/'reports/local_transfer_v3/qc/acetone'),
        ('formaldehyde_hydration',ROOT/'reports/bimolecular_v5/qc/formaldehyde_water'),
        ('formaldehyde_dimerization',ROOT/'reports/bimolecular_v5/qc/formaldehyde_dimer_refined')]
    for name,folder in case_specs:
        q=json.loads((folder/'verification.json').read_text());ts=read(folder/'ts.xyz')
        ends=[read(folder/f'minimum_{d}.xyz') for d in ('forward','reverse')]
        case=dict(case=name,DFT_method=q['method']+'/'+q['basis'],models=[])
        for family,model in families.items():
            if not set(ts.numbers)<=model.supported_species:continue
            positions=np.stack([ts.positions,ends[0].positions,ends[1].positions])
            values=model.evaluate_many(ts.numbers,positions)
            h=finite_hessian(ts,step=.005,
                batch_forces=lambda n,x:model.evaluate_many(n,x)['forces'],batch_size=32)
            vib=vibrational_analysis(ts.positions,ts.get_masses(),h)
            case['models'].append(dict(family=family,
                force_max_eV_A=float(np.linalg.norm(values['forces'][0],axis=1).max()),
                imaginary_count=vib['imaginary_count'],lowest_frequencies_cm_1=vib['frequencies_cm-1'][:3],
                barriers_from_DFT_endpoints_eV=[float(values['energy'][0]-v) for v in values['energy'][1:]]))
        ts_cases.append(case)
    # The pinned HF rxn artifact and current official registry artifact should
    # agree numerically; a failure means this is not a clean model-family comparison.
    consistency=None
    if a.device=='cpu' and 'aimnet2-rxn' in families:
        from mechbridge.aimnet_backend import ReactionPotential
        old=ReactionPotential(ROOT/'models/aimnet2-rxn');numbers,xyz=systems['chno_11']
        old_result=old.evaluate_many(numbers,xyz[None]);official=families['aimnet2-rxn'].evaluate_many(numbers,xyz[None])
        no_d3=AIMNetCentralPotential('aimnet2-rxn',ROOT/'models/aimnetcentral',device='cpu',threads=2,
                                     needs_dispersion=False)
        matched=no_d3.evaluate_many(numbers,xyz[None])
        consistency=dict(
            official_registry_default_vs_pinned=dict(
                energy_abs_difference_eV=float(abs(old_result['energy'][0]-official['energy'][0])),
                force_max_abs_difference_eV_A=float(np.max(abs(old_result['forces']-official['forces']))),
                cause='Current registry family default adds external D3; pinned historical artifact explicitly suppresses it'),
            registry_weights_without_D3_vs_pinned=dict(
                energy_abs_difference_eV=float(abs(old_result['energy'][0]-matched['energy'][0])),
                force_max_abs_difference_eV_A=float(np.max(abs(old_result['forces']-matched['forces'])))))
    files=sorted((ROOT/'models/aimnetcentral').glob('*.pt'))
    result=dict(device=a.device,compile_model=a.compile,python=sys.version,platform=platform.platform(),
        package_versions={},models=models,reaction_TS_diagnostics=ts_cases,
        current_rxn_artifact_consistency=consistency,
        model_files=[dict(file=f.relative_to(ROOT).as_posix(),bytes=f.stat().st_size,
            sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files],
        timing_scope='Warm energy+force inference; current finite-difference Hessian; two CPU threads when CPU',
        evidence_limits='Timing and local PES diagnostics, not a chemical accuracy benchmark or TS success rate')
    import importlib.metadata
    for pkg in ['aimnet','torch','ase','numpy','rdkit','nvalchemi-toolkit-ops','warp-lang']:
        try:result['package_versions'][pkg]=importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:result['package_versions'][pkg]=None
    temporary=a.outdir/'results.json.tmp'
    temporary.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    temporary.replace(a.outdir/'results.json')
    print(json.dumps(dict(device=a.device,models=[dict(family=r['family'],load=r['load_seconds'],
        single=r['systems']['chno_11']['single']['median_seconds'],
        batch32=r['systems']['chno_11']['batch32']['median_seconds'],
        hessian=r['finite_difference_hessian']['seconds']) for r in models],consistency=consistency),indent=2))


if __name__=='__main__':main()
