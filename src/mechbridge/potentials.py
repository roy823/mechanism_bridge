"""Select a pinned reaction-search potential and return explicit provenance."""
import hashlib,importlib.metadata
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_potential(name,root,device='cpu',compile_model=False,threads=2):
    root=Path(root)
    if name=='aimnet2-rxn':
        from .aimnet_backend import ReactionPotential
        folder=root/'models/aimnet2-rxn';backend=ReactionPotential(folder)
        model_file=folder/'ensemble_0.safetensors'
        details=dict(reference='wB97M-V/def2-TZVPP (pinned HF artifact)',
            elements=[1,6,7,8],official_registry_D3_default=False,
            note='Historical project baseline explicitly suppresses D3 according to its pinned artifact card')
    elif name in ('aimnet2','aimnet2-2025','aimnet2-nse'):
        from .aimnetcentral_backend import AIMNetCentralPotential,FAMILIES
        folder=root/'models/aimnetcentral';backend=AIMNetCentralPotential(name,folder,device=device,
            compile_model=compile_model,threads=threads)
        model_file=folder/FAMILIES[name]['file']
        details=dict(reference={'aimnet2':'wB97M-D3/def2-TZVPP','aimnet2-2025':'B97-3c',
            'aimnet2-nse':'wB97M-D3 open-shell'}[name],elements=sorted(backend.supported_species),
            official_registry_D3_default=backend.predict.external_dftd3 is not None,
            note='Official AIMNetCentral metadata defaults retained')
    else:raise ValueError(f'Unknown potential: {name}')
    packages={}
    for package in ['aimnet','torch','ase','numpy','rdkit','scipy']:
        packages[package]=importlib.metadata.version(package)
    provenance=dict(model=name,device=backend.device if hasattr(backend,'device') else 'cpu',
        compile_model=bool(compile_model),model_file=model_file.relative_to(root).as_posix(),
        weights_sha256=sha256(model_file),packages=packages,**details)
    return backend,provenance
