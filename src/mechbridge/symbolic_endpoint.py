"""Construct a three-dimensional endpoint from an atom-indexed arrow proposal.

The endpoint is a hypothesis generated from the symbolic graph.  It never uses a
reference product conformer or transition state.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from .event_graph import graph_smiles
from .symbolic_library import replay


def product_geometry(numbers, reactant_positions, reactant_mol, proposal, random_seed):
    """Embed the replayed product while fixing its unchanged heavy-atom core."""
    product, edits = replay(reactant_mol, proposal['arrows'])
    if graph_smiles(product) != proposal['predicted_graph']:
        raise ValueError('Arrow replay does not match the proposed product graph')
    product = Chem.Mol(product)
    product.RemoveAllConformers()
    active = {i for edit in edits for i in edit['atoms']}
    core = [i for i, z in enumerate(numbers) if z > 1 and i not in active]
    # A coordinate map keeps the largest unchanged part in the observed frame.
    # With fewer than two anchors ETKDG is more stable without a partial map.
    coord_map = ({i: Point3D(*map(float, reactant_positions[i])) for i in core}
                 if len(core) >= 2 else {})
    conformer = AllChem.EmbedMolecule(product, randomSeed=int(random_seed & 0x7fffffff),
        useRandomCoords=True, coordMap=coord_map, enforceChirality=True)
    if conformer < 0:
        raise ValueError('RDKit could not embed the symbolically replayed product')
    positions = np.asarray(product.GetConformer(conformer).GetPositions(), dtype=float)
    if not np.isfinite(positions).all():
        raise ValueError('Embedded product contains nonfinite coordinates')
    return positions, dict(predicted_graph=proposal['predicted_graph'], active_atoms=sorted(active),
        fixed_core_atoms=core, generator='RDKit_ETKDG_from_arrow_replay',
        reference_product_geometry_used=False)
