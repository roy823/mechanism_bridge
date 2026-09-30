"""Classify actual atom-indexed endpoints, independently of seed intent."""
from rdkit import Chem
from .event_graph import bond_orders,graph_smiles,resonance_equivalent


def matches_reference_pair(endpoints,reactant,products):
    """Both sides of one observed event must match; seeing a product alone is insufficient."""
    return len(endpoints)==2 and any(resonance_equivalent(endpoints[i],reactant) and
        any(resonance_equivalent(endpoints[1-i],p) for p in products) for i in (0,1))


def classify_event(left,right):
    fragments=[Chem.GetMolFrags(m) for m in (left,right)]
    labels=[{i:k for k,f in enumerate(fs) for i in f} for fs in fragments]
    bonds=[bond_orders(m) for m in (left,right)]
    cross=[]
    for side in (0,1):
        for (i,j),order in bonds[1-side].items():
            if (i,j) not in bonds[side] and labels[side][i]!=labels[side][j]:
                cross.append(dict(from_endpoint=side,atoms=[i,j],
                    heavy_atom_bond=left.GetAtomWithIdx(i).GetAtomicNum()>1 and
                                   left.GetAtomWithIdx(j).GetAtomicNum()>1))
    equivalent=resonance_equivalent(left,right)
    if equivalent:kind='same_species_geometry_or_resonance'
    elif any(p['heavy_atom_bond'] for p in cross):kind='intermolecular_heavy_atom_bond'
    elif cross:kind='intermolecular_hydrogen_transfer'
    elif max(map(len,fragments))>1:kind='intramolecular_with_spectator'
    else:kind='intramolecular'
    return dict(classification=kind,cross_fragment_bonds=cross,
        fragment_counts=[len(f) for f in fragments],resonance_equivalent=equivalent,
        endpoint_components=[[graph_smiles(m) for m in Chem.GetMolFrags(mol,asMols=True)]
                             for mol in (left,right)])
