"""Graph bookkeeping, not mechanistic inference. No templates or reagent completion."""
from collections import Counter
import hashlib
import json
from rdkit import Chem

def mol(smiles):
    p = Chem.SmilesParserParams()
    p.removeHs = False
    m = Chem.MolFromSmiles(smiles, p)
    if m is None:
        raise ValueError("Invalid SMILES")
    return m

def canonical(smiles):
    m = mol(smiles)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    # Ordinary explicit H and implicit H share an identity; isotope/stereo H retained by RDKit.
    m = Chem.RemoveHs(m)
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)

def identity(smiles):
    m = mol(smiles)
    h = Chem.AddHs(m)
    return {"canonical_smiles": canonical(smiles),
            "composition": dict(sorted(Counter(str(a.GetAtomicNum()) for a in h.GetAtoms()).items())),
            "charge": sum(a.GetFormalCharge() for a in m.GetAtoms()),
            "radical_electrons": sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()),
            "all_h_explicit": all(a.GetNumImplicitHs() + a.GetNumExplicitHs() == 0 for a in m.GetAtoms())}

def split_reaction(rxn):
    fields = rxn.split(">")
    if len(fields) != 3 or not fields[0] or not fields[2]:
        raise ValueError("Expected nonempty reactant>agents>product SMILES")
    return fields[0], fields[2], fields[1]

def graph_pair_key(reactant, product):
    pair = sorted([canonical(reactant), canonical(product)])
    return hashlib.sha256(json.dumps(pair).encode()).hexdigest()

def audit_reaction(reactant, product):
    r, p = identity(reactant), identity(product)
    def maps(s):
        atoms = list(mol(s).GetAtoms())
        ids = [a.GetAtomMapNum() for a in atoms]
        return ids, {a.GetAtomMapNum(): (a.GetAtomicNum(), a.GetIsotope()) for a in atoms}
    rm, ra = maps(reactant); pm, pa = maps(product)
    complete = bool(rm and pm) and min(rm + pm) > 0 and len(set(rm)) == len(rm) and len(set(pm)) == len(pm)
    return {"reactant": r, "product": p, "composition_conserved": r["composition"] == p["composition"],
            "charge_conserved": r["charge"] == p["charge"], "complete_unique_atom_maps": complete,
            "mapped_atom_identity_conserved": complete and ra == pa,
            "graph_pair_key": graph_pair_key(reactant, product)}

def net_changes(reactant, product):
    """Return net graph edits only. These do NOT define unique curved arrows."""
    a = audit_reaction(reactant, product)
    if not a["mapped_atom_identity_conserved"]:
        raise ValueError("Net changes require complete conserved unique atom maps")
    def bonds(s):
        m = mol(s)
        return {tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()))):
                b.GetBondTypeAsDouble() for b in m.GetBonds()}
    rb, pb = bonds(reactant), bonds(product)
    return [{"atom_maps": list(k), "before": rb.get(k, 0), "after": pb.get(k, 0)}
            for k in sorted(rb.keys() | pb.keys()) if rb.get(k, 0) != pb.get(k, 0)]
