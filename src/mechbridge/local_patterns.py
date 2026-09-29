"""Compile electron actions into valence- and charge-constrained local SMARTS."""
from rdkit import Chem


def compile_pattern(mol, arrows):
    """Keep arrow atoms and one heavy-atom shell; allow external substitution.

    Only arrow-participating hydrogens are included. Valence and formal charge
    are explicit queries; an omitted spectator H can therefore become an alkyl
    substituent without changing the electron inventory of the active atom.
    """
    active = {i for a in arrows for site in ('source', 'sink') for i in a[site]}
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in active):
        raise ValueError('Aromatic active sites require a separate Kekule action model')
    selected = sorted(active | {n.GetIdx() for i in active
        for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() != 1})
    ids = {old: new for new, old in enumerate(selected)}
    query = Chem.RWMol()
    for i in selected:
        atom = mol.GetAtomWithIdx(i)
        aromatic = 'a' if atom.GetIsAromatic() else 'A'
        query.AddAtom(Chem.AtomFromSmarts(
            f'[#{atom.GetAtomicNum()};{atom.GetFormalCharge():+d};{aromatic};v{atom.GetTotalValence()}]'))
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in ids and j in ids:
            query.AddBond(ids[i], ids[j], bond.GetBondType())
    local_arrows = [dict(a, source=[ids[i] for i in a['source']],
                        sink=[ids[i] for i in a['sink']]) for a in arrows]
    return query.GetMol(), local_arrows


def matches_pattern(mol, query, match):
    """Reject extra bonds inside the matched core; external substituents survive."""
    for i in range(len(match)):
        for j in range(i):
            expected = query.GetBondBetweenAtoms(i, j)
            actual = mol.GetBondBetweenAtoms(match[i], match[j])
            if (expected is None) != (actual is None):
                return False
    return True
