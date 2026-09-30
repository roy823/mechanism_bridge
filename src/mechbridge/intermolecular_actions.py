"""Small, explicit electron-action grammar for neutral two-fragment encounters.

These are analyst-specified hypotheses, not published arrow labels or learned
predictions. Only the input graph is read; no reference product or TS is used.
"""
from rdkit import Chem


def electron_actions(mol):
    fragments = Chem.GetMolFrags(mol)
    if len(fragments) != 2:
        return
    component = {i: k for k, f in enumerate(fragments) for i in f}

    def matches(smarts):
        return mol.GetSubstructMatches(Chem.MolFromSmarts(smarts), uniquify=False)

    def action(name, flows):
        return name, [dict(source=list(s), sink=list(t), electrons=2,
                           kind='reactant_graph_grammar') for s, t in flows]

    carbonyls = matches('[C;+0]=[O;+0]')
    donors = matches('[O,N;+0]-[H]')
    enols = matches('[C;+0]=[C;+0]-[O;+0]-[H]')
    # Addition with an explicit donor proton; no invented proton or catalyst.
    for c, o in carbonyls:
        for nu, h in donors:
            if component[c] != component[nu]:
                yield action('carbonyl_addition_with_proton_transfer',
                    [((nu,), (nu,c)), ((c,o), (o,h)), ((nu,h), (nu,))])
        for ca, cb, oe, h in enols:
            if component[c] != component[ca]:
                yield action('neutral_enol_aldol_addition',
                    [((ca,cb), (ca,c)), ((c,o), (o,h)), ((oe,h), (oe,cb))])
        # A possible concerted carbonyl C-H coupling, not a low-barrier claim.
        for donor, od in carbonyls:
            if component[c] == component[donor]:
                continue
            for atom in mol.GetAtomWithIdx(donor).GetNeighbors():
                if atom.GetAtomicNum() == 1:
                    h = atom.GetIdx()
                    yield action('carbonyl_C_H_coupled_addition',
                        [((donor,h), (donor,c)), ((c,o), (o,h))])
    # Covers localized nitrone/nitrile-ylide/nitrile-imine representations.
    for a, b, c in matches('[C,N;+0]=,#[N,O;+1]-,=[C,N,O;-1]'):
        for d, e in matches('[C;+0]=,#[C;+0]'):
            if component[a] != component[d]:
                yield action('dipolar_3_plus_2_cycloaddition',
                    [((c,), (c,d)), ((d,e), (e,a)), ((a,b), (b,))])


def crosses_components(mol, edits):
    component = {i: k for k, f in enumerate(Chem.GetMolFrags(mol)) for i in f}
    return any(e['before'] == 0 and component[e['atoms'][0]] != component[e['atoms'][1]]
               for e in edits)


def diverse_proposals(proposals, limit):
    """Round-robin product hypotheses before symmetry mappings; inter first."""
    groups = {}
    for p in proposals:
        groups.setdefault((not p['intermolecular'], p['predicted_graph']), []).append(p)
    result = []
    keys = sorted(groups)
    while keys and len(result) < limit:
        for key in keys:
            result.append(groups[key].pop(0))
            if len(result) == limit:
                return result
        keys = [key for key in keys if groups[key]]
    return result
