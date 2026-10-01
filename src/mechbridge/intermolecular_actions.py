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
    dienes = matches('[C;+0]=[C;+0]-[C;+0]=[C;+0]')
    alkenes = matches('[C;+0]=[C;+0]')
    peracids = matches('[C;+0](=[O;+0])-[O;+0]-[O;+0]-[H]')
    # Both the symbolic zwitterion step and a concerted proton-transfer route
    # are hypotheses. The PES decides whether the zwitterion is a minimum.
    for c, other, o in matches('[C;+0]1[C;+0][O;+0]1'):
        for (nu,) in matches('[N;+0;v3]'):
            if component[c] == component[nu]:continue
            yield action('epoxide_amine_opening_stepwise',
                [((nu,), (nu,c)), ((c,o), (o,))])
            for atom in mol.GetAtomWithIdx(nu).GetNeighbors():
                if atom.GetAtomicNum()!=1:continue
                h=atom.GetIdx()
                yield action('epoxide_amine_opening_with_proton_transfer',
                    [((nu,), (nu,c)), ((c,o), (o,h)), ((nu,h), (nu,))])
    # Addition with an explicit donor proton; no invented proton or catalyst.
    for c, o in carbonyls:
        for (nu,) in matches('[N;+0;v3]'):
            if component[c] != component[nu]:
                yield action('carbonyl_amine_addition_stepwise',
                    [((nu,), (nu,c)), ((c,o), (o,))])
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
    # FlowER Fig. 4 / Fig. S16 reaction classes. These are minimal, explicit
    # reactant-only hypotheses; the downstream saddle search decides whether a
    # matching physical event exists for the supplied encounter geometry.
    for a, b, c, d in dienes:
        for e, f in alkenes:
            if component[a] != component[e]:
                yield action('diels_alder_4_plus_2_cycloaddition',
                    [((a,b), (a,e)), ((e,f), (f,d)), ((c,d), (b,c))])
    for e, f in alkenes:
        for acyl, carbonyl_o, proximal_o, terminal_o, h in peracids:
            if component[e] != component[acyl]:
                yield action('prilezhaev_epoxidation',
                    [((e,f), (e,terminal_o)),
                     ((proximal_o,terminal_o), (acyl,proximal_o)),
                     ((acyl,carbonyl_o), (carbonyl_o,h)),
                     ((terminal_o,h), (f,terminal_o))])


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
