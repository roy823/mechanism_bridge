"""Small, explicit electron-action grammar for closed-shell CHNOP encounters.

These are analyst-specified hypotheses, not published arrow labels or learned
predictions. Only the input graph is read; no reference product or TS is used.
"""
from rdkit import Chem


def electron_actions(mol):
    fragments = Chem.GetMolFrags(mol)
    if not 1 <= len(fragments) <= 3:
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
    aldotetroses = matches('[C;H1;+0](=[O;+0])-[C;H1;+0](-[O;+0]-[H])-[C;H1;+0](-[O;+0]-[H])-[C;H2;+0]-[O;+0]-[H]')
    glycolaldehydes = matches('[C;H1;+0](=[O;+0])-[C;H2;+0]-[O;+0]-[H]')
    phosphate_esters = matches('[C;+0]-[O;+0]-[P;+0](=[O;+0])(-[O;+0]-[H])-[O;+0]-[H]')
    phosphoric_acids = matches('[P;+0](=[O;+0])(-[O;+0]-[H])(-[O;+0]-[H])-[O;+0]-[H]')
    waters = matches('[O;+0](-[H])-[H]')
    for carboxyl_c, carbonyl_o, carboxylate_o in matches('[C;+0](=[O;+0])-[O;-1]'):
        for phosphorus, phosphoryl_o, acid_o1, acid_h1, acid_o2, acid_h2, acid_o3, acid_h3 in phosphoric_acids:
            if component[carboxylate_o] == component[phosphorus]:continue
            for acid_o, proton in ((acid_o1,acid_h1),(acid_o2,acid_h2),(acid_o3,acid_h3)):
                yield action('phosphoric_acid_to_carboxylate_proton_transfer',
                    [((acid_o,proton), (acid_o,)),
                     ((carboxylate_o,), (carboxylate_o,proton))])
    # Enol pyruvate -> keto pyruvate after PEP phosphate hydrolysis.
    for terminal_c, enol_c, enol_o, enol_h, carboxyl_c, carboxyl_o, acid_o, acid_h in matches(
            '[C;H2;+0]=[C;+0](-[O;+0]-[H])-[C;+0](=[O;+0])-[O;+0]-[H]'):
        yield action('pyruvate_enol_to_keto',
            [((enol_o,enol_h), (enol_o,enol_c)),
             ((terminal_c,enol_c), (terminal_c,enol_h))])
    for methyl_c, keto_c, keto_o, carboxyl_c, carboxyl_o, acid_o, acid_h in matches(
            '[C;H3;+0]-[C;+0](=[O;+0])-[C;+0](=[O;+0])-[O;+0]-[H]'):
        methyl_hydrogens=[a.GetIdx() for a in mol.GetAtomWithIdx(methyl_c).GetNeighbors()
                          if a.GetAtomicNum()==1]
        for proton in methyl_hydrogens:
            yield action('pyruvate_keto_to_enol',
                [((methyl_c,proton), (methyl_c,keto_c)),
                 ((keto_c,keto_o), (keto_o,proton))])
    # Carbonyl/enediol tautomerization. These explicit proton-coupled arrows
    # cover the shared G6P/F6P and GAP/DHAP enediols without using a product
    # geometry. Both alpha sides are retained when a ketone has two choices.
    for oxygen, carbonyl_c, alpha_c, alpha_o, alpha_h in matches(
            '[O;+0]=[C;+0]-[C;H1,H2;+0](-[O;+0]-[H])'):
        alpha_hydrogens = [a.GetIdx() for a in mol.GetAtomWithIdx(alpha_c).GetNeighbors()
                           if a.GetAtomicNum() == 1]
        for proton in alpha_hydrogens:
            yield action('alpha_hydroxy_carbonyl_to_enediol',
                [((alpha_c,proton), (carbonyl_c,alpha_c)),
                 ((carbonyl_c,oxygen), (oxygen,proton))])
    for left_o, left_h, left_c, right_c, right_o, right_h in matches(
            '[O;+0](-[H])-[C;+0]=[C;+0]-[O;+0]-[H]'):
        yield action('enediol_to_left_carbonyl',
            [((left_o,left_h), (left_o,left_c)),
             ((left_c,right_c), (right_c,left_h))])
        yield action('enediol_to_right_carbonyl',
            [((right_o,right_h), (right_o,right_c)),
             ((left_c,right_c), (left_c,right_h))])
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
        for nu, _ in matches('[N;+0;v3]-[H]'):
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
    # Generic neutral phosphate monoester substitution.  In the enzyme-free
    # glycolysis pilot H3PO4/H2O are explicit atom-conserving surrogates; these
    # actions do not claim to reproduce ATP/ADP catalysis.
    for carbon, ester_o, phosphorus, phosphoryl_o, acid_o1, acid_h1, acid_o2, acid_h2 in phosphate_esters:
        for water_o, water_h1, water_h2 in waters:
            if component[ester_o] == component[water_o]:continue
            yield action('phosphate_monoester_hydrolysis',
                [((water_o,), (water_o,phosphorus)),
                 ((phosphorus,ester_o), (ester_o,)),
                 ((water_o,water_h1), (water_o,)),
                 ((ester_o,), (ester_o,water_h1))])
        # A neighboring alcohol can exchange positions with the ester oxygen.
        attached_carbon=carbon
        for neighbour in mol.GetAtomWithIdx(attached_carbon).GetNeighbors():
            if neighbour.GetAtomicNum()!=6:continue
            for oxygen in neighbour.GetNeighbors():
                if oxygen.GetAtomicNum()!=8 or oxygen.GetIdx()==ester_o:continue
                hydrogens=[a.GetIdx() for a in oxygen.GetNeighbors() if a.GetAtomicNum()==1]
                if not hydrogens:continue
                acceptor_o=oxygen.GetIdx();proton=hydrogens[0]
                yield action('vicinal_phosphate_migration',
                    [((acceptor_o,), (acceptor_o,phosphorus)),
                     ((phosphorus,ester_o), (ester_o,)),
                     ((acceptor_o,proton), (acceptor_o,)),
                     ((ester_o,), (ester_o,proton))])
    for phosphorus, phosphoryl_o, acid_o1, acid_h1, acid_o2, acid_h2, acid_o3, acid_h3 in phosphoric_acids:
        for alcohol_c, alcohol_o, alcohol_h in matches('[C;+0]-[O;+0]-[H]'):
            if component[phosphorus] == component[alcohol_o]:continue
            for leaving_o in (acid_o1,acid_o2,acid_o3):
                yield action('phosphate_monoester_condensation',
                    [((alcohol_o,), (alcohol_o,phosphorus)),
                     ((phosphorus,leaving_o), (leaving_o,)),
                     ((alcohol_o,alcohol_h), (alcohol_o,)),
                     ((leaving_o,), (leaving_o,alcohol_h))])
    # Site-specific glucose C6 phosphorylation. Internal secondary alcohols are
    # deliberately excluded from this reconstruction action.
    for c1,o1,c2,o2,h2,c3,o3,h3,c4,o4,h4,c5,o5,h5,c6,o6,h6 in matches(
            '[C;H1;+0](=[O;+0])-[C;H1;+0](-[O;+0]-[H])-[C;H1;+0](-[O;+0]-[H])-[C;H1;+0](-[O;+0]-[H])-[C;H1;+0](-[O;+0]-[H])-[C;H2;+0]-[O;+0]-[H]'):
        for phosphorus, phosphoryl_o, anion_o, acid_o1, acid_h1, acid_o2, acid_h2 in matches(
                '[P;+0](=[O;+0])(-[O;-1])(-[O;+0]-[H])-[O;+0]-[H]'):
            if component[phosphorus] == component[o6]:continue
            for leaving_o in (acid_o1,acid_o2):
                yield action('glucose_c6_phosphate_condensation',
                    [((o6,), (o6,phosphorus)),
                     ((phosphorus,leaving_o), (leaving_o,)),
                     ((o6,h6), (o6,)),
                     ((leaving_o,), (leaving_o,h6))])
    # Site-specific 3PG -> 2PG migration before the generic charged action.
    for carboxyl_c, carbonyl_o, acid_o, acid_h, c2, o2, h2, c3, o3, phosphorus, phosphoryl_o, anion_o, phosphate_o, phosphate_h in matches(
            '[C;+0](=[O;+0])(-[O;+0]-[H])-[C;H1;+0](-[O;+0]-[H])-[C;H2;+0]-[O;+0]-[P;+0](=[O;+0])(-[O;-1])-[O;+0]-[H]'):
        yield action('phosphoglycerate_3_to_2_site_specific',
            [((o2,), (o2,phosphorus)),
             ((phosphorus,o3), (o3,)),
             ((o2,h2), (o2,)),
             ((o3,), (o3,h2))])
    # Site-specific F6P/FBP pair. FBP has two terminal phosphate esters; only
    # hydrolysis at the CH2 group adjacent to the C2 carbonyl yields F6P.
    # Keep these before the generic charged actions so the named site-specific
    # proposal survives equivalent-arrow deduplication.
    for phosphorus, phosphoryl_o, anion_o, acid_o, acid_h, ester_o, c1, c2, carbonyl_o in matches(
            '[P;+0](=[O;+0])(-[O;-1])(-[O;+0]-[H])-[O;+0]-[C;H2;+0]-[C;+0](=[O;+0])'):
        for water_o, water_h1, water_h2 in waters:
            if component[ester_o] == component[water_o]:continue
            yield action('fbp_c1_phosphate_hydrolysis_to_f6p',
                [((water_o,), (water_o,phosphorus)),
                 ((phosphorus,ester_o), (ester_o,)),
                 ((water_o,water_h1), (water_o,)),
                 ((ester_o,), (ester_o,water_h1))])
    for alcohol_o, alcohol_h, c1, c2, carbonyl_o in matches(
            '[O;+0](-[H])-[C;H2;+0]-[C;+0](=[O;+0])'):
        for phosphorus, phosphoryl_o, anion_o, acid_o1, acid_h1, acid_o2, acid_h2 in matches(
                '[P;+0](=[O;+0])(-[O;-1])(-[O;+0]-[H])-[O;+0]-[H]'):
            if component[phosphorus] == component[alcohol_o]:continue
            for leaving_o in (acid_o1,acid_o2):
                yield action('f6p_c1_phosphate_condensation_to_fbp',
                    [((alcohol_o,), (alcohol_o,phosphorus)),
                     ((phosphorus,leaving_o), (leaving_o,)),
                     ((alcohol_o,alcohol_h), (alcohol_o,)),
                     ((leaving_o,), (leaving_o,alcohol_h))])
    # Hydrate the metaphosphate endpoint actually observed after C1 phosphate
    # cleavage. One water is already present in the fixed inventory; no extra
    # solvent molecule is introduced.
    for phosphorus, acceptor_o1, acceptor_o2, anion_o in matches(
            '[P;+0](=[O;+0])(=[O;+0])-[O;-1]'):
        for water_o, water_h1, water_h2 in waters:
            if component[phosphorus] == component[water_o]:continue
            for acceptor_o in (acceptor_o1,acceptor_o2):
                yield action('metaphosphate_water_hydration_to_h2po4',
                    [((water_o,), (water_o,phosphorus)),
                     ((phosphorus,acceptor_o), (acceptor_o,water_h1)),
                     ((water_o,water_h1), (water_o,))])
    # Singly deprotonated phosphate microstates used by the charged glycolysis
    # benchmark. The same atom-conserving substitution connects H2PO4- plus an
    # alcohol to a monoanionic phosphate ester plus water.
    for phosphorus, phosphoryl_o, anion_o, acid_o1, acid_h1, acid_o2, acid_h2 in matches(
            '[P;+0](=[O;+0])(-[O;-1])(-[O;+0]-[H])-[O;+0]-[H]'):
        for alcohol_c, alcohol_o, alcohol_h in matches('[C;+0]-[O;+0]-[H]'):
            if component[phosphorus] == component[alcohol_o]:continue
            for leaving_o in (acid_o1,acid_o2):
                yield action('phosphate_anion_monoester_condensation',
                    [((alcohol_o,), (alcohol_o,phosphorus)),
                     ((phosphorus,leaving_o), (leaving_o,)),
                     ((alcohol_o,alcohol_h), (alcohol_o,)),
                     ((leaving_o,), (leaving_o,alcohol_h))])
    for carbon, ester_o, phosphorus, phosphoryl_o, anion_o, acid_o, acid_h in matches(
            '[C;+0]-[O;+0]-[P;+0](=[O;+0])(-[O;-1])-[O;+0]-[H]'):
        for water_o, water_h1, water_h2 in waters:
            if component[ester_o] == component[water_o]:continue
            yield action('phosphate_anion_monoester_hydrolysis',
                [((water_o,), (water_o,phosphorus)),
                 ((phosphorus,ester_o), (ester_o,)),
                 ((water_o,water_h1), (water_o,)),
                 ((ester_o,), (ester_o,water_h1))])
        for neighbour in mol.GetAtomWithIdx(carbon).GetNeighbors():
            if neighbour.GetAtomicNum()!=6:continue
            for oxygen in neighbour.GetNeighbors():
                if oxygen.GetAtomicNum()!=8 or oxygen.GetIdx()==ester_o:continue
                hydrogens=[a.GetIdx() for a in oxygen.GetNeighbors() if a.GetAtomicNum()==1]
                if not hydrogens:continue
                acceptor_o=oxygen.GetIdx();proton=hydrogens[0]
                yield action('vicinal_phosphate_anion_migration',
                    [((acceptor_o,), (acceptor_o,phosphorus)),
                     ((phosphorus,ester_o), (ester_o,)),
                     ((acceptor_o,proton), (acceptor_o,)),
                     ((ester_o,), (ester_o,proton))])
    # Direct, atom-conserving retro-aldol grammar for a beta-hydroxy carbonyl.
    # The beta-OH proton terminates the carbon fragment created by C-C cleavage.
    for carbonyl_c, carbonyl_o, alpha_c, alpha_o, alpha_h, beta_c, beta_o, beta_h in matches(
            '[C;+0](=[O;+0])-[C;H1;+0](-[O;+0]-[H])-[C;H1;+0]-[O;+0]-[H]'):
        yield action('beta_hydroxy_carbonyl_retro_aldol',
            [((alpha_c,beta_c), (alpha_c,beta_h)),
             ((beta_o,beta_h), (beta_o,beta_c))])
    # Reverse aldol from an alpha-hydroxy carbonyl donor and an aldehyde/ketone.
    for donor_c, donor_o, alpha_c, alpha_o, alpha_h in matches(
            '[C;+0](=[O;+0])-[C;H1,H2;+0]-[O;+0]-[H]'):
        for acceptor_c, acceptor_o in carbonyls:
            if component[donor_c] == component[acceptor_c]:continue
            carbon_h=[a.GetIdx() for a in mol.GetAtomWithIdx(alpha_c).GetNeighbors()
                      if a.GetAtomicNum()==1]
            for proton in carbon_h:
                yield action('alpha_hydroxy_carbonyl_aldol_addition',
                    [((alpha_c,proton), (alpha_c,acceptor_c)),
                     ((acceptor_c,acceptor_o), (acceptor_o,proton))])
    # Enzyme-free 2-phosphoglycerate dehydration and its exact hydration reverse.
    for carboxyl_c, carboxyl_o, acid_o, acid_h, alpha_c, phosphate_o, phosphorus, beta_c, beta_o, beta_h in matches(
            '[C;+0](=[O;+0])(-[O;+0]-[H])-[C;H1;+0](-[O;+0]-[P;+0])-[C;H2;+0]-[O;+0]-[H]'):
        alpha_h=[a.GetIdx() for a in mol.GetAtomWithIdx(alpha_c).GetNeighbors()
                 if a.GetAtomicNum()==1]
        for proton in alpha_h:
            yield action('phosphoglycerate_dehydration',
                [((alpha_c,proton), (alpha_c,beta_c)),
                 ((beta_c,beta_o), (beta_o,proton))])
    for carboxyl_c, carboxyl_o, acid_o, acid_h, alpha_c, phosphate_o, phosphorus, beta_c in matches(
            '[C;+0](=[O;+0])(-[O;+0]-[H])-[C;+0](-[O;+0]-[P;+0])=[C;H2;+0]'):
        for water_o, water_h1, water_h2 in waters:
            if component[alpha_c] == component[water_o]:continue
            yield action('phosphoenolpyruvate_hydration',
                [((alpha_c,beta_c), (alpha_c,water_h1)),
                 ((water_o,water_h1), (beta_c,water_o))])
    # Canonical formose closure.  Retro-aldol transfers the beta-OH proton
    # while cleaving the central C-C bond, yielding two neutral glycolaldehydes.
    for c1,o1,c2,o2,h2,c3,o3,h3,c4,o4,h4 in aldotetroses:
        yield action('formose_retro_aldol_tetrose_to_2go',
            [((c2,c3),(c2,h3)),((o3,h3),(o3,)),((o3,),(o3,c3))])
    # Exact reverse action from two glycolaldehydes. Symmetry-related choices
    # and the two alpha hydrogens remain distinct three-dimensional seeds.
    for acceptor in glycolaldehydes:
        for donor in glycolaldehydes:
            ca,oa,alpha,oh_a,h_oh_a=acceptor
            cd,od,alpha_d,oh_d,h_oh_d=donor
            if component[ca]==component[cd]:continue
            for atom in mol.GetAtomWithIdx(alpha).GetNeighbors():
                if atom.GetAtomicNum()!=1:continue
                h_alpha=atom.GetIdx()
                yield action('formose_inverse_aldol_2go_to_tetrose',
                    [((alpha,h_alpha),(alpha,cd)),((od,),(od,h_alpha)),((cd,od),(od,))])
    # Transamidation continuation from a neutral tetrahedral intermediate.
    # Identify the original NH2 leaving group and the carbon-substituted incoming
    # amine directly from the observed graph; no expected product geometry enters.
    for center in mol.GetAtoms():
        if center.GetAtomicNum()!=6:continue
        neighbours=list(center.GetNeighbors())
        oxygens=[a for a in neighbours if a.GetAtomicNum()==8 and a.GetFormalCharge()==0]
        nitrogens=[a for a in neighbours if a.GetAtomicNum()==7]
        if len(oxygens)!=1 or len(nitrogens)!=2:continue
        oxygen=oxygens[0]
        oxygen_h=[a for a in oxygen.GetNeighbors() if a.GetAtomicNum()==1]
        if len(oxygen_h)!=1:continue
        leaving=[];incoming=[]
        for nitrogen in nitrogens:
            heavy=[a for a in nitrogen.GetNeighbors() if a.GetAtomicNum()!=1 and a.GetIdx()!=center.GetIdx()]
            hydrogens=[a for a in nitrogen.GetNeighbors() if a.GetAtomicNum()==1]
            if not heavy and len(hydrogens)>=2:leaving.append((nitrogen,hydrogens))
            if len(heavy)==1:incoming.append((nitrogen,hydrogens))
        if len(leaving)!=1 or len(incoming)!=1:continue
        leave,leave_h=leaving[0];arrive,arrive_h=incoming[0]
        c,o,nl,na,ho=center.GetIdx(),oxygen.GetIdx(),leave.GetIdx(),arrive.GetIdx(),oxygen_h[0].GetIdx()
        if leave.GetFormalCharge()==0 and arrive.GetFormalCharge()==0 and arrive_h:
            hn=arrive_h[0].GetIdx()
            yield action('transamidation_amine_to_amine_proton_transfer',
                [((na,hn),(na,)),((nl,),(nl,hn))])
            yield action('transamidation_proton_coupled_amine_elimination',
                [((o,),(o,c)),((c,nl),(nl,ho)),((o,ho),(o,))])
        elif leave.GetFormalCharge()==1 and arrive.GetFormalCharge()==-1 and not arrive_h:
            yield action('transamidation_tetrahedral_collapse',
                [((o,),(o,c)),((c,nl),(nl,)),((o,ho),(o,)),((na,),(na,ho))])


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
