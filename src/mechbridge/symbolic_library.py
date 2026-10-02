"""Replay and transfer published two-electron arrows onto an observed reactant.

No reference TS, product geometry, or measured barrier is accepted by this API.
Local patterns retain charge, valence and one heavy-atom context shell.
"""
from collections import Counter
import json
from rdkit import Chem
from .event_graph import graph_smiles, bond_orders
from .local_patterns import compile_pattern, matches_pattern
from .intermolecular_actions import electron_actions, crosses_components, diverse_proposals


def parse_explicit(smiles):
    p = Chem.SmilesParserParams()
    p.removeHs = False
    mol = Chem.MolFromSmiles(smiles, p)
    if mol is None:
        raise ValueError('Cannot parse reaction graph')
    return Chem.AddHs(mol)


def replay(mol, arrows):
    """Apply all arrows simultaneously; formal electron ownership fixes charges."""
    delta = Counter()
    charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
    old = bond_orders(mol)
    outgoing = Counter(tuple(sorted(a['source'])) for a in arrows)
    for site, count in outgoing.items():
        if len(site) == 1:
            i = site[0]
            if i < 0 or i >= mol.GetNumAtoms():
                raise ValueError('Atom index outside molecule')
            atom = mol.GetAtomWithIdx(i)
            valence = Chem.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())
            capacity = (valence - atom.GetFormalCharge() - atom.GetTotalValence()) / 2
        elif len(site) == 2:
            capacity = old.get(site, 0)
        else:
            raise ValueError('Invalid electron source')
        if count > capacity:
            raise ValueError('Arrow source lacks the requested occupied electron pairs')
    for arrow in arrows:
        if arrow['electrons'] != 2:
            raise ValueError('Only paired-electron arrows are supported')
        for sign, site in [(-1, arrow['source']), (1, arrow['sink'])]:
            if len(site) not in (1, 2) or len(set(site)) != len(site):
                raise ValueError('Invalid electron site')
            if min(site) < 0 or max(site) >= mol.GetNumAtoms():
                raise ValueError('Atom index outside molecule')
            if len(site) == 2:
                delta[tuple(sorted(site))] += sign
            for i in site:
                charges[i] -= sign * (2 // len(site))
    out = Chem.RWMol(mol)
    types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    edits = []
    for (i, j), change in sorted(delta.items()):
        if not change:
            continue
        before = old.get((i, j), 0)
        after = before + change
        if after not in (0, 1, 2, 3) or before not in (0, 1, 2, 3):
            raise ValueError('Unsupported bond order in arrow replay')
        if before:
            out.RemoveBond(i, j)
        if after:
            out.AddBond(i, j, types[after])
        edits.append(dict(atoms=[i, j], before=before, after=after))
    for i, q in enumerate(charges):
        out.GetAtomWithIdx(i).SetFormalCharge(q)
    out = out.GetMol()
    Chem.SanitizeMol(out)
    if any(a.GetNumRadicalElectrons() for a in out.GetAtoms()):
        raise ValueError('Replay generated an open-shell graph')
    if Chem.GetFormalCharge(out) != Chem.GetFormalCharge(mol):
        raise ValueError('Electron replay changed total charge')
    return out, edits


class ArrowLibrary:
    policy = 'published_local_arrows_plus_explicit_reactant_grammar_v8_charged_CHNOP'

    def __init__(self, path):
        self.templates = []
        self.source_graphs = set()
        self.audit = Counter()
        self.rejections = []
        data = json.loads(open(path, encoding='utf-8').read())
        for rec in data['records']:
            self.audit['source_records'] += 1
            try:
                r, p = [parse_explicit(s) for s in rec['rsmi'].split('>>')]
                if (any(a.GetAtomicNum() not in (1, 6, 7, 8) for a in r.GetAtoms())
                    or Chem.GetFormalCharge(r) != 0
                    or Chem.GetFormalCharge(p) != 0
                    or len(Chem.GetMolFrags(r)) > 2
                    or any(a.GetNumRadicalElectrons() for a in r.GetAtoms())):
                    self.audit['outside_pilot_domain'] += 1
                    continue
                maps = {a.GetAtomMapNum(): a.GetIdx() for a in r.GetAtoms() if a.GetAtomMapNum()}
                arrows = [dict(source=[maps[i] for i in src], sink=[maps[i] for i in dst],
                               electrons=2, kind=kind) for kind, src, dst in rec['epd']]
                product, edits = replay(r, arrows)
                if graph_smiles(product) != graph_smiles(p):
                    raise ValueError('Published arrow replay differs from recorded product')
                if not edits or graph_smiles(r) == graph_smiles(product):
                    self.audit['no_distinct_chemical_graph'] += 1
                    continue
                for direction, mol, arr in [('forward', r, arrows), ('reverse', product,
                        [dict(a, source=a['sink'], sink=a['source']) for a in arrows])]:
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    self.source_graphs.add(graph_smiles(mol))
                    try:
                        replay(mol, arr)
                        query, local_arrows = compile_pattern(mol, arr)
                    except (ValueError, RuntimeError) as exc:
                        self.audit['unsupported_local_directions'] += 1
                        self.rejections.append(dict(source_id=rec['id'],direction=direction,
                                                    stage='local_pattern',reason=str(exc)))
                        continue
                    self.templates.append(dict(template_id=f"SynEPD:{rec['id']}:{direction}",
                        query=query, arrows=local_arrows, source_graph=graph_smiles(mol),
                        smarts=Chem.MolToSmarts(query), name=rec['reaction_name'], source_id=rec['id']))
                self.audit['replay_valid_records'] += 1
                self.audit['bimolecular_records' if len(Chem.GetMolFrags(r)) == 2
                           else 'unimolecular_records'] += 1
            except (ValueError, KeyError, RuntimeError) as exc:
                self.audit['replay_or_parse_rejected'] += 1
                self.rejections.append(dict(source_id=rec['id'],stage='source_replay',reason=str(exc)))

    def propose(self, mol, limit=24):
        if limit < 1:
            raise ValueError('Proposal limit must be positive')
        if any(a.GetAtomicNum() not in (1, 6, 7, 8, 15) or a.GetNumRadicalElectrons()
               or a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in mol.GetAtoms()):
            raise ValueError('Local pilot requires explicit-H closed-shell CHNOP')
        key = graph_smiles(mol)
        result, seen = [], set()
        for template in self.templates:
            query = template['query']
            # Keep symmetry-related atom assignments as geometrically different proposals.
            for match in mol.GetSubstructMatches(query, uniquify=False, useChirality=True, maxMatches=64):
                if not matches_pattern(mol, query, match):
                    continue
                arrows = [dict(a, source=[match[i] for i in a['source']],
                               sink=[match[i] for i in a['sink']]) for a in template['arrows']]
                try:
                    product, edits = replay(mol, arrows)
                    product_graph = graph_smiles(product)
                except (ValueError, RuntimeError):
                    continue
                if not edits or product_graph == key:
                    continue
                # Different electron source/sink assignments can have identical net
                # edits. Preserve these alternatives for the cross-representation task.
                signature = (tuple((tuple(e['atoms']), e['after']) for e in edits),
                             tuple(sorted((tuple(sorted(a['source'])),tuple(sorted(a['sink'])),a['electrons'])
                                          for a in arrows)))
                if signature in seen:
                    continue
                seen.add(signature)
                result.append(dict(template_id=template['template_id'], name=template['name'],
                                   arrows=arrows, edits=edits, predicted_graph=product_graph,
                                   origin='local_transfer_of_published_SynEPD_arrows',
                                   template_source_graph=template['source_graph'],
                                   transferred=key != template['source_graph'],
                                   pattern_smarts=template['smarts'], matched_atoms=list(match),
                                   intermolecular=crosses_components(mol, edits),
                                   stereo_policy='No new stereochemistry inferred by symbolic replay'))
        for name, arrows in electron_actions(mol):
            try:
                product, edits = replay(mol, arrows)
                product_graph = graph_smiles(product)
            except (ValueError, RuntimeError):
                continue
            signature = (tuple((tuple(e['atoms']), e['after']) for e in edits),
                         tuple(sorted((tuple(sorted(a['source'])), tuple(sorted(a['sink'])), 2)
                                      for a in arrows)))
            if not edits or signature in seen or product_graph == key:
                continue
            seen.add(signature)
            result.append(dict(template_id='grammar:'+name, name=name, arrows=arrows,
                edits=edits, predicted_graph=product_graph,
                origin='analyst_defined_reactant_only_electron_action_grammar',
                independently_annotated=False, intermolecular=crosses_components(mol,edits),
                stereo_policy='No new stereochemistry inferred by symbolic replay'))
        return diverse_proposals(result, limit)


class FilteredArrowLibrary:
    """Expose only explicitly requested reviewed actions for targeted searches."""
    def __init__(self,library,template_ids):
        self.library=library;self.template_ids=tuple(template_ids)
        if not self.template_ids:raise ValueError('At least one template ID is required')
        self.policy=library.policy+'; filter='+','.join(self.template_ids)
        self.audit=library.audit

    def propose(self,mol,limit=24):
        # Filtering must happen before truncation: large polyfunctional molecules
        # can otherwise fill the diverse proposal budget before the requested
        # reviewed action is reached.
        proposals=self.library.propose(mol,limit=max(1024,limit))
        return [proposal for proposal in proposals
                if proposal['template_id'] in self.template_ids][:limit]


class ResonanceAwareLibrary:
    """Also propose from enumerated resonance forms of the node graph.

    Geometry perception picks one Lewis structure, but templates and grammar
    rules are written for particular forms (e.g. a nitrile ylide rather than its
    allenyl-enolate form), so a reaction can be missed only because of the
    drawing. Resonance forms keep atom order and connectivity, so their arrows
    and edits apply to the same geometry. Proposals are deduplicated by edits
    and arrows and record the form they came from.
    """
    def __init__(self, library, max_forms=8):
        if max_forms < 1:
            raise ValueError('max_forms must be positive')
        self.library, self.max_forms = library, int(max_forms)
        self.policy = library.policy + f'; resonance_forms={self.max_forms}'
        self.audit = library.audit

    def forms(self, mol):
        flags = (Chem.UNCONSTRAINED_ANIONS | Chem.UNCONSTRAINED_CATIONS |
                 Chem.ALLOW_CHARGE_SEPARATION)
        key, out = graph_smiles(mol), []
        try:
            structures = Chem.ResonanceMolSupplier(mol, flags=flags, maxStructs=4*self.max_forms)
        except RuntimeError:
            return out
        seen = {key}
        for form in structures:
            if form is None or len(out) >= self.max_forms:
                continue
            try:
                form = Chem.Mol(form)
                Chem.SanitizeMol(form)
                smiles = graph_smiles(form)
            except (ValueError, RuntimeError):
                continue
            if smiles not in seen and Chem.GetFormalCharge(form) == Chem.GetFormalCharge(mol):
                seen.add(smiles)
                out.append(form)
        return out

    def propose(self, mol, limit=24):
        base = self.library.propose(mol, limit=max(limit, 64))
        signature = lambda p: (tuple((tuple(e['atoms']), e['after']) for e in p['edits']),
                               tuple(sorted((tuple(sorted(a['source'])), tuple(sorted(a['sink'])))
                                            for a in p['arrows'])))
        seen = {signature(p) for p in base}
        result = list(base)
        forms = self.forms(mol)
        # A "product" that is only another drawing of the node is not a reaction.
        drawings = {graph_smiles(mol)} | {graph_smiles(f) for f in forms}
        for form in forms:
            try:
                proposals = self.library.propose(form, limit=max(limit, 64))
            except ValueError:
                continue
            for p in proposals:
                s = signature(p)
                if s in seen or p['predicted_graph'] in drawings:
                    continue
                seen.add(s)
                result.append(dict(p, origin=p.get('origin', '')+'_via_resonance_form',
                                   resonance_form=graph_smiles(form)))
        return diverse_proposals(result, limit)
