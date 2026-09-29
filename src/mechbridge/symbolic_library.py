"""Replay and transfer published two-electron arrows onto an observed reactant.

No reference TS, product geometry, or measured barrier is accepted by this API.
Whole-system graph matches are deliberately conservative; no template invention.
"""
from collections import Counter, defaultdict
import json
from rdkit import Chem
from .event_graph import graph_smiles, bond_orders


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
    old = bond_orders(mol)
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
    return out, edits


class ArrowLibrary:
    def __init__(self, path):
        self.by_reactant = defaultdict(list)
        self.audit = Counter()
        data = json.loads(open(path, encoding='utf-8').read())
        for rec in data['records']:
            self.audit['source_records'] += 1
            try:
                r, p = [parse_explicit(s) for s in rec['rsmi'].split('>>')]
                if (any(a.GetAtomicNum() not in (1, 6, 7, 8) for a in r.GetAtoms())
                    or Chem.GetFormalCharge(r) != 0
                    or Chem.GetFormalCharge(p) != 0
                    or len(Chem.GetMolFrags(r)) != 1
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
                    self.by_reactant[graph_smiles(mol)].append(dict(
                        template_id=f"SynEPD:{rec['id']}:{direction}", mol=mol, arrows=arr,
                        name=rec['reaction_name'], source_id=rec['id']))
                self.audit['replay_valid_records'] += 1
            except (ValueError, KeyError, RuntimeError):
                self.audit['replay_or_parse_rejected'] += 1

    def propose(self, mol, limit=24):
        key = graph_smiles(mol)
        result, seen = [], set()
        for template in self.by_reactant.get(key, []):
            query = template['mol']
            if query.GetNumAtoms() != mol.GetNumAtoms():
                continue
            # Keep symmetry-related atom assignments as geometrically different proposals.
            for match in mol.GetSubstructMatches(query, uniquify=False, useChirality=True, maxMatches=64):
                arrows = [dict(a, source=[match[i] for i in a['source']],
                               sink=[match[i] for i in a['sink']]) for a in template['arrows']]
                try:
                    product, edits = replay(mol, arrows)
                except (ValueError, RuntimeError):
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
                                   arrows=arrows, edits=edits, predicted_graph=graph_smiles(product),
                                   origin='published_SynEPD_arrows_replayed_on_observed_state'))
                if len(result) >= limit:
                    return result
        return result
