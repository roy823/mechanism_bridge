import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load():
    path = ROOT/'scripts/diagnostics/audit_benchmark_overlap.py'
    spec = importlib.util.spec_from_file_location('audit_benchmark_overlap', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_overlap_labels_follow_the_key_hierarchy(tmp_path):
    audit = load()
    rgd1 = tmp_path/'rgd1.csv'
    rgd1.write_text('reaction,reactant,product\n'
                    'MR_1,[CH3:1][CH2:2][OH:3],[CH3:1][CH:2]=[O:3].[H][H]\n'
                    'MR_2,C[C@@H]1CO1,CCC=O\n')
    index, stats = audit.rgd1_index(rgd1)
    assert stats == dict(rows=2, parsed=2, parse_errors=0)
    cases = [('seen_reaction', 'CC=O.[H][H]', 'CCO'),      # reversed direction, maps dropped
             ('seen_reaction', 'CC1CO1', 'CCC=O'),         # stereo dropped
             ('seen_reactant', 'CCO', 'C=C.O'),
             ('seen_formula', 'COC', 'C=O.C'),
             ('unseen', 'CCCCO', 'CCCC=O.[H][H]')]
    for expected, reactant, product in cases:
        result = audit.label(dict(reactant_key=reactant, product_key=product), index)
        assert result['label'] == expected, (expected, reactant, result)
    assert audit.label(dict(reactant_key='CC1CO1', product_key='CCC=O'), index)['k1_rgd1_ids'] == ['MR_2']
