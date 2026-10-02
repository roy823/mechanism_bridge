"""Prepare targeted roots for the four missing sparse-endpoint glycolysis segments."""
import json
from pathlib import Path
import sys

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
from mechbridge.event_graph import graph_smiles
from mechbridge.symbolic_library import parse_explicit
from scripts.data.prepare_glycolysis_sparse_endpoints import SPECIES

SOURCE_STARTS = ROOT / 'data/processed/glycolysis_sparse_endpoint_starts.jsonl'
SOURCE_REPORT = ROOT / 'reports/glycolysis_sparse_endpoint_v17/runs'
OUTPUT = ROOT / 'data/processed/glycolysis_gap_recovery_starts.jsonl'


def skeleton_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def main():
    source_starts = [json.loads(line) for line in SOURCE_STARTS.read_text(encoding='utf-8').splitlines()]
    starts = []
    enol_target = skeleton_key(graph_smiles(parse_explicit('C=C(O)C(=O)O.' + SPECIES['pi'])))
    enol_nodes = []
    for path in sorted(SOURCE_REPORT.rglob('network.json')):
        network = json.loads(path.read_text(encoding='utf-8'))
        if network['start']['provenance']['system'] != 'pep_hydrolysis':
            continue
        for node in network['nodes']:
            if skeleton_key(node['graph_smiles']) == enol_target:
                enol_nodes.append((path, network, node))
    if not enol_nodes:
        raise ValueError('No observed enol-pyruvate endpoint')
    selected_enols = (enol_nodes[:2] if len(enol_nodes) >= 2 else enol_nodes * 2)
    for index, (path, network, node) in enumerate(selected_enols):
        scale = (.8, 1.2)[index]
        starts.append(dict(id=f'pep_enol_continuation_c{index}',
            atomic_numbers=network['start']['atomic_numbers'], positions_A=node['positions_A'],
            charge=network['start']['charge'], multiplicity=network['start']['multiplicity'],
            provenance=dict(recovery_group='pep_enol_to_pyruvate', seed_scale=scale,
                source_network=path.relative_to(ROOT).as_posix(), source_node=node['id'],
                input_smiles=node['graph_smiles'], actual_descent_endpoint=True,
                target_product_geometry_used=False)))
    copies = {
        'glucose_phosphorylation': [('glucose_c6', 1.0)],
        'phosphoglycerate_shift': [('pg_shift_soft', .7), ('pg_shift_strong', 1.3)],
        'fbp_cleavage': [('fbp_cleavage_strong', 1.4), ('fbp_cleavage_very_strong', 1.8)],
    }
    for source in source_starts:
        system = source['provenance']['system']
        if system not in copies:
            continue
        conformer = source['id'].rsplit('_c', 1)[-1]
        for label, scale in copies[system]:
            row = json.loads(json.dumps(source))
            row['id'] = f'{label}_c{conformer}'
            row['provenance'] = dict(row['provenance'], recovery_group=label,
                seed_scale=scale, source_start=source['id'], target_product_geometry_used=False)
            starts.append(row)
    OUTPUT.write_text(''.join(json.dumps(row) + '\n' for row in starts), encoding='utf-8')
    print(json.dumps(dict(starts=len(starts), groups={group:sum(
        row['provenance']['recovery_group']==group for row in starts)
        for group in sorted({row['provenance']['recovery_group'] for row in starts})}), indent=2))


if __name__ == '__main__':
    main()
