"""Export observed F6P + water + metaphosphate minima as hydration roots."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NETWORK = (ROOT / 'reports/glycolysis_charged_sparse_v16/recovery_multistep/'
           'fbp_water_c1/fbp_water_c1/arrows/network.json')
OUTPUT = ROOT / 'data/processed/glycolysis_metaphosphate_continuations.jsonl'


def main():
    network = json.loads(NETWORK.read_text(encoding='utf-8'))
    rows = []
    for node in network['nodes']:
        if 'O=P(=O)[O-]' not in node['graph_smiles']:
            continue
        rows.append(dict(id=f"f6p_water_metaphosphate_n{node['id']}",
            atomic_numbers=network['start']['atomic_numbers'], positions_A=node['positions_A'],
            charge=network['start']['charge'], multiplicity=network['start']['multiplicity'],
            provenance=dict(system='f6p_water_metaphosphate', source_network=NETWORK.relative_to(ROOT).as_posix(),
                source_node=node['id'], input_smiles=node['graph_smiles'],
                role='observed two-sided-descent endpoint; hydration continuation',
                added_small_molecules=0, reference_TS_or_product_geometry_used=False)))
    OUTPUT.write_text(''.join(json.dumps(row) + '\n' for row in rows), encoding='utf-8')
    print(json.dumps(dict(starts=len(rows), ids=[row['id'] for row in rows]), indent=2))


if __name__ == '__main__':
    main()
