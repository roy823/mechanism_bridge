"""Export the observed pyruvate anion + phosphoric acid endpoint."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NETWORK = (ROOT / 'reports/glycolysis_gap_recovery_v18/runs/pep_enol_continuation_c0/'
           'pep_enol_continuation_c0/arrows/network.json')
OUTPUT = ROOT / 'data/processed/pyruvate_protonation_continuations.jsonl'


def main():
    network = json.loads(NETWORK.read_text(encoding='utf-8'))
    node = next(node for node in network['nodes']
                if 'CC(=O)C(=O)[O-].O=P(O)(O)O' in node['graph_smiles'])
    rows = []
    for index, scale in enumerate((.8, 1.2)):
        rows.append(dict(id=f'pyruvate_protonation_c{index}',
            atomic_numbers=network['start']['atomic_numbers'], positions_A=node['positions_A'],
            charge=network['start']['charge'], multiplicity=network['start']['multiplicity'],
            provenance=dict(system='pyruvate_protonation', seed_scale=scale,
                source_network=NETWORK.relative_to(ROOT).as_posix(), source_node=node['id'],
                actual_descent_endpoint=True, target_product_geometry_used=False)))
    OUTPUT.write_text(''.join(json.dumps(row) + '\n' for row in rows), encoding='utf-8')
    print(json.dumps(dict(starts=len(rows), source_node=node['id']), indent=2))


if __name__ == '__main__':
    main()
