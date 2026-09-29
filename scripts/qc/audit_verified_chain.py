"""Compare actual DFT endpoints of the two preselected network edges."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol
from mechbridge.reaction_network import molecular_rmsd,atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    root=p.parse_args().root.resolve()
    selection=json.loads((root/'DFT_selection.json').read_text())
    if selection is None:raise ValueError('No chain selected')
    network_path=ROOT/selection['run']/selection['start']/selection['strategy']/'network.json'
    network=json.loads(network_path.read_text(encoding='utf-8'))
    records=[];events=[];feedback=[];shared=[]
    for edge in selection['edges_to_verify']:
        folder=root/'qc'/f"{selection['start']}_edge{edge}"
        path=folder/'verification.json';q=json.loads(path.read_text())
        if q['status'] in ('started','refining_source_ts','integrating_irc','computing_electronic_path'):
            raise ValueError('DFT calculation still running')
        records.append((folder,q))
        source_edge=next(e for e in network['edges'] if e['id']==edge)
        source_attempt=network['attempts'][source_edge['attempt']]
        proposing_node=network['nodes'][source_attempt['source_node']]
        feedback.append(dict(event_id=q['event_id'],MLIP_edge_id=edge,
            status=q['status'],physical_event_verified=q.get('physical_event_verified',False),
            original_pair_preserved=q.get('expected_endpoint_match'),
            actual_graphs=[e['graph_smiles'] for e in q.get('endpoints',[])],
            model_weights_updated=False,source=str(path.relative_to(ROOT)),
            meaning='Event-specific physical evidence; failure is not proof of reaction impossibility'))
        if q.get('physical_event_verified'):
            ts=read(folder/'ts.xyz')
            events.append(dict(event_id=q['event_id'],atomic_numbers=ts.numbers.tolist(),
                positions_A=dict(ts=ts.positions.tolist(),endpoints=[
                    read(folder/f'minimum_{side}.xyz').positions.tolist() for side in ('forward','reverse')]),
                endpoint_graphs=[e['graph_smiles'] for e in q['endpoints']],
                method=q['method'],basis=q['basis'],charge=q['charge'],multiplicity=q['multiplicity'],
                physical_event_verified=True,verified_pair=False,is_IRC=True,
                original_MLIP_pair_preserved=q.get('expected_endpoint_match'),
                proposed_from_graph=proposing_node['graph_smiles'],
                proposed_from_positions_A=proposing_node['positions_A'],
                source_network=str(network_path.relative_to(ROOT)),
                symbolic_status='proposal_only_not_independently_reviewed',
                symbolic_proposal=q['source']['original_symbolic_proposal'],
                source_verification=str(path.relative_to(ROOT))))
    both=all(q.get('physical_event_verified') for _,q in records)
    same_level=len({(q['method'],q['basis'],q['charge'],q['multiplicity'],q['environment']) for _,q in records})==1
    if both and same_level:
        (f1,q1),(f2,q2)=records
        for i,a in enumerate(q1['endpoints']):
            for j,b in enumerate(q2['endpoints']):
                if a['graph_smiles']!=b['graph_smiles']:continue
                x=read(f1/f"minimum_{('forward','reverse')[i]}.xyz")
                y=read(f2/f"minimum_{('forward','reverse')[j]}.xyz")
                rmsd=molecular_rmsd(geometry_mol(x.numbers,x.positions,0),x.positions,
                    geometry_mol(y.numbers,y.positions,0),y.positions)
                delta=abs(a['energy_eV']-b['energy_eV'])
                shared.append(dict(graph=a['graph_smiles'],RMSD_A=float(rmsd),energy_difference_eV=delta,
                                   matches_same_minimum=bool(rmsd<.15 and delta<.03)))
    distinct_graphs={g for row in feedback for g in row['actual_graphs']}
    actual_chain=bool(both and same_level and len(distinct_graphs)>=3 and any(s['matches_same_minimum'] for s in shared))
    verdict=dict(selection=selection,events=feedback,shared_intermediates=shared,
        actual_two_step_DFT_chain=actual_chain,
        original_MLIP_chain_preserved=bool(actual_chain and all(q.get('expected_endpoint_match') for _,q in records)),
        model_weights_updated=False,independent_arrow_gold=False)
    atomic_json(root/'DFT_chain_audit.json',verdict)
    for name,rows in [('DFT_events',events),('DFT_feedback',feedback)]:
        (root/(name+'.jsonl')).write_text(''.join(json.dumps(r)+'\n' for r in rows),encoding='utf-8')
    print(json.dumps(verdict,indent=2))


if __name__=='__main__':main()
