"""Give events from repeated campaigns unique IDs inside a molecular report."""
import argparse
from collections import Counter
import json
from pathlib import Path


def campaign_tag(source):
    if 'glycolysis_relaxed_recovery_v15_resume' in source:
        return 'v15r'
    if 'glycolysis_relaxed_recovery_v15' in source:
        return 'v15'
    if 'glycolysis_targeted_refine_v14' in source:
        return 'v14t'
    return 'v14'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html',type=Path)
    args=parser.parse_args()
    html=args.html.resolve();text=html.read_text(encoding='utf-8')
    marker='const DATA=';start=text.index(marker)+len(marker)
    payload,consumed=json.JSONDecoder().raw_decode(text[start:])
    events=payload['events'];edges=payload['networks'][0]['edges']
    if len(events)!=len(edges):
        raise ValueError('Event/edge ordering invariant is unavailable')
    before=len(events)-len({event['id'] for event in events});seen=Counter()
    for event,edge in zip(events,edges):
        base=f"{campaign_tag(event.get('source_result',''))}_{event['id']}"
        occurrence=seen[base];seen[base]+=1
        identifier=base if occurrence==0 else f"{base}_{occurrence}"
        event['id']=identifier;event['aggregate_edge_id']=identifier
        event['downloads']=dict(path=identifier+'_path.xyz',ts=identifier+'_TS.xyz')
        edge['id']=identifier;edge['event_id']=identifier
    encoded=json.dumps(payload,ensure_ascii=False,separators=(',',':'))
    html.write_text(text[:start]+encoded+text[start+consumed:],encoding='utf-8')
    network_path=html.parent/'global_network.json'
    network=json.loads(network_path.read_text(encoding='utf-8'))
    if len(network['edges'])!=len(edges):
        raise ValueError('Global network edge count differs from HTML payload')
    for target,source in zip(network['edges'],edges):
        target['id']=source['id'];target['event_id']=source['event_id']
    network_path.write_text(json.dumps(network,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(events=len(events),duplicate_ids_before=before,
                          duplicate_ids_after=len(events)-len({e['id'] for e in events}))))


if __name__=='__main__':
    main()
