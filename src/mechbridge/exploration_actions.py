"""Finite, diverse action scheduling on observed minima; no imagined graph edges."""
import json


def action_key(proposal):
    return json.dumps(dict(edits=proposal['edits'],arrows=proposal['arrows']),sort_keys=True)


def choose_action(proposals,used,known_graphs,strategy,geometry_limit=9,symbolic_variant_limit=3):
    """Spread trials across product hypotheses before symmetry-related mappings.

    Each atom-mapped action has three geometric variants. Successful observations
    update known_graphs; this is adaptive scheduling, not neural-model training.
    """
    geometry_count=sum(n for key,n in used.items() if key=='geometry')
    trials=sum(used.values())
    symbolic=[]
    product_trials={}
    for p in proposals:
        count=used.get(action_key(p),0)
        product_trials[p['predicted_graph']]=product_trials.get(p['predicted_graph'],0)+count
    for index,p in enumerate(proposals):
        key=action_key(p);count=used.get(key,0)
        if count<symbolic_variant_limit:
            symbolic.append(((p['predicted_graph'] in known_graphs,
                product_trials[p['predicted_graph']],count,index),p,key,count))
    geometry_due=(strategy=='geometry' or (strategy=='hybrid' and (trials%4==3 or not symbolic)))
    if geometry_due and geometry_count<geometry_limit:
        return None,'geometry',geometry_count,'geometry'
    if symbolic and strategy!='geometry':
        _,proposal,key,variant=min(symbolic,key=lambda v:v[0])
        return proposal,key,variant,'arrows' if strategy=='hybrid' else strategy
    return None
