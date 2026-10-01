"""Shared-state CPU scheduler for GPA-style reaction-network exploration."""
from concurrent.futures import FIRST_COMPLETED,ProcessPoolExecutor,wait
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import time

import networkx as nx
import numpy as np
from ase import Atoms
from rdkit import RDLogger

from .event_classification import classify_event
from .event_graph import geometry_mol,graph_smiles
from .exploration_actions import action_key,choose_action
from .potentials import load_potential
from .reaction_network import (BudgetExceeded,CountedCalculator,SearchProtocol,aligned_rmsd,
                               atomic_json,initialize_root,molecular_rmsd,search_connection)
from .search_seeds import make_seed
from .species_network import project_species_network


_WORKER_BACKEND=None


def _initialize_worker(potential_name,root,device,compile_model,threads):
    global _WORKER_BACKEND
    os.environ['OMP_NUM_THREADS']=str(threads)
    os.environ['MKL_NUM_THREADS']=str(threads)
    RDLogger.DisableLog('rdApp.*')
    _WORKER_BACKEND,_=load_potential(potential_name,root,device,compile_model,threads=threads)


def _run_task(task):
    """Worker boundary: generate one seed and return one physical search result."""
    backend=_WORKER_BACKEND
    validator=getattr(backend,'validate_system',None)
    if validator is not None:
        validator(task['numbers'],task['charge'],task['multiplicity'])
    protocol=SearchProtocol(**task['protocol'])
    destination=Path(task['destination'])
    state=Atoms(numbers=task['numbers'],positions=task['positions_A'])
    mol=geometry_mol(task['numbers'],task['positions_A'],task['charge'])
    try:
        started=time.perf_counter()
        x,direction,meta=make_seed(state,mol,task['seed_strategy'],task['proposal'],
                                   task['variant'],task['random_seed'])
        meta.update(random_seed=task['random_seed'],seed_strategy=task['seed_strategy'],
                    generation_seconds=time.perf_counter()-started)
    except (ValueError,RuntimeError) as exc:
        destination.mkdir(parents=True,exist_ok=False)
        result=dict(status='seed_generation_failed',error=str(exc),evaluations=0,seconds=0.,
                    evidence='MLIP_two_sided_mode_displacement_descent',is_IRC=False,
                    DFT_verified=False,endpoints=[])
        atomic_json(destination/'result.json',result)
        return dict(result=result,proposal=task['proposal'],meta=task['proposal'])
    calculator=CountedCalculator(backend,task['evaluation_limit'])
    trial=state.copy();trial.positions=x
    result=search_connection(trial,direction,calculator,destination,protocol,task['charge'])
    return dict(result=result,proposal=task['proposal'],meta=meta)


def reservation_identity(node,proposal,variant,seed_strategy,conformer_cluster='default'):
    payload=dict(species=node['graph_smiles'],conformer_cluster=conformer_cluster,
                 action=action_key(proposal) if proposal is not None else 'geometry',
                 variant=variant,seed_strategy=seed_strategy)
    encoded=json.dumps(payload,sort_keys=True,separators=(',',':'))
    return hashlib.sha256(encoded.encode()).hexdigest(),payload


def conformer_frontier(nodes,mols,tolerance_A=.5):
    """Cluster physical minima within species while retaining reactive conformers."""
    clusters=[]
    for node,mol in zip(nodes,mols):
        match=None
        for cluster in clusters:
            if cluster['graph_smiles']!=node['graph_smiles']:continue
            representative=cluster['representative_node']
            if molecular_rmsd(mols[representative],np.asarray(nodes[representative]['positions_A']),
                              mol,np.asarray(node['positions_A']))<tolerance_A:
                match=cluster;break
        if match is None:
            match=dict(id=node['id'],graph_smiles=node['graph_smiles'],physical_nodes=[],
                       representative_node=node['id'],energy_eV=node['energy_eV'])
            clusters.append(match)
        match['physical_nodes'].append(node['id'])
        if node['energy_eV']<match['energy_eV']:
            match['representative_node']=node['id'];match['energy_eV']=node['energy_eV']
    return clusters


def explore_shared(start,library,root_backend,strategy,outdir,potential_name,root,
                   protocol=SearchProtocol(),workers=2,threads_per_worker=2,
                   device='cpu',compile_model=False):
    """Explore one shared TransitionNet while independent workers evaluate seeds."""
    if workers<1 or threads_per_worker<1:
        raise ValueError('workers and threads_per_worker must be positive')
    if strategy not in ('geometry','center_random','bond_edits','arrows','hybrid'):
        raise ValueError('Unknown search strategy')
    if device!='cpu' and workers>1:
        raise ValueError('Shared multi-worker scheduling currently targets CPU execution')
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=False)
    numbers=start['atomic_numbers']
    validator=getattr(root_backend,'validate_system',None)
    if validator is not None:validator(numbers,start['charge'],start['multiplicity'])
    if start['multiplicity']!=1:
        raise ValueError('Current Lewis-graph registry supports closed-shell singlet exploration only')
    if strategy!='geometry' and (start['charge']!=0 or not set(numbers)<={1,6,7,8,15}):
        raise ValueError('Current symbolic proposal layer only supports neutral CHNOP')
    calculator=CountedCalculator(root_backend,protocol.total_evaluations)
    nodes=[];mols=[];edges=[];attempts=[];reservations=[]
    actions_used={};proposal_cache={};unsupported=set();exhausted=set();visits={}
    report=dict(start=start,strategy=strategy,protocol=asdict(protocol),nodes=nodes,edges=edges,
        attempts=attempts,reservations=reservations,evidence='MLIP_descents_not_DFT_IRC',
        library_policy=library.policy if library is not None else 'geometry_only',status='running',
        scheduler=dict(kind='central_species_registry_process_workers',workers=workers,
            threads_per_worker=threads_per_worker,dispatch='completion_driven',
            reservation_key='species_graph+conformer_cluster+action+geometry_variant',
            conformer_policy='symmetry-aware molecular RMSD clusters; lowest-energy representative per cluster',
            conformer_cluster_tolerance_A=.5),
        unsupported_conformer_clusters=[],exhausted_conformer_clusters=[])
    started=time.perf_counter();worker_evaluations=0;worker_model_seconds=0.

    def register(end,depth):
        x=np.asarray(end['positions_A']);mol=geometry_mol(numbers,x,start['charge'])
        for node,oldmol in zip(nodes,mols):
            if (abs(node['energy_eV']-end['energy_eV'])<=protocol.energy_tolerance_eV and
                molecular_rmsd(oldmol,np.asarray(node['positions_A']),mol,x)<protocol.geometry_tolerance_A):
                return node['id'],False
        idx=len(nodes);record=dict(id=idx,depth_discovered=depth,**end)
        record['graph_smiles']=graph_smiles(mol)
        nodes.append(record);mols.append(mol)
        return idx,True

    def update_projection():
        graph=nx.Graph();graph.add_nodes_from(range(len(nodes)));graph.add_edges_from(e['nodes'] for e in edges)
        roots=sorted(nx.node_connected_component(graph,0)) if nodes else []
        projection=project_species_network(nodes,edges,roots)
        report['species_nodes']=projection['nodes']
        report['physical_to_species']={str(k):v for k,v in projection['physical_to_species'].items()}
        report['conformer_clusters']=conformer_frontier(nodes,mols)
        report['root_component_nodes']=roots

    def save():
        update_projection()
        report['evaluations']=calculator.calls+worker_evaluations
        report['model_calls']=calculator.model_calls+sum(a.get('model_calls',0) for a in attempts)
        report['model_seconds']=calculator.model_seconds+worker_model_seconds
        report['elapsed_seconds']=time.perf_counter()-started
        report['unsupported_conformer_clusters']=sorted(unsupported)
        report['exhausted_conformer_clusters']=sorted(exhausted)
        atomic_json(outdir/'network.json',report)

    def select_task(reserved_evaluations):
        remaining=protocol.total_evaluations-(calculator.calls+worker_evaluations+reserved_evaluations)
        if len(attempts)>=protocol.max_attempts or remaining<=0:return None
        clusters=conformer_frontier(nodes,mols)
        candidates=[]
        for cluster in clusters:
            cluster_key=f"{cluster['graph_smiles']}|c{cluster['id']}"
            if cluster_key not in unsupported|exhausted:
                candidates.append((cluster['representative_node'],cluster_key))
        candidates.sort(key=lambda item:(visits.get(item[1],0),nodes[item[0]]['depth_discovered'],item[0]))
        known_graphs={node['graph_smiles'] for node in nodes}
        for node_id,cluster_key in candidates:
            if node_id not in proposal_cache:
                try:proposal_cache[node_id]=library.propose(mols[node_id]) if strategy!='geometry' else []
                except ValueError:proposal_cache[node_id]=[]
            proposals=proposal_cache[node_id]
            if strategy not in ('geometry','hybrid') and not proposals:
                unsupported.add(cluster_key);continue
            species_graph=nodes[node_id]['graph_smiles']
            used=actions_used.setdefault(cluster_key,{})
            choice=choose_action(proposals,used,known_graphs,strategy,protocol.geometry_seeds_per_node)
            if choice is None:
                exhausted.add(cluster_key);continue
            proposal,key,variant,seed_strategy=choice
            reservation_id,payload=reservation_identity(nodes[node_id],proposal,variant,seed_strategy,cluster_key)
            if any(item['id']==reservation_id for item in reservations):
                used[key]=max(used.get(key,0),variant+1);continue
            used[key]=used.get(key,0)+1;visits[cluster_key]=visits.get(cluster_key,0)+1
            identity=json.dumps(dict(graph=nodes[node_id]['graph_smiles'],
                positions=np.round(nodes[node_id]['positions_A'],6).tolist(),
                edits=proposal['edits'] if proposal else None,seed=protocol.random_seed,
                variant=variant),sort_keys=True)
            seed_rng=int.from_bytes(hashlib.sha256(identity.encode()).digest()[:4],'little')
            limit=min(protocol.evaluations_per_attempt,remaining)
            attempt_id=len(attempts)
            reservation=dict(id=reservation_id,identity=payload,status='running',attempt=attempt_id,
                             source_node=node_id,reserved_evaluations=limit)
            reservations.append(reservation)
            attempts.append(dict(id=attempt_id,source_node=node_id,status='running',
                reservation_id=reservation_id,proposal=proposal,evaluations=0,
                artifact=f'attempt_{attempt_id:03d}/result.json'))
            task=dict(numbers=numbers,positions_A=nodes[node_id]['positions_A'],charge=start['charge'],
                multiplicity=start['multiplicity'],proposal=proposal,variant=variant,
                seed_strategy=seed_strategy,random_seed=seed_rng,protocol=asdict(protocol),
                evaluation_limit=limit,destination=str((outdir/f'attempt_{attempt_id:03d}').resolve()))
            return task,reservation,attempt_id
        return None

    try:
        root_info,checks,initialization_evaluations=initialize_root(start,calculator,outdir,protocol)
        report['initial_relaxation_checks']=checks;report['initial_point']=root_info.copy()
        report['initialization_evaluations']=initialization_evaluations
        root_valid=(root_info['full_system_minimum'] if protocol.endpoint_acceptance=='full_system'
                    else root_info['carbon_skeleton_force_converged'])
        if not root_valid:
            report['status']='initial_minimum_unresolved';return report
        register(root_info,0);save()
        pending={};reserved_evaluations=0
        with ProcessPoolExecutor(max_workers=workers,initializer=_initialize_worker,
                initargs=(potential_name,str(Path(root).resolve()),device,compile_model,threads_per_worker)) as pool:
            while True:
                while len(pending)<workers:
                    selected=select_task(reserved_evaluations)
                    if selected is None:break
                    task,reservation,attempt_id=selected
                    future=pool.submit(_run_task,task);pending[future]=(reservation,attempt_id)
                    reserved_evaluations+=reservation['reserved_evaluations'];save()
                if not pending:break
                completed,_=wait(pending,return_when=FIRST_COMPLETED)
                for future in completed:
                    reservation,attempt_id=pending.pop(future)
                    reserved_evaluations-=reservation['reserved_evaluations']
                    attempt=attempts[attempt_id]
                    try:payload=future.result()
                    except Exception as exc:
                        payload=dict(result=dict(status='calculation_failed',evaluations=0,seconds=0.,
                            error=f'{type(exc).__name__}: {exc}',endpoints=[]),proposal=attempt['proposal'],
                            meta=attempt['proposal'])
                    result=payload['result'];attempt['proposal']=payload['meta']
                    attempt.update(status=result['status'],evaluations=result.get('evaluations',0),
                                   seconds=result.get('seconds',0.))
                    attempt['model_calls']=result.get('model_calls',0)
                    worker_evaluations+=attempt['evaluations'];worker_model_seconds+=result.get('model_seconds',0.)
                    reservation['status']='completed';reservation['actual_evaluations']=attempt['evaluations']
                    if result['status'] in ('validated_descents','validated_core_descents'):
                        source=nodes[attempt['source_node']]
                        registered=[register(endpoint,source['depth_discovered']+1) for endpoint in result['endpoints']]
                        ends=[item[0] for item in registered]
                        attempt['observed_nodes']=ends;attempt['new_nodes']=[i for i,new in registered if new]
                        attempt['source_is_endpoint']=attempt['source_node'] in ends
                        if ends[0]==ends[1]:attempt['status']='same_basin_return'
                        else:
                            duplicate=any(sorted(edge['nodes'])==sorted(ends) and
                                abs(edge['ts_energy_eV']-result['ts']['energy_eV'])<.03 and
                                aligned_rmsd(np.asarray(edge['ts_positions_A']),
                                             np.asarray(result['ts_positions_A']))<.15 for edge in edges)
                            attempt['status']='duplicate_connection' if duplicate else 'new_connection'
                            if not duplicate:
                                chemistry=classify_event(mols[ends[0]],mols[ends[1]])
                                edges.append(dict(id=len(edges),nodes=ends,attempt=attempt_id,
                                    proposed_from=attempt['source_node'],ts_energy_eV=result['ts']['energy_eV'],
                                    ts_positions_A=result['ts_positions_A'],
                                    barriers_eV=[e['barrier_eV'] for e in result['endpoints']],
                                    source_connected=attempt['source_node'] in ends,
                                    endpoint_acceptance=result['status'],
                                    endpoint_chemistry=chemistry,
                                    kind='conformational' if chemistry['resonance_equivalent'] else 'chemical'))
                    save()
                    print(json.dumps(dict(start=start['id'],strategy=strategy,attempt=attempt_id,
                        status=attempt['status'],nodes=len(nodes),species=len(report['species_nodes']),
                        edges=len(edges),evaluations=report['evaluations'])),flush=True)
        report['status']='completed'
        report['searched_nodes']=sorted({a['source_node'] for a in attempts})
        report['expanded_species']=sorted({nodes[i]['graph_smiles'] for i in report['searched_nodes']})
        report['conformer_cluster_batches']=dict(visits)
        report['stop_reason']=('evaluation_budget' if report['evaluations']>=protocol.total_evaluations else
            'attempt_budget' if len(attempts)>=protocol.max_attempts else 'frontier_exhausted')
    except BudgetExceeded:
        report['status']='initialization_budget_exhausted'
    finally:save()
    return report
