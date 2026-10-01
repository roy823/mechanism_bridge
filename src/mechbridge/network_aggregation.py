"""Merge compatible physical minimum/TS networks across strategies and orientations."""
import re
import numpy as np
from ase import Atoms

from .event_graph import geometry_mol,graph_smiles
from .reaction_network import molecular_rmsd,aligned_rmsd


def system_name(start_id):
    return re.sub(r'_[oc]\d+$','',start_id)


def aggregate_network_records(records,geometry_tolerance_A=.15,energy_tolerance_eV=.03,
                              ts_geometry_tolerance_A=.15,ts_energy_tolerance_eV=.03,group_by='system'):
    """Merge records on the same PES while preserving distinct minima and TS."""
    groups={}
    for record in records:
        network=record['network'];start=network['start']
        source_system=start.get('provenance',{}).get('aggregate_system',system_name(start['id']))
        if group_by=='inventory':name=Atoms(numbers=start['atomic_numbers']).get_chemical_formula(mode='hill')
        elif group_by=='system':name=source_system
        else:raise ValueError('group_by must be system or inventory')
        composition=tuple(sorted(start['atomic_numbers']));key=(name,composition,start['charge'],start['multiplicity'])
        group=groups.setdefault(key,dict(system=name,composition=composition,charge=start['charge'],
            multiplicity=start['multiplicity'],source_systems=set(),nodes=[],edges=[],root_nodes=set(),node_mols=[],node_map={},edge_map={}))
        group['source_systems'].add(source_system);numbers=start['atomic_numbers']
        local={}
        for node in network['nodes']:
            positions=np.asarray(node['positions_A']);mol=geometry_mol(numbers,positions,group['charge']);match=None
            for aggregate,oldmol in zip(group['nodes'],group['node_mols']):
                if (abs(aggregate['energy_eV']-node['energy_eV'])<=energy_tolerance_eV and
                    molecular_rmsd(oldmol,np.asarray(aggregate['positions_A']),mol,positions)<geometry_tolerance_A):
                    match=aggregate['id'];break
            if match is None:
                match=len(group['nodes']);group['nodes'].append(dict(id=match,atomic_numbers=numbers,positions_A=positions.tolist(),
                    energy_eV=node['energy_eV'],graph_smiles=graph_smiles(mol),origins=[]));group['node_mols'].append(mol)
            group['nodes'][match]['origins'].append(dict(record=record['id'],node=node['id']))
            local[node['id']]=match;group['node_map'][(record['id'],node['id'])]=match
        group['root_nodes'].add(local[0])
        for edge in network['edges']:
            nodes=[local[i] for i in edge['nodes']];ts=np.asarray(edge['ts_positions_A']);duplicate=None
            for aggregate in group['edges']:
                if (tuple(aggregate['atomic_numbers'])==tuple(numbers) and sorted(aggregate['nodes'])==sorted(nodes) and
                    abs(aggregate['ts_energy_eV']-edge['ts_energy_eV'])<ts_energy_tolerance_eV and
                    aligned_rmsd(np.asarray(aggregate['ts_positions_A']),ts)<ts_geometry_tolerance_A):
                    duplicate=aggregate;break
            origin=dict(record=record['id'],edge=edge['id'],attempt=edge['attempt'])
            if duplicate is None:
                duplicate=dict(id=len(group['edges']),nodes=nodes,atomic_numbers=numbers,ts_energy_eV=edge['ts_energy_eV'],
                    ts_positions_A=ts.tolist(),barriers_eV=edge['barriers_eV'],origins=[origin])
                group['edges'].append(duplicate)
            else:duplicate['origins'].append(origin)
            group['edge_map'][(record['id'],edge['id'])]=duplicate['id']
    for group in groups.values():
        group['root_nodes']=sorted(group['root_nodes']);group['source_systems']=sorted(group['source_systems']);group.pop('node_mols')
    return groups
