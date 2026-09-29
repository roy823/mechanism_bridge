"""Restricted unimolecular TST utilities. Never accept electronic barriers as free energies."""
import numpy as np
from scipy.constants import Boltzmann, Planck, Avogadro

def tst_unimolecular(delta_g_kcal_mol, temperature_K, kappa=1.):
    if temperature_K<=0 or kappa<=0 or not np.isfinite([delta_g_kcal_mol,temperature_K,kappa]).all():
        raise ValueError('Finite values, positive temperature and transmission factor required')
    return kappa*Boltzmann*temperature_K/Planck*np.exp(-delta_g_kcal_mol*4184/(Avogadro*Boltzmann*temperature_K))

def rate_matrix(states, edges, temperature_K):
    """Row-generator Q for unimolecular/reversibly lumped states with explicit G barriers.

    Not a solution-reaction mass-action solver; no concentration or standard-state guessing.
    """
    indices={name:i for i,name in enumerate(states)}
    if len(indices)!=len(states):raise ValueError('State names must be unique')
    q=np.zeros((len(states),len(states)))
    for edge in edges:
        if edge.get('barrier_kind')!='activation_gibbs_free_energy' or edge.get('molecularity')!=1:
            raise ValueError('Only explicitly labelled unimolecular Gibbs barriers supported')
        if not edge.get('connectivity_verified',False):raise ValueError('Connectivity evidence required')
        i,j=indices[edge['from']],indices[edge['to']]
        if i==j:raise ValueError('Self-events do not enter the state generator')
        q[i,j]+=tst_unimolecular(edge['delta_g_kcal_mol'],temperature_K,edge.get('kappa',1.))
    q[np.diag_indices_from(q)]=-q.sum(axis=1)
    return q
