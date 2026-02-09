#!/usr/bin/env python3
# ============================================================================
# boltz_extract_metrics.py
# Extract interface quality metrics from Boltz-1 predictions
# Based on colabfold_extract_metrics.py with added PDE-based metrics
# ============================================================================

import argparse
import sys
import csv
from pathlib import Path

import json
import math
import numpy as np

from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

def parse_pdb_structure(structure_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', str(structure_file))
    model = list(structure.get_models())[0]
    chains = []
    chain_residues = {}
    for chain in model:
        chains.append(chain.id)
        residues = [res for res in chain if is_aa(res, standard=True)]
        chain_residues[chain.id] = residues
    #
    return structure, chains, chain_residues

def get_interface(chain_residues, chain_a, chain_b, interface_threshold):
    #
    residues_a = chain_residues[chain_a]
    residues_b = chain_residues[chain_b]
    #
    atoms_a = []
    atom_to_residue_a = {}
    #
    for i, res in enumerate(residues_a):
        for atom in res:
            if atom.element != 'H':
                atoms_a.append(atom)
                atom_to_residue_a[atom] = i
    #
    atoms_b = []
    atom_to_residue_b = {}
    #
    for i, res in enumerate(residues_b):
        for atom in res:
            if atom.element != 'H':
                atoms_b.append(atom)
                atom_to_residue_b[atom] = i
    #
    ns = NeighborSearch(atoms_a + atoms_b)
    #
    interface_a = set()
    interface_b = set()
    contact_pairs = set()
    #
    for atom_a in atoms_a:
        close_atoms = ns.search(atom_a.coord, interface_threshold, 'A')
        for atom_b in close_atoms:
            if atom_b in atoms_b:
                interface_a.add(atom_to_residue_a[atom_a])
                interface_b.add(atom_to_residue_b[atom_b])
                contact_pairs.add((atom_to_residue_a[atom_a], atom_to_residue_b[atom_b]))
    
    return interface_a, interface_b, contact_pairs

def get_chain_mapping(chains, chain_residues):
    #
    mapping = {}
    current_idx = 0
    #
    for chain_id in chains:
        n_res = len(chain_residues[chain_id])
        mapping[chain_id] = (current_idx, current_idx + n_res)
        current_idx += n_res
    #
    return mapping

def get_cb_positions(chain_residues, chain_id):
    """
    Extract CB positions (CA for GLY) and pLDDT values
    """
    residues = chain_residues[chain_id]
    cb_coords = []
    plddt_values = []
    #
    for res in residues:
        # Get CB atom, or CA for glycine
        if res.get_resname() == 'GLY':
            if 'CA' in res:
                atom = res['CA']
            else:
                continue
        else:
            if 'CB' in res:
                atom = res['CB']
            elif 'CA' in res:
                atom = res['CA']
            else:
                continue
        #
        cb_coords.append(atom.get_coord())
        plddt_values.append(atom.get_bfactor())
    #
    return np.array(cb_coords), np.array(plddt_values)


def compute_cb_distances(cb_coords_a, cb_coords_b):
    """
    Calculate pairwise Euclidean distances between CB atoms
    """
    # Use broadcasting: (n, 1, 3) - (1, m, 3) = (n, m, 3)
    diff = cb_coords_a[:, np.newaxis, :] - cb_coords_b[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=2))
    return distances


def calc_d0(n_residues):
    """
    Calculate d0 normalization factor (Yang & Skolnick 2004)
    """
    n_residues = float(n_residues)
    if n_residues > 26:
        d0 = 1.24 * (n_residues - 15) ** (1.0/3.0) - 1.8
    else:
        d0 = 1.0
    return max(1.0, d0)


def ptm_transform(pae_values, d0):
    """
    Apply PTM transformation: 1/(1 + (pae/d0)^2)
    """
    return 1.0 / (1.0 + (pae_values / d0) ** 2.0)

# ============================================================================
# Compute ipLDDT, interchain_PAE, interface_PAE metrics
# ============================================================================

def compute_interface_metrics(plddt_array, pae_matrix,
                              interface_a, interface_b,
                              contact_pairs,
                              chain_a_start, chain_b_start,
                              chain_a_end, chain_b_end, pae_threshold):
    metrics = {}
    #
    # ipLDDT
    if len(interface_a) > 0 or len(interface_b) > 0:
        interface_indices = []
        for idx in interface_a:
            interface_indices.append(chain_a_start + idx)
        for idx in interface_b:
            interface_indices.append(chain_b_start + idx)
        #
        if interface_indices:
            metrics['ipLDDT_mean'] = float(np.mean(plddt_array[interface_indices]))
            metrics['ipLDDT_min'] = float(np.min(plddt_array[interface_indices]))
            metrics['ipLDDT_max'] = float(np.max(plddt_array[interface_indices]))
        else:
            metrics['ipLDDT_mean'] = np.nan
            metrics['ipLDDT_min'] = np.nan
            metrics['ipLDDT_max'] = np.nan
    else:
        metrics['ipLDDT_mean'] = np.nan
        metrics['ipLDDT_min'] = np.nan
        metrics['ipLDDT_max'] = np.nan
    #
    # interchain PAE
    chain_a_indices = list(range(chain_a_start, chain_a_end))
    chain_b_indices = list(range(chain_b_start, chain_b_end))
    #
    if len(chain_a_indices) > 0 and len(chain_b_indices) > 0:
        pae_ab = pae_matrix[np.ix_(chain_a_indices, chain_b_indices)]
        pae_ba = pae_matrix[np.ix_(chain_b_indices, chain_a_indices)]
        all_cross_pae = np.concatenate([pae_ab.flatten(), pae_ba.flatten()])
        metrics['interchain_pAE_mean'] = float(np.mean(all_cross_pae))
        metrics['interchain_pAE_min'] = float(np.min(all_cross_pae))
        metrics['interchain_pAE_max'] = float(np.max(all_cross_pae))
    else:
        metrics['interchain_pAE_mean'] = np.nan
        metrics['interchain_pAE_min'] = np.nan
        metrics['interchain_pAE_max'] = np.nan
    #
    # interface_iPAE
    if len(contact_pairs) > 0:
        pae_values = []
        for (idx_a, idx_b) in contact_pairs:
            matrix_idx_a = chain_a_start + idx_a
            matrix_idx_b = chain_b_start + idx_b
            pae_values.append(pae_matrix[matrix_idx_a, matrix_idx_b])
            pae_values.append(pae_matrix[matrix_idx_b, matrix_idx_a])
        #
        pae_values = np.array(pae_values)
        metrics['ipAE_mean'] = float(np.mean(pae_values))
        metrics['ipAE_min'] = float(np.min(pae_values))
        metrics['ipAE_max'] = float(np.max(pae_values))
        metrics['ipAE_pct_lt_threshold'] = float(np.mean(pae_values < pae_threshold) * 100)
    else:
        metrics['ipAE_mean'] = np.nan
        metrics['ipAE_min'] = np.nan
        metrics['ipAE_max'] = np.nan
        metrics['ipAE_pct_lt_threshold'] = np.nan
    #
    return metrics


# ============================================================================
# Compute ipSAE, pDockQ, pDockQ2, LIS
# ============================================================================

def compute_ipsae_variants(pae_matrix, chain_a_start, chain_a_end, 
                          chain_b_start, chain_b_end, pae_cutoff):
    """
    Compute all ipSAE variants (d0res, d0dom, d0chn) with max and min aggregations
    """
    metrics = {}
    #
    n_res_a = chain_a_end - chain_a_start
    n_res_b = chain_b_end - chain_b_start
    #
    # Extract chain indices
    chain_a_indices = np.arange(chain_a_start, chain_a_end)
    chain_b_indices = np.arange(chain_b_start, chain_b_end)
    #
    # --- Variant 1: ipSAE_d0chn (chain-length normalization) ---
    n0chn = n_res_a + n_res_b
    d0chn = calc_d0(n0chn)
    ptm_matrix_d0chn = ptm_transform(pae_matrix, d0chn)
    #
    # A -> B direction
    scores_AB_d0chn = []
    for i in chain_a_indices:
        valid_mask = pae_matrix[i, chain_b_indices] < pae_cutoff
        if np.any(valid_mask):
            score = np.mean(ptm_matrix_d0chn[i, chain_b_indices[valid_mask]])
            scores_AB_d0chn.append(score)
    #
    if len(scores_AB_d0chn) > 0:
        metrics['ipSAE_d0chn_AB_max'] = float(np.max(scores_AB_d0chn))
        metrics['ipSAE_d0chn_AB_min'] = float(np.min(scores_AB_d0chn))
    else:
        metrics['ipSAE_d0chn_AB_max'] = 0.0
        metrics['ipSAE_d0chn_AB_min'] = 0.0
    #
    # B -> A direction
    scores_BA_d0chn = []
    for i in chain_b_indices:
        valid_mask = pae_matrix[i, chain_a_indices] < pae_cutoff
        if np.any(valid_mask):
            score = np.mean(ptm_matrix_d0chn[i, chain_a_indices[valid_mask]])
            scores_BA_d0chn.append(score)
    #
    if len(scores_BA_d0chn) > 0:
        metrics['ipSAE_d0chn_BA_max'] = float(np.max(scores_BA_d0chn))
        metrics['ipSAE_d0chn_BA_min'] = float(np.min(scores_BA_d0chn))
    else:
        metrics['ipSAE_d0chn_BA_max'] = 0.0
        metrics['ipSAE_d0chn_BA_min'] = 0.0
    #
    # --- Variant 2: ipSAE_d0dom (domain-size normalization) ---
    # Find interface residues for domain size calculation
    interface_residues_a = set()
    interface_residues_b = set()
    #
    for i in chain_a_indices:
        valid_mask = pae_matrix[i, chain_b_indices] < pae_cutoff
        if np.any(valid_mask):
            interface_residues_a.add(i)
            for j in chain_b_indices[valid_mask]:
                interface_residues_b.add(j)
    #
    for i in chain_b_indices:
        valid_mask = pae_matrix[i, chain_a_indices] < pae_cutoff
        if np.any(valid_mask):
            interface_residues_b.add(i)
            for j in chain_a_indices[valid_mask]:
                interface_residues_a.add(j)
    #
    n0dom = len(interface_residues_a) + len(interface_residues_b)
    if n0dom > 0:
        d0dom = calc_d0(n0dom)
        ptm_matrix_d0dom = ptm_transform(pae_matrix, d0dom)
        #
        # A -> B direction
        scores_AB_d0dom = []
        for i in chain_a_indices:
            valid_mask = pae_matrix[i, chain_b_indices] < pae_cutoff
            if np.any(valid_mask):
                score = np.mean(ptm_matrix_d0dom[i, chain_b_indices[valid_mask]])
                scores_AB_d0dom.append(score)
        #
        if len(scores_AB_d0dom) > 0:
            metrics['ipSAE_d0dom_AB_max'] = float(np.max(scores_AB_d0dom))
            metrics['ipSAE_d0dom_AB_min'] = float(np.min(scores_AB_d0dom))
        #
        else:
            metrics['ipSAE_d0dom_AB_max'] = 0.0
            metrics['ipSAE_d0dom_AB_min'] = 0.0
        #
        # B -> A direction
        scores_BA_d0dom = []
        for i in chain_b_indices:
            valid_mask = pae_matrix[i, chain_a_indices] < pae_cutoff
            if np.any(valid_mask):
                score = np.mean(ptm_matrix_d0dom[i, chain_a_indices[valid_mask]])
                scores_BA_d0dom.append(score)
        #
        if len(scores_BA_d0dom) > 0:
            metrics['ipSAE_d0dom_BA_max'] = float(np.max(scores_BA_d0dom))
            metrics['ipSAE_d0dom_BA_min'] = float(np.min(scores_BA_d0dom))
        else:
            metrics['ipSAE_d0dom_BA_max'] = 0.0
            metrics['ipSAE_d0dom_BA_min'] = 0.0
    else:
        metrics['ipSAE_d0dom_AB_max'] = 0.0
        metrics['ipSAE_d0dom_AB_min'] = 0.0
        metrics['ipSAE_d0dom_BA_max'] = 0.0
        metrics['ipSAE_d0dom_BA_min'] = 0.0
    #
    # --- Variant 3: ipSAE_d0res (per-residue normalization) ---
    # A -> B direction
    scores_AB_d0res = []
    for i in chain_a_indices:
        valid_mask = pae_matrix[i, chain_b_indices] < pae_cutoff
        if np.any(valid_mask):
            n0res = np.sum(valid_mask)
            d0res = calc_d0(n0res)
            pae_vals = pae_matrix[i, chain_b_indices[valid_mask]]
            ptm_vals = ptm_transform(pae_vals, d0res)
            score = np.mean(ptm_vals)
            scores_AB_d0res.append(score)
    #
    if len(scores_AB_d0res) > 0:
        metrics['ipSAE_d0res_AB_max'] = float(np.max(scores_AB_d0res))
        metrics['ipSAE_d0res_AB_min'] = float(np.min(scores_AB_d0res))
    else:
        metrics['ipSAE_d0res_AB_max'] = 0.0
        metrics['ipSAE_d0res_AB_min'] = 0.0
    #
    # B -> A direction
    scores_BA_d0res = []
    for i in chain_b_indices:
        valid_mask = pae_matrix[i, chain_a_indices] < pae_cutoff
        if np.any(valid_mask):
            n0res = np.sum(valid_mask)
            d0res = calc_d0(n0res)
            pae_vals = pae_matrix[i, chain_a_indices[valid_mask]]
            ptm_vals = ptm_transform(pae_vals, d0res)
            score = np.mean(ptm_vals)
            scores_BA_d0res.append(score)
    #
    if len(scores_BA_d0res) > 0:
        metrics['ipSAE_d0res_BA_max'] = float(np.max(scores_BA_d0res))
        metrics['ipSAE_d0res_BA_min'] = float(np.min(scores_BA_d0res))
    else:
        metrics['ipSAE_d0res_BA_max'] = 0.0
        metrics['ipSAE_d0res_BA_min'] = 0.0
    #
    return metrics


def compute_pdockqs(cb_coords_a, cb_coords_b, plddt_a, plddt_b, 
                   pae_matrix, chain_a_start, chain_b_start, dockq_threshold):
    """
    Compute pDockQ score
    """
    # Compute distance matrix
    distances = compute_cb_distances(cb_coords_a, cb_coords_b)
    # Find contacts
    contacts = distances <= dockq_threshold
    #
    if not np.any(contacts):
        return 0.0, 0.0, 0.0
    #
    # Find unique interface residues
    interface_residues_a = set(np.where(np.any(contacts, axis=1))[0])
    interface_residues_b = set(np.where(np.any(contacts, axis=0))[0])
    #
    interface_residues = list(interface_residues_a) + [i + len(plddt_a) for i in interface_residues_b]
    #
    # Combine pLDDT arrays
    combined_plddt = np.concatenate([plddt_a, plddt_b])
    #
    mean_plddt = np.mean(combined_plddt[interface_residues])
    #
    npairs = np.sum(contacts)    
    if npairs > 0:
        x = mean_plddt * math.log10(npairs)
        pdockq = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
    else:
        pdockq = 0.0
    #
    # Calculate mean PTM for contacts
    pae_values_AB = []
    pae_values_BA = []
    #
    contact_indices = np.where(contacts)
    #
    for i, j in zip(contact_indices[0], contact_indices[1]):
        matrix_idx_a = chain_a_start + i
        matrix_idx_b = chain_b_start + j
        pae_values_AB.append(pae_matrix[matrix_idx_a, matrix_idx_b])
        pae_values_BA.append(pae_matrix[matrix_idx_b, matrix_idx_a])
    #
    if len(pae_values_AB) > 0:
        pae_values_AB = np.array(pae_values_AB)
        ptm_values_AB = ptm_transform(pae_values_AB, 10.0)
        mean_ptm_AB = np.mean(ptm_values_AB)
        x_AB = mean_plddt * mean_ptm_AB
        pdockq2_AB = 1.31 / (1 + math.exp(-0.075 * (x_AB - 84.733))) + 0.005
    else:
        pdockq2_AB = 0.0
    #
    if len(pae_values_BA) > 0:
        pae_values_BA = np.array(pae_values_BA)
        ptm_values_BA = ptm_transform(pae_values_BA, 10.0)
        mean_ptm_BA = np.mean(ptm_values_BA)
        x_BA = mean_plddt * mean_ptm_BA
        pdockq2_BA = 1.31 / (1 + math.exp(-0.075 * (x_BA - 84.733))) + 0.005
    else:
        pdockq2_BA = 0.0
    #
    return pdockq, pdockq2_AB, pdockq2_BA

def compute_lis(pae_matrix, chain_a_start, chain_a_end, 
               chain_b_start, chain_b_end, pae_threshold):
    """
    Compute LIS (Local Interaction Score)
    """
    # A -> B direction
    pae_ab = pae_matrix[chain_a_start:chain_a_end, chain_b_start:chain_b_end]
    valid_pae_ab = pae_ab[pae_ab < pae_threshold]
    
    if len(valid_pae_ab) > 0:
        scores_ab = (pae_threshold - valid_pae_ab) / pae_threshold
        lis_ab = float(np.mean(scores_ab))
    else:
        lis_ab = 0.0
    
    # B -> A direction
    pae_ba = pae_matrix[chain_b_start:chain_b_end, chain_a_start:chain_a_end]
    valid_pae_ba = pae_ba[pae_ba < pae_threshold]
    
    if len(valid_pae_ba) > 0:
        scores_ba = (pae_threshold - valid_pae_ba) / pae_threshold
        lis_ba = float(np.mean(scores_ba))
    else:
        lis_ba = 0.0
    
    lis_avg = (lis_ab + lis_ba) / 2.0
    
    return lis_ab, lis_ba, lis_avg

# ============================================================================
# NEW: PDE-based metrics (analogous to PAE-based metrics)
# ============================================================================

def compute_interface_metrics_pde(pde_matrix, interface_a, interface_b,
                                  contact_pairs,
                                  chain_a_start, chain_b_start,
                                  chain_a_end, chain_b_end):
    """
    Compute PDE-based interface metrics (analogous to PAE metrics)
    
    Returns:
    - interchain_pDE_mean/min/max
    - ipDE_mean/min/max
    """
    metrics = {}
    #
    # interchain PDE
    chain_a_indices = list(range(chain_a_start, chain_a_end))
    chain_b_indices = list(range(chain_b_start, chain_b_end))
    #
    if len(chain_a_indices) > 0 and len(chain_b_indices) > 0:
        pde_ab = pde_matrix[np.ix_(chain_a_indices, chain_b_indices)]
        pde_ba = pde_matrix[np.ix_(chain_b_indices, chain_a_indices)]
        all_cross_pde = np.concatenate([pde_ab.flatten(), pde_ba.flatten()])
        metrics['interchain_pDE_mean'] = float(np.mean(all_cross_pde))
        metrics['interchain_pDE_min'] = float(np.min(all_cross_pde))
        metrics['interchain_pDE_max'] = float(np.max(all_cross_pde))
    else:
        metrics['interchain_pDE_mean'] = np.nan
        metrics['interchain_pDE_min'] = np.nan
        metrics['interchain_pDE_max'] = np.nan
    #
    # interface_ipDE
    if len(contact_pairs) > 0:
        pde_values = []
        for (idx_a, idx_b) in contact_pairs:
            matrix_idx_a = chain_a_start + idx_a
            matrix_idx_b = chain_b_start + idx_b
            pde_values.append(pde_matrix[matrix_idx_a, matrix_idx_b])
            pde_values.append(pde_matrix[matrix_idx_b, matrix_idx_a])
        #
        pde_values = np.array(pde_values)
        metrics['ipDE_mean'] = float(np.mean(pde_values))
        metrics['ipDE_min'] = float(np.min(pde_values))
        metrics['ipDE_max'] = float(np.max(pde_values))
    else:
        metrics['ipDE_mean'] = np.nan
        metrics['ipDE_min'] = np.nan
        metrics['ipDE_max'] = np.nan
    #
    return metrics


def compute_ipsde_variants(pde_matrix, chain_a_start, chain_a_end,
                           chain_b_start, chain_b_end, pde_cutoff):
    """
    Compute ipSDE (Interface Predicted Structural Distance Error)
    Analogous to ipSAE but using PDE matrix
    Using d0res normalization (most sensitive)
    
    Returns:
    - ipSDE_d0res_AB_max/min
    - ipSDE_d0res_BA_max/min
    """
    metrics = {}
    #
    chain_a_indices = np.arange(chain_a_start, chain_a_end)
    chain_b_indices = np.arange(chain_b_start, chain_b_end)
    #
    # A -> B direction
    scores_AB = []
    for i in chain_a_indices:
        valid_mask = pde_matrix[i, chain_b_indices] < pde_cutoff
        if np.any(valid_mask):
            n0res = np.sum(valid_mask)
            d0res = calc_d0(n0res)
            pde_vals = pde_matrix[i, chain_b_indices[valid_mask]]
            ptm_vals = ptm_transform(pde_vals, d0res)
            score = np.mean(ptm_vals)
            scores_AB.append(score)
    #
    if len(scores_AB) > 0:
        metrics['ipSDE_d0res_AB_max'] = float(np.max(scores_AB))
        metrics['ipSDE_d0res_AB_min'] = float(np.min(scores_AB))
    else:
        metrics['ipSDE_d0res_AB_max'] = 0.0
        metrics['ipSDE_d0res_AB_min'] = 0.0
    #
    # B -> A direction
    scores_BA = []
    for i in chain_b_indices:
        valid_mask = pde_matrix[i, chain_a_indices] < pde_cutoff
        if np.any(valid_mask):
            n0res = np.sum(valid_mask)
            d0res = calc_d0(n0res)
            pde_vals = pde_matrix[i, chain_a_indices[valid_mask]]
            ptm_vals = ptm_transform(pde_vals, d0res)
            score = np.mean(ptm_vals)
            scores_BA.append(score)
    #
    if len(scores_BA) > 0:
        metrics['ipSDE_d0res_BA_max'] = float(np.max(scores_BA))
        metrics['ipSDE_d0res_BA_min'] = float(np.min(scores_BA))
    else:
        metrics['ipSDE_d0res_BA_max'] = 0.0
        metrics['ipSDE_d0res_BA_min'] = 0.0
    #
    return metrics


def compute_lis_pde(pde_matrix, chain_a_start, chain_a_end,
                    chain_b_start, chain_b_end, pde_threshold):
    """
    Compute LIS using PDE instead of PAE
    PDE-LIS = (threshold - PDE) / threshold for PDE < threshold
    
    Returns:
    - lis_pde_AB, lis_pde_BA, lis_pde_avg
    """
    # A -> B direction
    pde_ab = pde_matrix[chain_a_start:chain_a_end, chain_b_start:chain_b_end]
    valid_pde_ab = pde_ab[pde_ab < pde_threshold]
    
    if len(valid_pde_ab) > 0:
        scores_ab = (pde_threshold - valid_pde_ab) / pde_threshold
        lis_pde_ab = float(np.mean(scores_ab))
    else:
        lis_pde_ab = 0.0
    
    # B -> A direction
    pde_ba = pde_matrix[chain_b_start:chain_b_end, chain_a_start:chain_a_end]
    valid_pde_ba = pde_ba[pde_ba < pde_threshold]
    
    if len(valid_pde_ba) > 0:
        scores_ba = (pde_threshold - valid_pde_ba) / pde_threshold
        lis_pde_ba = float(np.mean(scores_ba))
    else:
        lis_pde_ba = 0.0
    
    lis_pde_avg = (lis_pde_ab + lis_pde_ba) / 2.0
    
    return lis_pde_ab, lis_pde_ba, lis_pde_avg


def compute_pdockq2_pde(cb_coords_a, cb_coords_b, plddt_a, plddt_b,
                        pde_matrix, chain_a_start, chain_b_start, distance_cutoff):
    """
    Compute pDockQ2 using PDE instead of PAE
    Same formula as pDockQ2 but with PDE at contact points
    
    Returns:
    - pDockQ2_pde_AB, pDockQ2_pde_BA
    """
    # Compute distance matrix
    distances = compute_cb_distances(cb_coords_a, cb_coords_b)
    contacts = distances <= distance_cutoff
    #
    if not np.any(contacts):
        return 0.0, 0.0
    #
    # Find unique interface residues
    interface_residues_a = set(np.where(np.any(contacts, axis=1))[0])
    interface_residues_b = set(np.where(np.any(contacts, axis=0))[0])
    interface_residues = list(interface_residues_a) + [i + len(plddt_a) for i in interface_residues_b]
    #
    # Combine pLDDT arrays
    combined_plddt = np.concatenate([plddt_a, plddt_b])
    mean_plddt = np.mean(combined_plddt[interface_residues])
    #
    # Calculate mean PTM for contacts using PDE
    pde_values_AB = []
    pde_values_BA = []
    #
    contact_indices = np.where(contacts)
    #
    for i, j in zip(contact_indices[0], contact_indices[1]):
        matrix_idx_a = chain_a_start + i
        matrix_idx_b = chain_b_start + j
        pde_values_AB.append(pde_matrix[matrix_idx_a, matrix_idx_b])
        pde_values_BA.append(pde_matrix[matrix_idx_b, matrix_idx_a])
    #
    if len(pde_values_AB) > 0:
        pde_values_AB = np.array(pde_values_AB)
        ptm_values_AB = ptm_transform(pde_values_AB, 10.0)
        mean_ptm_AB = np.mean(ptm_values_AB)
        x_AB = mean_plddt * mean_ptm_AB
        pdockq2_pde_AB = 1.31 / (1 + math.exp(-0.075 * (x_AB - 84.733))) + 0.005
    else:
        pdockq2_pde_AB = 0.0
    #
    if len(pde_values_BA) > 0:
        pde_values_BA = np.array(pde_values_BA)
        ptm_values_BA = ptm_transform(pde_values_BA, 10.0)
        mean_ptm_BA = np.mean(ptm_values_BA)
        x_BA = mean_plddt * mean_ptm_BA
        pdockq2_pde_BA = 1.31 / (1 + math.exp(-0.075 * (x_BA - 84.733))) + 0.005
    else:
        pdockq2_pde_BA = 0.0
    #
    return float(pdockq2_pde_AB), float(pdockq2_pde_BA)


# ============================================================================
# Main extraction function
# ============================================================================

def extract_sample_metrics(sample_dir, args):
    """
    Extract all metrics for a single Boltz-1 sample directory
    """
    #
    sample_dir = Path(sample_dir)
    sample_name = sample_dir.name
    #
    # Step 1: Find files (MODIFIED for Boltz format)
    # Find confidence JSON
    json_files = list(sample_dir.glob("*.json"))
    if not json_files:
        return None
    json_file = json_files[0]
    #
    # Find PDB
    pdb_files = list(sample_dir.glob("*_model_*.pdb"))
    if not pdb_files:
        return None
    pdb_file = pdb_files[0]
    #
    # Find pLDDT NPZ
    plddt_files = list(sample_dir.glob("plddt_*_model_*.npz"))
    if not plddt_files:
        return None
    plddt_file = plddt_files[0]
    #
    # Find PAE NPZ
    pae_files = list(sample_dir.glob("pae_*_model_*.npz"))
    if not pae_files:
        return None
    pae_file = pae_files[0]
    #
    # Find PDE NPZ (optional)
    pde_files = list(sample_dir.glob("pde_*_model_*.npz"))
    has_pde = len(pde_files) > 0
    pde_file = pde_files[0] if has_pde else None
    #
    # Step 2: Load data (MODIFIED for Boltz format)
    # Load confidence JSON
    with open(json_file) as f:
        confidence = json.load(f)
    #
    # Basic metrics from JSON
    metrics = {
        'sample_name': sample_name,
        'confidence_score': confidence.get('confidence_score'),
        'ptm': confidence.get('ptm'),
        'iptm': confidence.get('iptm'),
        'ligand_iptm': confidence.get('ligand_iptm'),
        'protein_iptm': confidence.get('protein_iptm'),
        'complex_plddt': confidence.get('complex_plddt'),
        'complex_iplddt': confidence.get('complex_iplddt'),
        'complex_pde': confidence.get('complex_pde'),
        'complex_ipde': confidence.get('complex_ipde'),
    }
    #
    # Load arrays from NPZ files
    plddt_data = np.load(plddt_file)
    plddt_array = plddt_data['plddt']
    #
    pae_data = np.load(pae_file)
    pae_matrix = pae_data['pae']
    #
    if has_pde:
        pde_data = np.load(pde_file)
        pde_matrix = pde_data['pde']
    else:
        pde_matrix = None
    #
    metrics['pLDDT_mean'] = float(np.mean(plddt_array))
    metrics['pAE_mean'] = float(np.mean(pae_matrix))
    if has_pde:
        metrics['pDE_mean'] = float(np.mean(pde_matrix))
    else:
        metrics['pDE_mean'] = np.nan
    #
    # Parse PDB structure
    structure, chains, chain_residues = parse_pdb_structure(pdb_file)
    #
    # Extract per-chain metrics from Boltz JSON (NEW)
    if len(chains) >= 2:
        chains_ptm = confidence.get('chains_ptm', {})
        metrics['chains_ptm_A'] = chains_ptm.get('0', np.nan)
        metrics['chains_ptm_B'] = chains_ptm.get('1', np.nan)
        #
        pair_chains_iptm = confidence.get('pair_chains_iptm', {})
        metrics['pair_chains_iptm_AA'] = pair_chains_iptm.get('0', {}).get('0', np.nan)
        metrics['pair_chains_iptm_AB'] = pair_chains_iptm.get('0', {}).get('1', np.nan)
        metrics['pair_chains_iptm_BA'] = pair_chains_iptm.get('1', {}).get('0', np.nan)
        metrics['pair_chains_iptm_BB'] = pair_chains_iptm.get('1', {}).get('1', np.nan)
    else:
        metrics['chains_ptm_A'] = np.nan
        metrics['chains_ptm_B'] = np.nan
        metrics['pair_chains_iptm_AA'] = np.nan
        metrics['pair_chains_iptm_AB'] = np.nan
        metrics['pair_chains_iptm_BA'] = np.nan
        metrics['pair_chains_iptm_BB'] = np.nan
    #
    # Return nan if no interface
    if len(chains) >= 2:
        #
        # Get chain mapping
        chain_a, chain_b = chains[0], chains[1]
        #
        chain_mapping = get_chain_mapping(chains, chain_residues)
        chain_a_start, chain_a_end = chain_mapping[chain_a]
        chain_b_start, chain_b_end = chain_mapping[chain_b]
        #
        # Find interface
        interface_a, interface_b, contact_pairs = get_interface(
            chain_residues, chain_a, chain_b, args.interface_threshold
        )
        metrics['num_interface_residues_A'] = len(interface_a)
        metrics['num_interface_residues_B'] = len(interface_b)
        #
        # Compute interface metrics (PAE-based - SAME as ColabFold)
        interface_metrics = compute_interface_metrics(
            plddt_array, pae_matrix,
            interface_a, interface_b,
            contact_pairs,
            chain_a_start, chain_b_start,
            chain_a_end, chain_b_end,
            args.pae_threshold
        )
        metrics.update(interface_metrics)
        #
        # Compute ipSAE variants (SAME as ColabFold)
        ipsae_metrics = compute_ipsae_variants(
            pae_matrix, chain_a_start, chain_a_end,
            chain_b_start, chain_b_end, args.pae_threshold
        )
        metrics.update(ipsae_metrics)
        #
        # Extract CB positions
        cb_coords_a, plddt_a = get_cb_positions(chain_residues, chain_a)
        cb_coords_b, plddt_b = get_cb_positions(chain_residues, chain_b)
        #
        # Compute pDockQ/pDockQ2 bidirectional (SAME as ColabFold)
        pdockq, pdockq2_AB, pdockq2_BA = compute_pdockqs(
            cb_coords_a, cb_coords_b, 
            plddt_a, plddt_b, pae_matrix, 
            chain_a_start, chain_b_start, args.pdockq_threshold
        )
        metrics['pDockQ'] = pdockq
        metrics['pDockQ2_AB'] = pdockq2_AB
        metrics['pDockQ2_BA'] = pdockq2_BA
        #
        # Compute LIS (SAME as ColabFold)
        lis_AB, lis_BA, lis_avg = compute_lis(
            pae_matrix, chain_a_start, chain_a_end, 
            chain_b_start, chain_b_end, args.pae_threshold
        )
        metrics['lis_AB'] = lis_AB
        metrics['lis_BA'] = lis_BA
        metrics['lis_avg'] = lis_avg
        #
        # Compute PDE-based metrics (NEW)
        if has_pde:
            # Interface metrics (PDE)
            pde_interface_metrics = compute_interface_metrics_pde(
                pde_matrix, interface_a, interface_b,
                contact_pairs,
                chain_a_start, chain_b_start,
                chain_a_end, chain_b_end
            )
            metrics.update(pde_interface_metrics)
            #
            # ipSDE variants
            ipsde_metrics = compute_ipsde_variants(
                pde_matrix, chain_a_start, chain_a_end,
                chain_b_start, chain_b_end, args.pde_threshold
            )
            metrics.update(ipsde_metrics)
            #
            # LIS-PDE
            lis_pde_AB, lis_pde_BA, lis_pde_avg = compute_lis_pde(
                pde_matrix, chain_a_start, chain_a_end,
                chain_b_start, chain_b_end, args.pde_threshold
            )
            metrics['lis_pde_AB'] = lis_pde_AB
            metrics['lis_pde_BA'] = lis_pde_BA
            metrics['lis_pde_avg'] = lis_pde_avg
            #
            # pDockQ2-PDE
            if len(cb_coords_a) > 0 and len(cb_coords_b) > 0:
                pdockq2_pde_AB, pdockq2_pde_BA = compute_pdockq2_pde(
                    cb_coords_a, cb_coords_b,
                    plddt_a, plddt_b, pde_matrix,
                    chain_a_start, chain_b_start, args.pdockq_threshold
                )
                metrics['pDockQ2_pde_AB'] = pdockq2_pde_AB
                metrics['pDockQ2_pde_BA'] = pdockq2_pde_BA
            else:
                metrics['pDockQ2_pde_AB'] = np.nan
                metrics['pDockQ2_pde_BA'] = np.nan
        else:
            # No PDE file - set all PDE metrics to nan
            metrics['interchain_pDE_mean'] = np.nan
            metrics['interchain_pDE_min'] = np.nan
            metrics['interchain_pDE_max'] = np.nan
            metrics['ipDE_mean'] = np.nan
            metrics['ipDE_min'] = np.nan
            metrics['ipDE_max'] = np.nan
            metrics['ipSDE_d0res_AB_max'] = np.nan
            metrics['ipSDE_d0res_AB_min'] = np.nan
            metrics['ipSDE_d0res_BA_max'] = np.nan
            metrics['ipSDE_d0res_BA_min'] = np.nan
            metrics['lis_pde_AB'] = np.nan
            metrics['lis_pde_BA'] = np.nan
            metrics['lis_pde_avg'] = np.nan
            metrics['pDockQ2_pde_AB'] = np.nan
            metrics['pDockQ2_pde_BA'] = np.nan
        #
        # Store thresholds
        metrics['pAE_threshold'] = args.pae_threshold
        metrics['pDE_threshold'] = args.pde_threshold
        metrics['pDockQ_threshold'] = args.pdockq_threshold
        metrics['interface_threshold'] = args.interface_threshold
    #
    else:
        # Set all to np.nan if insufficient chains
        metrics['num_interface_residues_A'] = np.nan
        metrics['num_interface_residues_B'] = np.nan
        
        # ipLDDT metrics
        metrics['ipLDDT_mean'] = np.nan
        metrics['ipLDDT_min'] = np.nan
        metrics['ipLDDT_max'] = np.nan
        
        # interchain PAE metrics
        metrics['interchain_pAE_mean'] = np.nan
        metrics['interchain_pAE_min'] = np.nan
        metrics['interchain_pAE_max'] = np.nan
        
        # ipAE metrics
        metrics['ipAE_mean'] = np.nan
        metrics['ipAE_min'] = np.nan
        metrics['ipAE_max'] = np.nan
        metrics['ipAE_pct_lt_threshold'] = np.nan
        
        # ipSAE variants (12 metrics)
        for variant in ['d0res', 'd0dom', 'd0chn']:
            for direction in ['AB', 'BA']:
                for agg in ['max', 'min']:
                    metrics[f'ipSAE_{variant}_{direction}_{agg}'] = np.nan
        
        # pDockQ variants
        metrics['pDockQ'] = np.nan
        metrics['pDockQ2_AB'] = np.nan
        metrics['pDockQ2_BA'] = np.nan
        
        # LIS variants
        metrics['lis_AB'] = np.nan
        metrics['lis_BA'] = np.nan
        metrics['lis_avg'] = np.nan
        
        # PDE-based metrics
        metrics['interchain_pDE_mean'] = np.nan
        metrics['interchain_pDE_min'] = np.nan
        metrics['interchain_pDE_max'] = np.nan
        metrics['ipDE_mean'] = np.nan
        metrics['ipDE_min'] = np.nan
        metrics['ipDE_max'] = np.nan
        metrics['ipSDE_d0res_AB_max'] = np.nan
        metrics['ipSDE_d0res_AB_min'] = np.nan
        metrics['ipSDE_d0res_BA_max'] = np.nan
        metrics['ipSDE_d0res_BA_min'] = np.nan
        metrics['lis_pde_AB'] = np.nan
        metrics['lis_pde_BA'] = np.nan
        metrics['lis_pde_avg'] = np.nan
        metrics['pDockQ2_pde_AB'] = np.nan
        metrics['pDockQ2_pde_BA'] = np.nan
        
        # Thresholds
        metrics['pAE_threshold'] = np.nan
        metrics['pDE_threshold'] = np.nan
        metrics['pDockQ_threshold'] = np.nan
        metrics['interface_threshold'] = np.nan
    #
    return metrics


# ============================================================================
# Main function
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Extract metrics from Boltz-2 predictions')
    parser.add_argument('--parent_dir', type=str, help='Directory containing sample_dirs')
    parser.add_argument('--prefix', type=str, help='Prefix to come before all keynames')
    parser.add_argument('--pae_threshold', type=float, default=10.0, 
                        help='PAE cutoff for ipSAE/LIS calculation (default: 10.0)')
    parser.add_argument('--pde_threshold', type=float, default=2.0,
                        help='PDE cutoff for ipSDE/LIS-PDE calculation (default: 2.0)')
    parser.add_argument('--interface_threshold', type=float, default=5.0,
                        help='Distance cutoff for interface residues (default: 5.0)')
    parser.add_argument('--pdockq_threshold', type=float, default=8.0,
                        help='Distance cutoff for pDockQ calculation (default: 8.0)')
    parser.add_argument('--output_filename', type=str, default="boltz_metrics.csv")

    args = parser.parse_args()
    #
    parent_dir = Path(args.parent_dir).expanduser()
        
    if not parent_dir.exists():
        print(f"Error: Directory {parent_dir} does not exist", file=sys.stderr)
        return
    #
    all_metrics = []
    #
    for sample_dir in sorted(parent_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        name = Path(sample_dir).name
        try:
            metrics = extract_sample_metrics(sample_dir, args)
            if metrics is not None:
                all_metrics.append(metrics)
            print(f"Extracted metrics from {name}")
        except Exception as e:
            print(f"Error processing {name}: {e}", file=sys.stderr)
    #
    if len(all_metrics) == 0:
        print("No metrics extracted.", file=sys.stderr)
        return
    
    all_metrics = [
        {f"{args.prefix}_{k}": v for k, v in m.items()}
        for m in all_metrics
    ] 
    #
    if all_metrics:
        keys = list(all_metrics[0].keys())
    else:
        keys = []
    output_file = Path(args.output_filename)
    #
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for m in all_metrics:
            writer.writerow(m)
    #
    print(f"Wrote CSV to {output_file}")


if __name__ == '__main__':
    main()