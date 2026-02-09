#!/usr/bin/env python3
# ============================================================================
# chai_extract_metrics.py
# Extract interface quality metrics from Chai-1 predictions
# Simplified version - only metrics that don't require PAE matrix
# ============================================================================

import argparse
import sys
import csv
from pathlib import Path

import math
import numpy as np

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBIO

# ============================================================================
# Helper functions - Structure parsing
# ============================================================================

def parse_cif_structure(structure_file):
    """
    Parse CIF structure file
    """
    parser = MMCIFParser(QUIET=True)
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


def extract_plddt_from_cif(structure, chains, chain_residues):
    """
    Extract pLDDT from CIF B-factors
    Returns array in same order as chain_residues
    """
    plddt_array = []
    #
    for chain_id in chains:
        residues = chain_residues[chain_id]
        for res in residues:
            # Get CA atom b-factor as pLDDT
            if 'CA' in res:
                plddt_array.append(res['CA'].bfactor)
            else:
                # Fallback to first atom
                atoms = list(res.get_atoms())
                if atoms:
                    plddt_array.append(atoms[0].bfactor)
                else:
                    plddt_array.append(0.0)
    #
    return np.array(plddt_array)


def get_interface(chain_residues, chain_a, chain_b, interface_threshold):
    """
    Find interface residues between two chains
    """
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
    """
    Map chain IDs to residue indices
    """
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


# ============================================================================
# Helper functions - Metric computation (simplified)
# ============================================================================

def compute_interface_metrics_simple(plddt_array, interface_a, interface_b,
                                     chain_a_start, chain_b_start):
    """
    Compute interface metrics - only ipLDDT (no PAE metrics)
    """
    metrics = {}
    #
    # ipLDDT only
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
    return metrics


def compute_pdockq_simple(cb_coords_a, cb_coords_b, plddt_a, plddt_b, distance_cutoff):
    """
    Compute pDockQ score (no pDockQ2 - that requires PAE)
    """
    # Compute distance matrix
    distances = compute_cb_distances(cb_coords_a, cb_coords_b)
    contacts = distances <= distance_cutoff
    #
    if not np.any(contacts):
        return 0.0
    #
    # Find unique interface residues
    interface_residues_a = set(np.where(np.any(contacts, axis=1))[0])
    interface_residues_b = set(np.where(np.any(contacts, axis=0))[0])
    interface_residues = list(interface_residues_a) + [i + len(plddt_a) for i in interface_residues_b]
    #
    # Combine pLDDT arrays
    combined_plddt = np.concatenate([plddt_a, plddt_b])
    mean_plddt = np.mean(combined_plddt[interface_residues])
    npairs = np.sum(contacts)
    #
    if npairs > 0:
        x = mean_plddt * math.log10(npairs)
        pdockq = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
        return float(pdockq)
    else:
        return 0.0


# ============================================================================
# Main extraction function
# ============================================================================

def extract_sample_metrics(sample_dir, args):
    """
    Extract metrics for a single Chai-1 sample directory
    """
    #
    sample_dir = Path(sample_dir)
    sample_name = sample_dir.name
    #
    # Step 1: Find files
    npz_files = list(sample_dir.glob("*.npz"))
    if not npz_files:
        return None
    npz_file = npz_files[0]
    #
    cif_files = list(sample_dir.glob("*.cif"))
    if not cif_files:
        return None
    cif_file = cif_files[0]
    #
    # Step 2: Load NPZ scores
    npz_data = np.load(npz_file)
    #
    metrics = {
        'sample_name': sample_name,
        'ptm': float(npz_data['ptm'][0]),
        'iptm': float(npz_data['iptm'][0]),
        'aggregate_score': float(npz_data['aggregate_score'][0]),
        'has_inter_chain_clashes': bool(npz_data['has_inter_chain_clashes'][0]),
    }
    #
    # Step 3: Parse CIF structure
    structure, chains, chain_residues = parse_cif_structure(cif_file)
    #
    # Extract pLDDT from CIF B-factors
    plddt_array = extract_plddt_from_cif(structure, chains, chain_residues)
    metrics['pLDDT_mean'] = float(np.mean(plddt_array))
    #
    # Step 4: Extract per-chain metrics
    if len(chains) >= 2:
        # Per-chain PTM scores
        per_chain_ptm = npz_data['per_chain_ptm'][0]
        metrics['per_chain_ptm_A'] = float(per_chain_ptm[0])
        metrics['per_chain_ptm_B'] = float(per_chain_ptm[1])
        #
        # Per-chain-pair iPTM scores (directional)
        per_chain_pair_iptm = npz_data['per_chain_pair_iptm'][0]
        metrics['per_chain_pair_iptm_AB'] = float(per_chain_pair_iptm[0, 1])
        metrics['per_chain_pair_iptm_BA'] = float(per_chain_pair_iptm[1, 0])
    else:
        metrics['per_chain_ptm_A'] = np.nan
        metrics['per_chain_ptm_B'] = np.nan
        metrics['per_chain_pair_iptm_AB'] = np.nan
        metrics['per_chain_pair_iptm_BA'] = np.nan
    #
    # Step 5: Compute interface metrics
    if len(chains) >= 2:
        #
        chain_a, chain_b = chains[0], chains[1]
        #
        # Get chain mapping
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
        # Compute ipLDDT metrics
        interface_metrics = compute_interface_metrics_simple(
            plddt_array, interface_a, interface_b,
            chain_a_start, chain_b_start
        )
        metrics.update(interface_metrics)
        #
        # Extract CB positions
        cb_coords_a, plddt_a = get_cb_positions(chain_residues, chain_a)
        cb_coords_b, plddt_b = get_cb_positions(chain_residues, chain_b)
        #
        # Compute pDockQ
        if len(cb_coords_a) > 0 and len(cb_coords_b) > 0:
            pdockq = compute_pdockq_simple(
                cb_coords_a, cb_coords_b,
                plddt_a, plddt_b,
                args.pdockq_threshold
            )
            metrics['pDockQ'] = pdockq
        else:
            metrics['pDockQ'] = np.nan
        #
        # Store thresholds
        metrics['pDockQ_threshold'] = args.pdockq_threshold
        metrics['interface_threshold'] = args.interface_threshold
    #
    else:
        # Set all to np.nan if insufficient chains
        metrics['num_interface_residues_A'] = np.nan
        metrics['num_interface_residues_B'] = np.nan
        metrics['ipLDDT_mean'] = np.nan
        metrics['ipLDDT_min'] = np.nan
        metrics['ipLDDT_max'] = np.nan
        metrics['pDockQ'] = np.nan
        metrics['pDockQ_threshold'] = np.nan
        metrics['interface_threshold'] = np.nan
    #
    return metrics


# ============================================================================
# Main function
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Extract metrics from Chai-1 predictions')
    parser.add_argument('--parent_dir', type=str, help='Directory containing sample_dirs')
    parser.add_argument('--prefix', type=str, help='Prefix to come before all keynames')
    parser.add_argument('--interface_threshold', type=float, default=5.0,
                        help='Distance cutoff for interface residues (default: 5.0)')
    parser.add_argument('--pdockq_threshold', type=float, default=8.0,
                        help='Distance cutoff for pDockQ calculation (default: 8.0)')
    parser.add_argument('--output_filename', type=str, default="chai_metrics.csv")
    args = parser.parse_args()
    #
    parent_dir = Path(args.parent_dir).expanduser()
    #
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
    #
    all_metrics = [
    {f"{args.prefix}_{k}": v for k, v in m.items()}
    for m in all_metrics
    ] 
    #
    # Write CSV
    if all_metrics:
        keys = list(all_metrics[0].keys())
    else:
        keys = []
    #
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