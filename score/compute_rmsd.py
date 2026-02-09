#!/usr/bin/env python3
"""
compute_rmsd.py

Compute pairwise binder RMSD across structure prediction methods.
Aligns all models to a common reference frame using target chain (A) CA atoms,
then computes chain B (binder) CA RMSD for each method pair.

Usage:
    python compute_rmsd.py --project /path/to/project

Inputs:
    $project/relaxed/{method}/*.pdb  (relaxed structures from compute_rosetta_metrics.py)

Outputs:
    $project/RMSD.csv
"""

import argparse
import sys
import csv
from pathlib import Path
from collections import defaultdict
from itertools import combinations

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa


# ============================================================================
# Configuration
# ============================================================================

METHODS = ['af2_multimer', 'af2_ptm_complex', 'boltz', 'chai']
REFERENCE_METHOD = 'af2_multimer'

PREFIX_MAP = {
    'af2_multimer': 'mul',
    'af2_ptm_complex': 'mon',
    'boltz': 'bol',
    'chai': 'cha',
}

# All ordered pairs for column output
METHOD_PAIRS = list(combinations(METHODS, 2))

TARGET_CHAIN = 'A'
BINDER_CHAIN = 'B'


# ============================================================================
# Helper functions
# ============================================================================

def discover_binders(relaxed_dir):
    """
    Scan relaxed directories for *.pdb files.
    
    Returns:
        dict: {binder_id: {method: Path, ...}}
    """
    binder_map = defaultdict(dict)
    #
    for method in METHODS:
        method_dir = relaxed_dir / method
        if not method_dir.is_dir():
            continue
        #
        for pdb_file in sorted(method_dir.glob("*.pdb")):
            binder_id = pdb_file.stem
            binder_map[binder_id][method] = pdb_file
    #
    # Filter to binders with >= 2 methods
    binder_map = {
        bid: methods 
        for bid, methods in binder_map.items() 
        if len(methods) >= 2
    }
    #
    return binder_map


def extract_ca_coords(structure, chain_id):
    """
    Extract CA atom coordinates from a chain.
    
    Returns:
        list of Bio.PDB.Atom objects (CA atoms in residue order)
    """
    model = list(structure.get_models())[0]
    #
    if chain_id not in [c.id for c in model]:
        return None
    #
    chain = model[chain_id]
    ca_atoms = []
    #
    for res in chain:
        if is_aa(res, standard=True) and 'CA' in res:
            ca_atoms.append(res['CA'])
    #
    return ca_atoms


def get_ca_coord_array(ca_atoms):
    """Convert list of CA atoms to numpy coordinate array."""
    return np.array([atom.get_coord() for atom in ca_atoms])


def compute_rmsd(coords1, coords2):
    """Compute RMSD between two coordinate arrays of equal length."""
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def align_to_reference(ref_structure, mobile_structure):
    """
    Align mobile_structure onto ref_structure using chain A CA atoms.
    Applies the transformation to the entire mobile structure in-place.
    
    Returns:
        bool: True if alignment succeeded, False otherwise
    """
    ref_ca = extract_ca_coords(ref_structure, TARGET_CHAIN)
    mob_ca = extract_ca_coords(mobile_structure, TARGET_CHAIN)
    #
    if ref_ca is None or mob_ca is None:
        return False
    #
    if len(ref_ca) != len(mob_ca):
        return False
    #
    if len(ref_ca) == 0:
        return False
    #
    # Superimpose using chain A CA atoms
    sup = Superimposer()
    sup.set_atoms(ref_ca, mob_ca)
    #
    # Apply rotation/translation to ALL atoms in mobile structure
    all_atoms = list(mobile_structure.get_atoms())
    sup.apply(all_atoms)
    #
    return True


def process_binder(binder_id, method_paths, parser):
    """
    Process a single binder: load structures, align, compute pairwise RMSDs.
    
    Returns:
        dict: row for CSV output
    """
    # Load all available structures
    structures = {}
    for method, pdb_path in method_paths.items():
        try:
            # Use unique structure ID to avoid BioPython caching issues
            struct_id = f"{binder_id}_{method}"
            structure = parser.get_structure(struct_id, str(pdb_path))
            structures[method] = structure
        except Exception as e:
            print(f"  [WARN] Failed to parse {method} for {binder_id}: {e}")
    #
    if len(structures) < 2:
        return None
    #
    # Choose reference: af2_multimer if available, else first available
    if REFERENCE_METHOD in structures:
        ref_method = REFERENCE_METHOD
    else:
        ref_method = list(structures.keys())[0]
    #
    ref_structure = structures[ref_method]
    #
    # Align all structures to reference using chain A CA atoms
    aligned = {ref_method: ref_structure}
    for method, structure in structures.items():
        if method == ref_method:
            continue
        success = align_to_reference(ref_structure, structure)
        if success:
            aligned[method] = structure
        else:
            print(f"  [WARN] Alignment failed for {method}/{binder_id} "
                  f"(chain A CA mismatch with {ref_method})")
    #
    if len(aligned) < 2:
        return None
    #
    # Extract chain B CA coords from all aligned structures
    binder_coords = {}
    for method, structure in aligned.items():
        ca_atoms = extract_ca_coords(structure, BINDER_CHAIN)
        if ca_atoms is None or len(ca_atoms) == 0:
            print(f"  [WARN] No chain B CA atoms for {method}/{binder_id}")
            continue
        binder_coords[method] = get_ca_coord_array(ca_atoms)
    #
    if len(binder_coords) < 2:
        return None
    #
    # Compute pairwise binder RMSDs
    row = {'binder': binder_id}
    pairwise_values = []
    #
    for method1, method2 in METHOD_PAIRS:
        prefix1 = PREFIX_MAP[method1]
        prefix2 = PREFIX_MAP[method2]
        col_name = f"rmsd_B_{prefix1}__{prefix2}"
        #
        if method1 in binder_coords and method2 in binder_coords:
            coords1 = binder_coords[method1]
            coords2 = binder_coords[method2]
            #
            if len(coords1) == len(coords2):
                rmsd = compute_rmsd(coords1, coords2)
                row[col_name] = round(rmsd, 4)
                pairwise_values.append(rmsd)
            else:
                print(f"  [WARN] Chain B length mismatch for {binder_id}: "
                      f"{method1}={len(coords1)}, {method2}={len(coords2)}")
                row[col_name] = np.nan
        else:
            row[col_name] = np.nan
    #
    # Summary stats
    if pairwise_values:
        row['rmsd_B_mean'] = round(float(np.mean(pairwise_values)), 4)
        row['rmsd_B_min'] = round(float(np.min(pairwise_values)), 4)
        row['rmsd_B_max'] = round(float(np.max(pairwise_values)), 4)
    else:
        row['rmsd_B_mean'] = np.nan
        row['rmsd_B_min'] = np.nan
        row['rmsd_B_max'] = np.nan
    #
    return row


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compute pairwise binder RMSD across prediction methods'
    )
    parser.add_argument('--project', type=str, required=True,
                        help='Project directory')
    parser.add_argument('--output', type=str, default='RMSD.csv',
                        help='Output filename (default: RMSD.csv)')
    args = parser.parse_args()
    #
    project = Path(args.project)
    relaxed_dir = project / 'relaxed'
    #
    if not relaxed_dir.is_dir():
        print(f"[ERROR] Relaxed directory not found: {relaxed_dir}", file=sys.stderr)
        return 1
    #
    # Step 1: Discover binders
    print("[INFO] Discovering binder structures...")
    binder_map = discover_binders(relaxed_dir)
    #
    if not binder_map:
        print("[ERROR] No binders found with >= 2 methods", file=sys.stderr)
        return 1
    #
    # Report discovery
    method_counts = defaultdict(int)
    for bid, methods in binder_map.items():
        for m in methods:
            method_counts[m] += 1
    print(f"[INFO] Found {len(binder_map)} binders with >= 2 methods")
    for method in METHODS:
        print(f"[INFO]   {method}: {method_counts.get(method, 0)} structures")
    #
    # Step 2-5: Process each binder
    print("[INFO] Computing pairwise RMSDs...")
    pdb_parser = PDBParser(QUIET=True)
    all_rows = []
    #
    for binder_id in sorted(binder_map.keys()):
        method_paths = binder_map[binder_id]
        row = process_binder(binder_id, method_paths, pdb_parser)
        if row is not None:
            all_rows.append(row)
            avail = [PREFIX_MAP[m] for m in METHODS if m in method_paths]
            print(f"  [OK] {binder_id} ({','.join(avail)})")
        else:
            print(f"  [SKIP] {binder_id} (insufficient data after alignment)")
    #
    if not all_rows:
        print("[ERROR] No RMSD results computed", file=sys.stderr)
        return 1
    #
    # Step 6: Write CSV
    output_path = project / args.output
    #
    # Define column order
    pairwise_cols = [
        f"rmsd_B_{PREFIX_MAP[m1]}__{PREFIX_MAP[m2]}" 
        for m1, m2 in METHOD_PAIRS
    ]
    fieldnames = ['binder'] + pairwise_cols + ['rmsd_B_mean', 'rmsd_B_min', 'rmsd_B_max']
    #
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)
    #
    print(f"[INFO] Wrote {len(all_rows)} rows to {output_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
