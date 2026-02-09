#!/usr/bin/env python3
"""
binder_extract_metrics.py

Extract metrics from single-chain binder structure predictions.
Processes a directory of sample subdirectories, each containing structure 
(CIF or PDB) and scores (NPZ or JSON) files for a single binder. 
Outputs a single CSV with one row per sample.

Supports multiple file formats:
- Structure files: .cif (preferred) or .pdb
- Scores files: .npz (preferred) or .json

Chain ID is fixed to "A" for single-chain binders.

Extracted metrics:
- ptm: Predicted TM-score (from both NPZ and JSON)
- aggregate_score: Aggregate confidence score (NPZ only)
- has_clashes: Inter-chain clash detection (NPZ only)
- plddt_min/mean/max: Per-residue confidence statistics
- pae_min/mean/max: Predicted aligned error statistics (JSON only)
- helix_pct/sheet_pct/loop_pct: Secondary structure composition
- Optional Rosetta metrics (with --compute_rosetta flag)

Usage:
    python binder_extract_metrics.py --parent_dir /path/to/predictions/
    python binder_extract_metrics.py --parent_dir /path/to/predictions/ --compute_rosetta
"""

import argparse
import sys
import csv
import json
from pathlib import Path
import numpy as np

from Bio.PDB import PDBParser, MMCIFParser, DSSP
from Bio.PDB.Polypeptide import is_aa



# ============================================================================
# Helper Functions
# ============================================================================

def parse_structure_file(structure_file):
    """
    Parse CIF or PDB file and extract structure information.
    Auto-detects format based on file extension.
    
    Parameters:
    -----------
    structure_file : str/Path
        Path to CIF or PDB file
    
    Returns:
    --------
    tuple : (structure, chains, chain_residues)
        - structure: Bio.PDB Structure object
        - chains: list of chain IDs
        - chain_residues: dict mapping chain_id -> list of residues
    """
    structure_file = Path(structure_file)
    
    # Select appropriate parser based on file extension
    if structure_file.suffix.lower() == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:  # .pdb or other
        parser = PDBParser(QUIET=True)
    
    structure = parser.get_structure('binder', str(structure_file))
    model = list(structure.get_models())[0]
    
    chains = []
    chain_residues = {}
    
    for chain in model:
        chains.append(chain.id)
        residues = [res for res in chain if is_aa(res, standard=True)]
        chain_residues[chain.id] = residues
    
    return structure, chains, chain_residues


def load_scores_file(scores_file):
    """
    Load scores from either NPZ or JSON file.
    Auto-detects format based on file extension.
    
    Parameters:
    -----------
    scores_file : str/Path
        Path to NPZ or JSON file
    
    Returns:
    --------
    dict : Dictionary with keys:
        - ptm: float
        - plddt: array or None
        - pae: 2D array or None
        - aggregate_score: float or None (NPZ only)
        - has_clashes: bool or None (NPZ only)
        - max_pae: float or None (JSON only)
    """
    scores_file = Path(scores_file)
    
    if scores_file.suffix.lower() == '.npz':
        # Load NPZ file
        data = np.load(scores_file)
        
        # Extract PTM (handle both scalar and array formats)
        ptm_value = data.get('ptm', np.nan)
        if isinstance(ptm_value, np.ndarray) and ptm_value.size > 0:
            ptm_value = float(ptm_value.flat[0])
        elif isinstance(ptm_value, np.ndarray):
            ptm_value = np.nan
        else:
            ptm_value = float(ptm_value)
        
        # Extract aggregate_score (NPZ-specific)
        aggregate_score = data.get('aggregate_score', None)
        if aggregate_score is not None:
            if isinstance(aggregate_score, np.ndarray) and aggregate_score.size > 0:
                aggregate_score = float(aggregate_score.flat[0])
            else:
                aggregate_score = float(aggregate_score) if aggregate_score is not None else None
        
        # Extract has_inter_chain_clashes (NPZ-specific)
        has_clashes = data.get('has_inter_chain_clashes', None)
        if has_clashes is not None:
            if isinstance(has_clashes, np.ndarray) and has_clashes.size > 0:
                has_clashes = bool(has_clashes.flat[0])
            else:
                has_clashes = bool(has_clashes) if has_clashes is not None else None
        
        # pLDDT and PAE not typically in NPZ, will be extracted from structure
        return {
            'ptm': ptm_value,
            'plddt': None,
            'pae': None,
            'aggregate_score': aggregate_score,
            'has_clashes': has_clashes,
            'max_pae': None
        }
    
    else:  # .json
        with open(scores_file) as f:
            scores = json.load(f)
        
        return {
            'ptm': scores.get('ptm', np.nan),
            'plddt': scores.get('plddt', None),
            'pae': scores.get('pae', None),
            'aggregate_score': None,
            'has_clashes': None,
            'max_pae': scores.get('max_pae', None)
        }


def extract_plddt_from_structure(structure, chain_id):
    """
    Extract pLDDT values from structure B-factors.
    AlphaFold stores pLDDT in the B-factor field.
    
    Parameters:
    -----------
    structure : Bio.PDB Structure object
        Parsed structure
    chain_id : str
        Chain ID to extract from
    
    Returns:
    --------
    list or None : List of pLDDT values, or None if extraction fails
    """
    try:
        plddt_list = []
        model = list(structure.get_models())[0]
        
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if is_aa(residue, standard=True):
                        # Extract B-factor from CA atom
                        if 'CA' in residue:
                            ca_atom = residue['CA']
                            plddt_list.append(ca_atom.get_bfactor())
        
        return plddt_list if plddt_list else None
    
    except Exception as e:
        print(f"Warning: Failed to extract pLDDT from structure: {e}")
        return None


def compute_secondary_structure(structure_file, chain_id='A'):
    """
    Compute secondary structure composition using DSSP.
    Handles both CIF and PDB files.
    
    Parameters:
    -----------
    structure_file : str/Path
        Path to CIF or PDB file
    chain_id : str
        Chain ID to analyze
    
    Returns:
    --------
    dict : Dictionary with helix_pct, sheet_pct, loop_pct
    """
    try:
        structure_file = Path(structure_file)
        
        # Select appropriate parser based on file extension
        if structure_file.suffix.lower() == '.cif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure('binder', str(structure_file))
        model = structure[0]
        
        # Run DSSP
        dssp = DSSP(model, str(structure_file), dssp='mkdssp')
        
        # Count secondary structure elements
        helix_count = 0
        sheet_count = 0
        loop_count = 0
        
        for key in dssp:
            ss = key[2]
            
            # H = alpha helix, G = 3-10 helix, I = pi helix
            if ss in ['H', 'G', 'I']:
                helix_count += 1
            # E = beta strand, B = beta bridge
            elif ss in ['E', 'B']:
                sheet_count += 1
            # Everything else is loop/coil
            else:
                loop_count += 1
        
        total = helix_count + sheet_count + loop_count
        
        if total == 0:
            return {
                'helix_pct': np.nan,
                'sheet_pct': np.nan,
                'loop_pct': np.nan
            }
        
        return {
            'helix_pct': round(100.0 * helix_count / total, 2),
            'sheet_pct': round(100.0 * sheet_count / total, 2),
            'loop_pct': round(100.0 * loop_count / total, 2)
        }
    
    except Exception as e:
        print(f"Warning: DSSP calculation failed: {e}")
        return {
            'helix_pct': np.nan,
            'sheet_pct': np.nan,
            'loop_pct': np.nan
        }


def clean_pdb(pdb_file):
    """
    Remove non-essential lines from PDB file.
    Modifies file in-place.
    """
    with open(pdb_file, 'r') as f_in:
        relevant = [ln for ln in f_in if ln.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END'))]
    with open(pdb_file, 'w') as f_out:
        f_out.writelines(relevant)

# ============================================================================
# Main Sample Extraction Function
# ============================================================================

def extract_sample_metrics(sample_dir, args):
    """
    Extract all metrics for a single sample directory.
    
    Parameters:
    -----------
    sample_dir : Path
        Path to sample directory containing structure and scores files
    args : argparse.Namespace
        Command line arguments
    
    Returns:
    --------
    dict : Dictionary with all metrics, or None if extraction fails
    """
    sample_dir = Path(sample_dir)
    sample_name = sample_dir.name
    
    # Step 1: Find structure file (CIF or PDB)
    cif_files = list(sample_dir.glob("*.cif"))
    pdb_files = list(sample_dir.glob("*.pdb"))
    
    if cif_files:
        structure_file = cif_files[0]
    elif pdb_files:
        structure_file = pdb_files[0]
    else:
        print(f"Warning: No CIF or PDB file found in {sample_name}")
        return None
    
    # Step 2: Find scores file (NPZ or JSON)
    npz_files = list(sample_dir.glob("*.npz"))
    json_files = list(sample_dir.glob("*_scores_*.json"))
    if not json_files:
        json_files = list(sample_dir.glob("*.json"))
    
    if npz_files:
        scores_file = npz_files[0]
    elif json_files:
        scores_file = json_files[0]
    else:
        print(f"Warning: No NPZ or JSON file found in {sample_name}")
        return None
    
    # Step 3: Load scores using unified loader
    try:
        scores = load_scores_file(scores_file)
    except Exception as e:
        print(f"Error: Failed to load scores for {sample_name}: {e}")
        return None
    
    metrics = {
        'sample_name': sample_name,
        'ptm': scores.get('ptm', np.nan),
    }
    
    # Add NPZ-specific metrics
    aggregate_score = scores.get('aggregate_score', None)
    metrics['aggregate_score'] = aggregate_score if aggregate_score is not None else np.nan
    
    has_clashes = scores.get('has_clashes', None)
    metrics['has_clashes'] = has_clashes if has_clashes is not None else np.nan
    
    # Step 4: Parse structure using unified parser
    # Chain ID is always "A" for single-chain binders
    chain_id = "A"
    
    try:
        structure, chains, chain_residues = parse_structure_file(structure_file)
        
        # Verify chain A exists
        if chain_id not in chains:
            print(f"Warning: Chain {chain_id} not found in {sample_name}. Available chains: {chains}")
            if chains:
                chain_id = chains[0]
                print(f"  Using chain {chain_id} instead")
            else:
                print(f"Error: No valid chain found in {sample_name}")
                return None
        
    except Exception as e:
        print(f"Error: Failed to parse structure for {sample_name}: {e}")
        return None
    
    # Step 5: Extract pLDDT statistics
    plddt_array = scores.get('plddt', None)
    
    # If pLDDT not in scores file, try extracting from structure B-factors
    if plddt_array is None:
        plddt_list = extract_plddt_from_structure(structure, chain_id)
        if plddt_list:
            plddt_array = np.array(plddt_list)
    
    if plddt_array is not None:
        plddt_array = np.array(plddt_array)
        metrics['plddt_min'] = round(float(np.min(plddt_array)), 2)
        metrics['plddt_mean'] = round(float(np.mean(plddt_array)), 2)
        metrics['plddt_max'] = round(float(np.max(plddt_array)), 2)
    else:
        metrics['plddt_min'] = np.nan
        metrics['plddt_mean'] = np.nan
        metrics['plddt_max'] = np.nan
    
    # Step 6: Extract PAE statistics
    # For JSON: calculate min/mean/max from PAE matrix
    # For NPZ: set to NaN (no PAE data)
    pae_matrix = scores.get('pae', None)
    if pae_matrix is not None:
        pae_matrix = np.array(pae_matrix)
        metrics['pae_min'] = round(float(np.min(pae_matrix)), 2)
        metrics['pae_mean'] = round(float(np.mean(pae_matrix)), 2)
        metrics['pae_max'] = round(float(np.max(pae_matrix)), 2)
    else:
        # NPZ case - no PAE data available
        metrics['pae_min'] = np.nan
        metrics['pae_mean'] = np.nan
        metrics['pae_max'] = np.nan
    
    # Step 7: Compute secondary structure composition
    ss_metrics = compute_secondary_structure(structure_file, chain_id)
    metrics.update(ss_metrics)
    
    return metrics


# ============================================================================
# Main Execution
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Extract metrics from single-chain binder predictions (Chain ID fixed to "A")'
    )
    parser.add_argument('--parent_dir', type=str, required=True,
                        help='Directory containing sample subdirectories')
    parser.add_argument('--prefix', type=str, help='Prefix to come before all keynames')
    parser.add_argument('--output_filename', type=str, default='binder_metrics.csv',
                        help='Output CSV filename (default: binder_metrics.csv)')
    
    args = parser.parse_args()
    
    # Validate parent directory
    parent_dir = Path(args.parent_dir).expanduser()
    if not parent_dir.exists():
        print(f"Error: Directory {parent_dir} does not exist", file=sys.stderr)
        return
    # Process all sample directories
    all_metrics = []
    
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
    
    if len(all_metrics) == 0:
        print("No metrics extracted.", file=sys.stderr)
        return
    #
    all_metrics = [
    {f"{args.prefix}_{k}": v for k, v in m.items()}
    for m in all_metrics
    ] 
    # Write CSV
    output_file = Path(args.output_filename)
    
    if all_metrics:
        keys = list(all_metrics[0].keys())
    else:
        keys = []
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for m in all_metrics:
            writer.writerow(m)
    
    print(f"\nWrote {len(all_metrics)} rows to {output_file}")


if __name__ == '__main__':
    main()