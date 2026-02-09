#!/usr/bin/env python3
"""
calculate_distances.py - Calculate distances from target residue to binder for all predictions

Processes predictions from:
- AF2-multimer (ColabFold)
- AF2-ptm complex (ColabFold)
- Chai-1
- Boltz-2

Note: Binder-only predictions are skipped (single chain, no distance calculation needed).

For each structure, calculates the minimum heavy-atom distance between a specified
target residue (looked up per-target from residue CSV) and any residue in the binder chain.

Usage:
    python calculate_distances.py \\
        --predictions_dir /path/to/predictions \\
        --residue_csv /path/to/residues.csv \\
        --csv /path/to/input.csv \\
        --output_csv /path/to/distances.csv

Residue CSV format:
    target_name,residue_number
    AQP,122
    SLP,85
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Optional

# Try to import structure parsing libraries
try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("WARNING: BioPython not available, distance calculations will fail")


def extract_target_name(binder_name: str) -> str:
    """
    Extract target name from binder name.
    
    Binder names must follow the format TARGET_NUMBER (e.g., AQP_01, AQP2_01, IL6R_123).
    The target name is everything before the last underscore, and the suffix must be
    purely numeric.
    
    # NOTE: keep in sync with generate_boltz_yaml.py
    
    Examples:
        AQP_01 -> AQP
        AQP2_01 -> AQP2
        IL6R_123 -> IL6R
        my_target_5 -> my_target
    """
    match = re.match(r'^(.+)_(\d+)$', binder_name)
    if match:
        return match.group(1)
    
    raise ValueError(
        f"Binder name '{binder_name}' does not match expected TARGET_NUMBER format "
        f"(e.g., AQP_01, AQP2_01, IL6R_123)"
    )


def load_residue_map(residue_csv: Path) -> dict[str, int]:
    """
    Load target name to residue number mapping from CSV.
    
    CSV format: target_name,residue_number
    Returns dict like {"AQP": 122, "SLP": 85}
    """
    residue_map = {}
    with open(residue_csv, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Handle various column name possibilities
            target = None
            residue = None
            
            for key in row:
                key_lower = key.lower().strip()
                if key_lower in ('target_name', 'target', 'name'):
                    target = row[key].strip()
                elif key_lower in ('residue_number', 'residue', 'resnum', 'res'):
                    try:
                        residue = int(row[key].strip())
                    except ValueError:
                        pass
            
            if target and residue is not None:
                residue_map[target] = residue
                # Also store lowercase version for case-insensitive lookup
                residue_map[target.lower()] = residue
    
    return residue_map


def load_structure(path: Path):
    """Load a structure from PDB or CIF file."""
    if not HAS_BIOPYTHON:
        raise ImportError("BioPython required for structure parsing")
    
    path_str = str(path).lower()
    if path_str.endswith(('.cif', '.mmcif')):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    
    return parser.get_structure("complex", str(path))


def get_chain(model: Model, chain_id: str) -> Optional[Chain]:
    """Get a chain by ID from a model."""
    for ch in model.get_chains():
        if ch.id == chain_id:
            return ch
    return None


def find_residue_by_number(chain: Chain, resnum: int) -> Optional[Residue]:
    """Find a residue by its sequence number."""
    for res in chain.get_residues():
        # Skip heteroatoms
        if res.id[0] != " ":
            continue
        if res.id[1] == resnum:
            return res
    return None


def is_hydrogen(atom) -> bool:
    """Check if an atom is hydrogen."""
    el = (atom.element or "").strip().upper()
    if el == "H":
        return True
    return atom.get_name().strip().upper().startswith("H")


def min_heavyatom_distance(res1: Residue, res2: Residue) -> float:
    """Calculate minimum heavy-atom distance between two residues."""
    min_d = float("inf")
    for a1 in res1.get_atoms():
        if is_hydrogen(a1):
            continue
        for a2 in res2.get_atoms():
            if is_hydrogen(a2):
                continue
            d = a1 - a2
            if d < min_d:
                min_d = d
    return min_d


def calculate_min_distance_to_binder(
    structure_path: Path,
    target_residue: int,
    target_chain: str = "A",
    binder_chain: str = "B",
) -> Optional[float]:
    """
    Calculate minimum distance from target residue to any binder residue.
    
    Returns None if structure can't be parsed or residue not found.
    """
    try:
        structure = load_structure(structure_path)
        model = next(structure.get_models())
        
        target_ch = get_chain(model, target_chain)
        binder_ch = get_chain(model, binder_chain)
        
        if target_ch is None or binder_ch is None:
            # Try to find chains by order if named differently
            chains = list(model.get_chains())
            if len(chains) >= 2:
                target_ch = chains[0]
                binder_ch = chains[1]
            else:
                return None
        
        target_res = find_residue_by_number(target_ch, target_residue)
        if target_res is None:
            return None
        
        min_dist = float("inf")
        for res in binder_ch.get_residues():
            if res.id[0] != " ":  # Skip heteroatoms
                continue
            d = min_heavyatom_distance(target_res, res)
            if d < min_dist:
                min_dist = d
        
        return min_dist if min_dist < float("inf") else None
    
    except Exception as e:
        print(f"WARNING: Failed to process {structure_path}: {e}")
        return None


def find_structures_af2(base_dir: Path, binder_names: list[str], method: str) -> list[dict]:
    """
    Find all AF2 structure files (af2_multimer or af2_ptm_complex).
    
    Directory: base_dir/ (flat)
    Pattern: {binder}_unrelaxed_rank_00{N}_*.pdb
    Rank: N (1-5)
    """
    results = []
    
    for binder in binder_names:
        # Find all unrelaxed PDB files for this binder
        pattern = f"{binder}_unrelaxed_rank_*.pdb"
        pdb_files = list(base_dir.glob(pattern))
        
        for pdb_file in sorted(pdb_files):
            # Extract rank from filename (e.g., rank_001 -> 1)
            match = re.search(r'rank_(\d+)', pdb_file.name)
            if match:
                rank = int(match.group(1))
                results.append({
                    "binder": binder,
                    "method": method,
                    "model_rank": rank,
                    "structure_file": str(pdb_file),
                })
    
    return results


def find_structures_chai(base_dir: Path, binder_names: list[str]) -> list[dict]:
    """
    Find all Chai-1 structure files.
    
    Directory: base_dir/{binder}/
    Pattern: pred.model_idx_{N}.cif
    Rank: N (0-4)
    """
    results = []
    
    for binder in binder_names:
        binder_dir = base_dir / binder
        if not binder_dir.exists():
            continue
        
        # Find all CIF files for this binder
        cif_files = list(binder_dir.glob("pred.model_idx_*.cif"))
        
        for cif_file in sorted(cif_files):
            # Extract model index from filename
            match = re.search(r'model_idx_(\d+)', cif_file.name)
            if match:
                model_idx = int(match.group(1))
                results.append({
                    "binder": binder,
                    "method": "chai",
                    "model_rank": model_idx,
                    "structure_file": str(cif_file),
                })
    
    return results


def find_structures_boltz(base_dir: Path, binder_names: list[str]) -> list[dict]:
    """
    Find all Boltz-2 structure files.
    
    Directory: base_dir/boltz_results_yaml/predictions/{binder}/
    Pattern: {binder}_model_{N}.pdb
    Rank: N (0-4)
    """
    results = []
    
    # Boltz outputs to boltz_results_yaml/predictions/
    predictions_dir = base_dir / "boltz_results_yaml" / "predictions"
    if not predictions_dir.exists():
        # Fallback: try base_dir/predictions
        predictions_dir = base_dir / "predictions"
        if not predictions_dir.exists():
            return results
    
    for binder in binder_names:
        binder_dir = predictions_dir / binder
        if not binder_dir.exists():
            continue
        
        # Find all PDB files for this binder
        pdb_files = list(binder_dir.glob(f"{binder}_model_*.pdb"))
        
        for pdb_file in sorted(pdb_files):
            # Extract model index from filename
            match = re.search(r'model_(\d+)', pdb_file.name)
            if match:
                model_idx = int(match.group(1))
                results.append({
                    "binder": binder,
                    "method": "boltz",
                    "model_rank": model_idx,
                    "structure_file": str(pdb_file),
                })
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Calculate distances from target residue to binder chain"
    )
    parser.add_argument("--predictions_dir", required=True, 
                        help="Directory containing prediction subdirectories")
    parser.add_argument("--residue_csv", required=True,
                        help="CSV file with target_name,residue_number columns")
    parser.add_argument("--csv", required=True,
                        help="Input CSV file to get binder names")
    parser.add_argument("--output_csv", required=True,
                        help="Output CSV file for distances")
    parser.add_argument("--target_chain", default="A",
                        help="Target chain ID (default: A)")
    parser.add_argument("--binder_chain", default="B",
                        help="Binder chain ID (default: B)")
    parser.add_argument("--attempt", type=int, default=1,
                        help="Attempt number (1 or 2) - added to output CSV")
    
    args = parser.parse_args()
    
    if not HAS_BIOPYTHON:
        print("ERROR: BioPython is required for distance calculations")
        return 1
    
    predictions_dir = Path(args.predictions_dir)
    residue_csv = Path(args.residue_csv)
    output_csv = Path(args.output_csv)
    
    # Load residue map
    residue_map = load_residue_map(residue_csv)
    if not residue_map:
        print(f"ERROR: No residue mappings found in {residue_csv}")
        return 1
    print(f"Loaded {len(residue_map) // 2} target residue mappings")  # div by 2 due to lowercase copies
    
    # Get binder names from CSV
    with open(args.csv, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        binder_names = [row['name'].strip() for row in reader if row.get('name', '').strip()]
    
    print(f"Found {len(binder_names)} binders in CSV")
    
    # Process each method (excluding binder - single chain, no distance calc)
    all_results = []
    
    # AF2-multimer
    af2_multimer_dir = predictions_dir / "af2_multimer"
    if af2_multimer_dir.exists():
        print("Processing af2_multimer...")
        structures = find_structures_af2(af2_multimer_dir, binder_names, "af2_multimer")
        print(f"  Found {len(structures)} structures")
        all_results.extend(structures)
    else:
        print("Skipping af2_multimer: directory not found")
    
    # AF2-ptm complex
    af2_ptm_dir = predictions_dir / "af2_ptm_complex"
    if af2_ptm_dir.exists():
        print("Processing af2_ptm_complex...")
        structures = find_structures_af2(af2_ptm_dir, binder_names, "af2_ptm_complex")
        print(f"  Found {len(structures)} structures")
        all_results.extend(structures)
    else:
        print("Skipping af2_ptm_complex: directory not found")
    
    # Chai-1
    chai_dir = predictions_dir / "chai"
    if chai_dir.exists():
        print("Processing chai...")
        structures = find_structures_chai(chai_dir, binder_names)
        print(f"  Found {len(structures)} structures")
        all_results.extend(structures)
    else:
        print("Skipping chai: directory not found")
    
    # Boltz-2
    boltz_dir = predictions_dir / "boltz"
    if boltz_dir.exists():
        print("Processing boltz...")
        structures = find_structures_boltz(boltz_dir, binder_names)
        print(f"  Found {len(structures)} structures")
        all_results.extend(structures)
    else:
        print("Skipping boltz: directory not found")
    
    # Calculate distances for all structures
    print(f"\nCalculating distances for {len(all_results)} structures...")
    
    for struct in all_results:
        binder_name = struct["binder"]
        
        # Look up target residue for this binder
        target_name = extract_target_name(binder_name)
        target_residue = residue_map.get(target_name) or residue_map.get(target_name.lower())
        
        if target_residue is None:
            print(f"  WARNING: No residue mapping for target '{target_name}' (binder: {binder_name})")
            struct["min_distance"] = None
            struct["target_residue"] = None
        else:
            # Calculate distance
            distance = calculate_min_distance_to_binder(
                Path(struct["structure_file"]),
                target_residue,
                args.target_chain,
                args.binder_chain,
            )
            struct["min_distance"] = distance
            struct["target_residue"] = target_residue
        
        # Add attempt number
        struct["attempt"] = args.attempt
    
    # Write output CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    
    fieldnames = ["binder", "method", "model_rank", "structure_file", "target_residue", "min_distance", "attempt"]
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    
    print(f"\nWrote {len(all_results)} rows to {output_csv}")
    
    # Summary statistics
    successful = sum(1 for r in all_results if r["min_distance"] is not None)
    print(f"Successfully calculated distances for {successful}/{len(all_results)} structures")
    
    # Per-method breakdown
    from collections import defaultdict
    method_counts = defaultdict(lambda: {"total": 0, "success": 0})
    for r in all_results:
        method_counts[r["method"]]["total"] += 1
        if r["min_distance"] is not None:
            method_counts[r["method"]]["success"] += 1
    
    print(f"\nPer-method breakdown:")
    for method in sorted(method_counts.keys()):
        stats = method_counts[method]
        print(f"  {method}: {stats['success']}/{stats['total']} distances calculated")
    
    return 0


if __name__ == "__main__":
    exit(main())
