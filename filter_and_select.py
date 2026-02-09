#!/usr/bin/env python3
"""
filter_and_select.py - Filter structures by distance and select best by rank

For each (binder, method) combination:
1. Filter structures where min_distance <= cutoff
2. Among passing structures, select the one with lowest rank (rank 1 for AF2, model 0 for others)
3. If NO structures pass cutoff, select the lowest-ranked structure as fallback
4. Copy selected structure AND all related files to output directory
5. Generate summary CSV with selection status

For binder-only predictions (single chain), always copy the highest-ranked model (model 0).

Supports merging attempt 1 and attempt 2 results:
- If attempt 1 passed, use attempt 1
- If attempt 1 failed and attempt 2 exists, use attempt 2 result

Usage:
    python filter_and_select.py \\
        --distances_csv /path/to/distances.csv \\
        --predictions_dir /path/to/predictions \\
        --cutoff 20.0 \\
        --output_dir /path/to/filtered_structures \\
        --summary_csv /path/to/summary.csv \\
        [--distances_csv_retry /path/to/distances_retry.csv] \\
        [--predictions_dir_retry /path/to/predictions_retry]
"""

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Optional
import re


def get_related_files_af2(structure_path: Path, binder_name: str, rank: int) -> list[Path]:
    """
    Get all files related to a ColabFold prediction for a specific rank.
    
    Per-rank files:
    - {binder}_unrelaxed_rank_00{N}_*.pdb
    - {binder}_scores_rank_00{N}_*.json
    
    Shared files:
    - {binder}.done.txt
    - {binder}_predicted_aligned_error_v1.json
    """
    parent = structure_path.parent
    related = [structure_path]
    
    # Format rank with leading zeros (e.g., 1 -> "001")
    rank_str = f"{rank:03d}"
    
    # Per-rank score file
    score_pattern = f"{binder_name}_scores_rank_{rank_str}_*.json"
    for f in parent.glob(score_pattern):
        if f not in related:
            related.append(f)
    
    # Shared files
    shared_files = [
        f"{binder_name}.done.txt",
        f"{binder_name}_predicted_aligned_error_v1.json",
    ]
    
    for filename in shared_files:
        f = parent / filename
        if f.exists() and f not in related:
            related.append(f)
    
    return related


def get_related_files_chai(structure_path: Path, model_idx: int) -> list[Path]:
    """
    Get all files related to a Chai-1 prediction for a specific model.
    
    Per-model files:
    - pred.model_idx_{N}.cif
    - scores.model_idx_{N}.npz
    """
    parent = structure_path.parent
    related = [structure_path]
    
    # Per-model score file
    scores_file = parent / f"scores.model_idx_{model_idx}.npz"
    if scores_file.exists() and scores_file not in related:
        related.append(scores_file)
    
    return related


def get_related_files_boltz(structure_path: Path, binder_name: str, model_idx: int) -> list[Path]:
    """
    Get all files related to a Boltz-2 prediction for a specific model.
    
    Per-model files:
    - {binder}_model_{N}.pdb
    - confidence_{binder}_model_{N}.json
    - pae_{binder}_model_{N}.npz
    - pde_{binder}_model_{N}.npz
    - plddt_{binder}_model_{N}.npz
    """
    parent = structure_path.parent
    related = [structure_path]
    
    # Per-model files
    per_model_patterns = [
        f"confidence_{binder_name}_model_{model_idx}.json",
        f"pae_{binder_name}_model_{model_idx}.npz",
        f"pde_{binder_name}_model_{model_idx}.npz",
        f"plddt_{binder_name}_model_{model_idx}.npz",
    ]
    
    for pattern in per_model_patterns:
        f = parent / pattern
        if f.exists() and f not in related:
            related.append(f)
    
    return related


def get_related_files_binder(structure_path: Path, model_idx: int) -> list[Path]:
    """
    Get all files related to a binder-only prediction for a specific model.
    
    Per-model files (Chai-1 format):
    - pred.model_idx_{N}.cif
    - scores.model_idx_{N}.npz
    """
    parent = structure_path.parent
    related = [structure_path]
    
    # Per-model score file
    scores_file = parent / f"scores.model_idx_{model_idx}.npz"
    if scores_file.exists() and scores_file not in related:
        related.append(scores_file)
    
    return related


def get_related_files(structure_path: Path, method: str, binder_name: str, model_rank: int, binder_msa_mode: str | None = None) -> list[Path]:
    """Get all files related to a prediction based on method."""
    if method in ("af2_multimer", "af2_ptm_complex"):
        return get_related_files_af2(structure_path, binder_name, model_rank)
    elif method == "chai":
        return get_related_files_chai(structure_path, model_rank)
    elif method == "boltz":
        return get_related_files_boltz(structure_path, binder_name, model_rank)
    elif method == "binder":
        if binder_msa_mode == "true":
            # AF2-ptm format for binder-only predictions
            return get_related_files_af2(structure_path, binder_name, model_rank)
        else:
            # Chai format for binder-only predictions
            return get_related_files_binder(structure_path, model_rank)
    else:
        return [structure_path]


def copy_files_to_output(files: list[Path], output_dir: Path) -> None:
    """Copy a list of files to output directory, preserving names."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for f in files:
        if f.exists():
            dest = output_dir / f.name
            shutil.copy2(f, dest)


def find_binder_structures(predictions_dir: Path, binder_names: list[str], binder_msa_mode: str | None = None) -> list[dict]:
    """
    Find binder-only structures.
    
    When binder_msa_mode='true' (AF2-ptm format):
        Directory: predictions_dir/binder/ (flat)
        Pattern: {binder}_unrelaxed_rank_{NNN}_*.pdb
        Rank: NNN (1-indexed, rank_001 = best)
    
    When binder_msa_mode='false' (Chai-1 format):
        Directory: predictions_dir/binder/{binder}/
        Pattern: pred.model_idx_{N}.cif
        Rank: N (0-indexed, model_idx_0 = best)
    """
    results = []
    binder_dir = predictions_dir / "binder"
    
    if not binder_dir.exists():
        return results
    
    if binder_msa_mode == "true":
        # AF2-ptm format: flat directory with ranked PDB files
        for binder in binder_names:
            pattern = f"{binder}_unrelaxed_rank_*.pdb"
            pdb_files = list(binder_dir.glob(pattern))
            
            for pdb_file in sorted(pdb_files):
                match = re.search(r'rank_(\d+)', pdb_file.name)
                if match:
                    rank = int(match.group(1))
                    results.append({
                        "binder": binder,
                        "method": "binder",
                        "model_rank": rank,
                        "structure_file": str(pdb_file),
                        "min_distance": None,  # No distance for single-chain
                        "attempt": 1,
                    })
    else:
        # Chai-1 format: per-binder subdirectories with CIF files
        for binder in binder_names:
            binder_subdir = binder_dir / binder
            if not binder_subdir.exists():
                continue
            
            cif_files = list(binder_subdir.glob("pred.model_idx_*.cif"))
            
            for cif_file in sorted(cif_files):
                match = re.search(r'model_idx_(\d+)', cif_file.name)
                if match:
                    model_idx = int(match.group(1))
                    results.append({
                        "binder": binder,
                        "method": "binder",
                        "model_rank": model_idx,
                        "structure_file": str(cif_file),
                        "min_distance": None,  # No distance for single-chain
                        "attempt": 1,
                    })
    
    return results


def load_distances_csv(csv_path: Path) -> list[dict]:
    """Load distances CSV and return list of dicts."""
    rows = []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def select_best_structure(candidates: list[dict], cutoff: float) -> tuple[dict, str]:
    """
    Select the best structure from candidates.
    
    Returns (best_candidate, status) where status is 'passed' or 'failed'.
    """
    # Parse distances
    valid_candidates = []
    for c in candidates:
        try:
            dist = float(c["min_distance"]) if c["min_distance"] else None
        except (ValueError, TypeError):
            dist = None
        
        valid_candidates.append({
            **c,
            "min_distance_float": dist,
        })
    
    # Filter by distance cutoff
    passing = [c for c in valid_candidates 
               if c["min_distance_float"] is not None and c["min_distance_float"] <= cutoff]
    
    if passing:
        # Select best by rank (lowest rank number)
        best = min(passing, key=lambda x: int(x["model_rank"]))
        return best, "passed"
    else:
        # No structures pass - select lowest-ranked as fallback
        ranked = sorted(valid_candidates, key=lambda x: int(x.get("model_rank", 0)))
        if ranked:
            return ranked[0], "failed"
        return None, "no_structures"


def main():
    parser = argparse.ArgumentParser(
        description="Filter structures by distance and select best by rank"
    )
    parser.add_argument("--distances_csv", required=True,
                        help="Input CSV with distance calculations (attempt 1)")
    parser.add_argument("--predictions_dir", required=True,
                        help="Directory containing prediction subdirectories (attempt 1)")
    parser.add_argument("--cutoff", type=float, required=True,
                        help="Distance cutoff in Angstroms")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for filtered structures")
    parser.add_argument("--summary_csv", required=True,
                        help="Output summary CSV file")
    parser.add_argument("--distances_csv_retry", default=None,
                        help="Optional: distances CSV from retry attempt (attempt 2)")
    parser.add_argument("--predictions_dir_retry", default=None,
                        help="Optional: predictions directory from retry attempt")
    parser.add_argument("--binder_msa", default=None, choices=["true", "false"],
                        help="Binder MSA mode: 'true' = AF2-ptm binder outputs, 'false' = Chai binder outputs")
    
    args = parser.parse_args()
    
    distances_csv = Path(args.distances_csv)
    predictions_dir = Path(args.predictions_dir)
    output_dir = Path(args.output_dir)
    summary_csv = Path(args.summary_csv)
    cutoff = args.cutoff
    
    # Load attempt 1 distances
    rows_attempt1 = load_distances_csv(distances_csv)
    print(f"Loaded {len(rows_attempt1)} rows from attempt 1: {distances_csv}")
    
    # Load attempt 2 distances if provided
    rows_attempt2 = []
    predictions_dir_retry = None
    if args.distances_csv_retry and Path(args.distances_csv_retry).exists():
        rows_attempt2 = load_distances_csv(Path(args.distances_csv_retry))
        print(f"Loaded {len(rows_attempt2)} rows from attempt 2: {args.distances_csv_retry}")
        if args.predictions_dir_retry:
            predictions_dir_retry = Path(args.predictions_dir_retry)
    
    print(f"Distance cutoff: {cutoff} Å")
    
    # Get binder names from attempt 1
    binder_names = list(set(row["binder"] for row in rows_attempt1))
    
    # Determine binder MSA mode
    binder_msa_mode = args.binder_msa
    
    # Find binder-only structures (not in distances CSV)
    binder_structures = find_binder_structures(predictions_dir, binder_names, binder_msa_mode)
    print(f"Found {len(binder_structures)} binder-only structures")
    
    # Group attempt 1 by (binder, method)
    groups_attempt1 = defaultdict(list)
    for row in rows_attempt1:
        key = (row["binder"], row["method"])
        groups_attempt1[key].append(row)
    
    # Group attempt 2 by (binder, method)
    groups_attempt2 = defaultdict(list)
    for row in rows_attempt2:
        key = (row["binder"], row["method"])
        groups_attempt2[key].append(row)
    
    # Group binder structures by binder
    binder_groups = defaultdict(list)
    for struct in binder_structures:
        binder_groups[struct["binder"]].append(struct)
    
    print(f"Found {len(groups_attempt1)} (binder, method) combinations")
    
    # Process each group
    summary_rows = []
    method_stats = defaultdict(lambda: {"passed": 0, "failed": 0, "total": 0, "retry_helped": 0})
    
    for (binder, method), candidates_a1 in sorted(groups_attempt1.items()):
        method_stats[method]["total"] += 1
        
        # Try attempt 1 first
        best_a1, status_a1 = select_best_structure(candidates_a1, cutoff)
        
        # Determine which attempt to use
        use_attempt = 1
        best = best_a1
        status = status_a1
        pred_dir = predictions_dir
        
        # If attempt 1 failed and attempt 2 exists, try attempt 2
        if status_a1 == "failed" and (binder, method) in groups_attempt2:
            candidates_a2 = groups_attempt2[(binder, method)]
            best_a2, status_a2 = select_best_structure(candidates_a2, cutoff)
            
            if status_a2 == "passed":
                # Attempt 2 helped!
                use_attempt = 2
                best = best_a2
                status = status_a2
                pred_dir = predictions_dir_retry
                method_stats[method]["retry_helped"] += 1
            elif best_a2 is not None:
                # Both failed, but use attempt 2's best (freshest)
                use_attempt = 2
                best = best_a2
                status = "failed"
                pred_dir = predictions_dir_retry
        
        if best is not None:
            # Get related files and copy
            structure_path = Path(best["structure_file"])
            model_rank = int(best["model_rank"])
            related_files = get_related_files(structure_path, method, binder, model_rank, binder_msa_mode)
            
            # Create output subdirectory: filtered_structures/{method}/{binder}
            subdir = output_dir / method / binder
            copy_files_to_output(related_files, subdir)
            
            summary_rows.append({
                "binder": binder,
                "method": method,
                "selected_structure": str(subdir / structure_path.name),
                "min_distance": best["min_distance"],
                "model_rank": best["model_rank"],
                "num_candidates": len(candidates_a1),
                "num_passing": len([c for c in candidates_a1 
                                   if c.get("min_distance") and 
                                   float(c["min_distance"]) <= cutoff]) if status_a1 == "passed" else 0,
                "status": status,
                "attempt": use_attempt,
            })
            
            if status == "passed":
                method_stats[method]["passed"] += 1
                print(f"  {binder}/{method}: PASSED (dist={best['min_distance']}, rank={best['model_rank']}, attempt={use_attempt})")
            else:
                method_stats[method]["failed"] += 1
                dist_str = f"{best['min_distance']}" if best.get("min_distance") else "N/A"
                print(f"  {binder}/{method}: FAILED (dist={dist_str}, rank={best['model_rank']}, attempt={use_attempt})")
        else:
            summary_rows.append({
                "binder": binder,
                "method": method,
                "selected_structure": "",
                "min_distance": "",
                "model_rank": "",
                "num_candidates": 0,
                "num_passing": 0,
                "status": "no_structures",
                "attempt": 1,
            })
            method_stats[method]["failed"] += 1
            print(f"  {binder}/{method}: NO STRUCTURES")
    
    # Process binder-only structures (no distance filtering, just copy best-ranked model)
    for binder, candidates in sorted(binder_groups.items()):
        # Select the best-ranked model (lowest rank number)
        # AF2-ptm: rank_001 = best (1-indexed), Chai: model_idx_0 = best (0-indexed)
        best = min(candidates, key=lambda c: c["model_rank"])
        
        structure_path = Path(best["structure_file"])
        model_rank = best["model_rank"]
        related_files = get_related_files(structure_path, "binder", binder, model_rank, binder_msa_mode)
        
        # Create output subdirectory
        subdir = output_dir / "binder" / binder
        copy_files_to_output(related_files, subdir)
        
        summary_rows.append({
            "binder": binder,
            "method": "binder",
            "selected_structure": str(subdir / structure_path.name),
            "min_distance": "N/A",
            "model_rank": model_rank,
            "num_candidates": len(candidates),
            "num_passing": "N/A",
            "status": "copied",
            "attempt": 1,
        })
        print(f"  {binder}/binder: COPIED (rank={model_rank})")
    
    # Write summary CSV
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    
    fieldnames = [
        "binder", "method", "selected_structure", "min_distance", "model_rank",
        "num_candidates", "num_passing", "status", "attempt"
    ]
    
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)
    
    # Print summary statistics
    print(f"\n{'='*70}")
    print("Summary Statistics")
    print(f"{'='*70}")
    print(f"Distance cutoff: {cutoff} Å")
    print(f"\n{'Method':<20} {'Passed':<10} {'Failed':<10} {'Total':<10} {'Retry Helped':<12}")
    print("-" * 70)
    
    total_passed = 0
    total_failed = 0
    total_count = 0
    total_retry_helped = 0
    
    for method in sorted(method_stats.keys()):
        stats = method_stats[method]
        passed = stats["passed"]
        failed = stats["failed"]
        total = stats["total"]
        retry_helped = stats["retry_helped"]
        total_passed += passed
        total_failed += failed
        total_count += total
        total_retry_helped += retry_helped
        print(f"{method:<20} {passed:<10} {failed:<10} {total:<10} {retry_helped:<12}")
    
    # Add binder stats (not filtered)
    binder_count = len(binder_groups)
    if binder_count > 0:
        print(f"{'binder':<20} {'N/A':<10} {'N/A':<10} {binder_count:<10} {'N/A':<12}")
    
    print("-" * 70)
    print(f"{'Total (filtered)':<20} {total_passed:<10} {total_failed:<10} {total_count:<10} {total_retry_helped:<12}")
    
    if total_count > 0:
        pass_rate = 100 * total_passed / total_count
        print(f"\nPass rate: {total_passed}/{total_count} ({pass_rate:.1f}%)")
    
    if total_retry_helped > 0:
        print(f"Retry improved: {total_retry_helped} binder-method pairs")
    
    print(f"\nOutput directory: {output_dir}")
    print(f"Summary CSV: {summary_csv}")
    
    return 0


if __name__ == "__main__":
    exit(main())
