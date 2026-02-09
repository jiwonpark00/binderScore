#!/usr/bin/env python3
"""
prepare_retry_inputs.py - Prepare inputs for retrying failed structure predictions

Reads the summary CSV from initial filtering and prepares input files for
only the failed (binder, method) pairs.

For each method:
- af2_multimer/af2_ptm_complex: Copy .a3m files from colabfold/msa_multimer/
- chai: Copy .fasta files from chai/fasta_multimer/ (MSA dir is reused - hash-based)
- boltz: Copy .yaml files from boltz/yaml/

Usage:
    python prepare_retry_inputs.py \\
        --summary_csv /path/to/summary_attempt1.csv \\
        --project_dir /path/to/project \\
        --retry_inputs_dir /path/to/retry_inputs

Outputs:
    retry_inputs/
    ├── colabfold/
    │   └── msa_multimer/     # .a3m files for failed af2_* binders
    ├── chai/
    │   └── fasta_multimer/   # .fasta files for failed chai binders
    └── boltz/
        └── yaml/             # .yaml files for failed boltz binders
"""

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Prepare inputs for retrying failed structure predictions"
    )
    parser.add_argument("--summary_csv", required=True,
                        help="Summary CSV from initial filtering (with status column)")
    parser.add_argument("--project_dir", required=True,
                        help="Project directory containing original inputs")
    parser.add_argument("--retry_inputs_dir", required=True,
                        help="Output directory for retry inputs")
    
    args = parser.parse_args()
    
    summary_csv = Path(args.summary_csv)
    project_dir = Path(args.project_dir)
    retry_inputs_dir = Path(args.retry_inputs_dir)
    
    # Read summary CSV and find failed pairs
    failed_by_method = defaultdict(set)
    
    with open(summary_csv, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip binder-only (no distance filtering)
            if row["method"] == "binder":
                continue
            # Check for failed status
            if row["status"] == "failed":
                failed_by_method[row["method"]].add(row["binder"])
    
    # Print summary of failures
    total_failures = sum(len(binders) for binders in failed_by_method.values())
    print(f"Found {total_failures} failed (binder, method) pairs:")
    for method, binders in sorted(failed_by_method.items()):
        print(f"  {method}: {len(binders)} binders")
        for b in sorted(binders):
            print(f"    - {b}")
    
    if total_failures == 0:
        print("\nNo failures to retry. Exiting.")
        return 0
    
    # Create retry inputs directory structure
    retry_inputs_dir.mkdir(parents=True, exist_ok=True)
    
    # Process AF2 methods (af2_multimer, af2_ptm_complex)
    af2_methods = {"af2_multimer", "af2_ptm_complex"}
    af2_failed_binders = set()
    for method in af2_methods:
        af2_failed_binders.update(failed_by_method.get(method, set()))
    
    if af2_failed_binders:
        # Source: colabfold/msa_multimer/
        src_msa_dir = project_dir / "colabfold" / "msa_multimer"
        dst_msa_dir = retry_inputs_dir / "colabfold" / "msa_multimer"
        dst_msa_dir.mkdir(parents=True, exist_ok=True)
        
        copied = 0
        for binder in af2_failed_binders:
            src_file = src_msa_dir / f"{binder}.a3m"
            if src_file.exists():
                shutil.copy2(src_file, dst_msa_dir / src_file.name)
                copied += 1
            else:
                print(f"  WARNING: MSA file not found: {src_file}")
        
        print(f"\nCopied {copied} .a3m files for AF2 retry")
    
    # Process Chai method
    chai_failed_binders = failed_by_method.get("chai", set())
    if chai_failed_binders:
        # Source: chai/fasta_multimer/
        src_fasta_dir = project_dir / "chai" / "fasta_multimer"
        dst_fasta_dir = retry_inputs_dir / "chai" / "fasta_multimer"
        dst_fasta_dir.mkdir(parents=True, exist_ok=True)
        
        copied = 0
        for binder in chai_failed_binders:
            src_file = src_fasta_dir / f"{binder}.fasta"
            if src_file.exists():
                shutil.copy2(src_file, dst_fasta_dir / src_file.name)
                copied += 1
            else:
                print(f"  WARNING: FASTA file not found: {src_file}")
        
        print(f"Copied {copied} .fasta files for Chai retry")
        print("  Note: Chai MSA directory will be reused (hash-based lookup)")
    
    # Process Boltz method
    boltz_failed_binders = failed_by_method.get("boltz", set())
    if boltz_failed_binders:
        # Source: boltz/yaml/
        src_yaml_dir = project_dir / "boltz" / "yaml"
        dst_yaml_dir = retry_inputs_dir / "boltz" / "yaml"
        dst_yaml_dir.mkdir(parents=True, exist_ok=True)
        
        copied = 0
        for binder in boltz_failed_binders:
            src_file = src_yaml_dir / f"{binder}.yaml"
            if src_file.exists():
                shutil.copy2(src_file, dst_yaml_dir / src_file.name)
                copied += 1
            else:
                print(f"  WARNING: YAML file not found: {src_file}")
        
        print(f"Copied {copied} .yaml files for Boltz retry")
    
    print(f"\nRetry inputs prepared in: {retry_inputs_dir}")
    
    # Write a manifest of what needs to be retried
    manifest_file = retry_inputs_dir / "retry_manifest.csv"
    with open(manifest_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["binder", "method"])
        for method, binders in sorted(failed_by_method.items()):
            for binder in sorted(binders):
                writer.writerow([binder, method])
    
    print(f"Retry manifest written to: {manifest_file}")
    
    return 0


if __name__ == "__main__":
    exit(main())
