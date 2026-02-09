#!/usr/bin/env python3
"""
generate_boltz_yaml.py - Generate Boltz-2 YAML input files from CSV

Creates one YAML file per binder with:
- Target protein (chain A) with MSA from target_msa_dir
- Binder protein (chain B) with MSA from binder_msa_dir (or 'empty' if binder_msa_mode=false)
- Template from target_pdb_dir assigned to chain A

Usage:
    python generate_boltz_yaml.py \\
        --csv input.csv \\
        --target_msa_dir /path/to/target/msa \\
        --target_pdb_dir /path/to/target/pdb \\
        --binder_msa_dir /path/to/binder/msa \\
        --binder_msa_mode true|false \\
        --output_dir /path/to/output
"""

import argparse
import csv
import re
from pathlib import Path


def extract_target_name(binder_name: str) -> str:
    """
    Extract target name from binder name.
    
    Binder names must follow the format TARGET_NUMBER (e.g., AQP_01, AQP2_01, IL6R_123).
    The target name is everything before the last underscore, and the suffix must be
    purely numeric.
    
    # NOTE: keep in sync with calculate_distances.py
    
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


def find_file_case_insensitive(directory: Path, basename: str, extensions: list) -> Path | None:
    """Find a file with case-insensitive matching."""
    for ext in extensions:
        # Try exact match first
        exact = directory / f"{basename}{ext}"
        if exact.exists():
            return exact
        
        # Try case-insensitive
        for f in directory.iterdir():
            if f.name.lower() == f"{basename.lower()}{ext.lower()}":
                return f
    
    return None


def generate_yaml(
    name: str,
    target_seq: str,
    binder_seq: str,
    target_msa_path: str,
    binder_msa_path: str,
    template_pdb_path: str | None,
) -> str:
    """Generate YAML content for a single complex."""
    
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {target_seq}",
        f"      msa: {target_msa_path}",
        "  - protein:",
        "      id: B",
        f"      sequence: {binder_seq}",
        f"      msa: {binder_msa_path}",
    ]
    
    # Add template section if template exists
    if template_pdb_path:
        lines.extend([
            "",
            "templates:",
            f"  - pdb: {template_pdb_path}",
            "    chain_id: A",
        ])
    
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Generate Boltz-2 YAML input files from CSV"
    )
    parser.add_argument("--csv", required=True, help="Input CSV file (name,target,binder)")
    parser.add_argument("--target_msa_dir", required=True, help="Directory with target MSAs")
    parser.add_argument("--target_pdb_dir", required=True, help="Directory with target PDBs")
    parser.add_argument("--binder_msa_dir", required=True, help="Directory with binder MSAs")
    parser.add_argument("--binder_msa_mode", required=True, choices=["true", "false"],
                        help="Whether binder MSAs exist (true) or should use empty (false)")
    parser.add_argument("--output_dir", required=True, help="Output directory for YAML files")
    
    args = parser.parse_args()
    
    csv_path = Path(args.csv)
    target_msa_dir = Path(args.target_msa_dir)
    target_pdb_dir = Path(args.target_pdb_dir)
    binder_msa_dir = Path(args.binder_msa_dir)
    output_dir = Path(args.output_dir)
    use_binder_msa = args.binder_msa_mode == "true"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read CSV
    with open(csv_path, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        print("ERROR: No data rows in CSV")
        return 1
    
    # Validate headers
    required_cols = {'name', 'target', 'binder'}
    headers = set(rows[0].keys())
    missing = required_cols - headers
    if missing:
        print(f"ERROR: Missing columns in CSV: {missing}")
        return 1
    
    created = 0
    errors = 0
    
    for row in rows:
        name = row['name'].strip()
        target_seq = row['target'].strip()
        binder_seq = row['binder'].strip()
        
        if not name or not target_seq or not binder_seq:
            print(f"WARNING: Skipping row with missing data: {row}")
            continue
        
        # Extract target name from binder name
        target_name = extract_target_name(name)
        
        # Find target MSA
        target_msa = find_file_case_insensitive(target_msa_dir, target_name, ['.a3m', '.a3m.gz'])
        if not target_msa:
            print(f"ERROR: Target MSA not found for {target_name} in {target_msa_dir}")
            errors += 1
            continue
        
        # Find target PDB template
        target_pdb = find_file_case_insensitive(target_pdb_dir, target_name, ['.pdb', '.cif'])
        if not target_pdb:
            print(f"WARNING: Target PDB not found for {target_name}, proceeding without template")
        
        # Determine binder MSA path
        if use_binder_msa:
            binder_msa = find_file_case_insensitive(binder_msa_dir, name, ['.a3m', '.a3m.gz'])
            if not binder_msa:
                print(f"WARNING: Binder MSA not found for {name}, using 'empty'")
                binder_msa_path = "empty"
            else:
                binder_msa_path = str(binder_msa.resolve())
        else:
            binder_msa_path = "empty"
        
        # Generate YAML
        yaml_content = generate_yaml(
            name=name,
            target_seq=target_seq,
            binder_seq=binder_seq,
            target_msa_path=str(target_msa.resolve()),
            binder_msa_path=binder_msa_path,
            template_pdb_path=str(target_pdb.resolve()) if target_pdb else None,
        )
        
        # Write YAML file
        yaml_path = output_dir / f"{name}.yaml"
        with open(yaml_path, 'w') as f:
            f.write(yaml_content)
        
        created += 1
    
    print(f"Created {created} YAML files in {output_dir}")
    if errors:
        print(f"Encountered {errors} errors")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
