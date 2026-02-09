#!/usr/bin/env python3
"""
generate_multimer_a3m.py

Convert target-only MSAs to target:binder multimer MSAs by appending binder
sequences and gaps according to ColabFold multimer format.

This script reads:
1. Target MSA files (e.g., AQP.a3m) from targets/msa/
2. Binder sequences from input CSV
3. Generates multimer MSAs (e.g., AQP_1.a3m, AQP_2.a3m) with proper formatting

Format transformation:
  Target MSA:
    #242	1
    >101
    TARGETSEQ...

  Multimer MSA:
    #242,100	1,1
    >101	102
    TARGETSEQ...BINDERSEQ...
    >101
    TARGETSEQ...[100 gaps]
    [MSA sequences with gaps appended]
    >102
    [242 gaps]BINDERSEQ...

Usage:
    python generate_multimer_a3m.py \
        --input_csv input.csv \
        --target_msa_dir targets/msa \
        --output_dir project/colabfold/msa_multimer
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def parse_target_a3m(a3m_path: Path) -> Tuple[int, str, List[Tuple[str, str]]]:
    """
    Parse target MSA file.
    
    Returns:
        (target_length, target_sequence, msa_entries)
        
        msa_entries: List of (header, sequence) tuples
    """
    with open(a3m_path) as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        raise ValueError(f"Invalid a3m file: {a3m_path} (too few lines)")
    
    # Line 1: #242	1
    header_line = lines[0].strip()
    if not header_line.startswith('#'):
        raise ValueError(f"Invalid header line: {header_line}")
    
    target_length = int(header_line.split()[0][1:])  # Remove '#' and parse
    
    # Line 2: >101
    chain_header = lines[1].strip()
    
    # Line 3: Target sequence
    target_sequence = lines[2].strip()
    
    if len(target_sequence) != target_length:
        raise ValueError(
            f"Target sequence length mismatch: expected {target_length}, "
            f"got {len(target_sequence)}"
        )
    
    # Remaining lines: MSA entries (header + sequence pairs)
    msa_entries = []
    i = 3
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        
        if line.startswith('>'):
            # Header line
            header = line
            # Next line should be sequence
            if i + 1 < len(lines):
                sequence = lines[i + 1].strip()
                msa_entries.append((header, sequence))
                i += 2
            else:
                i += 1
        else:
            i += 1
    
    return target_length, target_sequence, msa_entries


def generate_multimer_a3m(
    target_a3m_path: Path,
    binder_name: str,
    binder_sequence: str,
    output_path: Path
) -> None:
    """
    Generate multimer a3m file from target MSA and binder sequence.
    
    Args:
        target_a3m_path: Path to target MSA file
        binder_name: Name for this binder (e.g., "AQP_1")
        binder_sequence: Binder amino acid sequence
        output_path: Output multimer a3m file path
    """
    # Parse target MSA
    target_len, target_seq, msa_entries = parse_target_a3m(target_a3m_path)
    binder_len = len(binder_sequence)
    
    # Create output
    with open(output_path, 'w') as f:
        # Line 1: #TARGET_LEN,BINDER_LEN	1,1
        f.write(f"#{target_len},{binder_len}\t1,1\n")
        
        # Line 2: >101	102
        f.write(">101\t102\n")
        
        # Line 3: TARGET_SEQ + BINDER_SEQ (no separator)
        f.write(f"{target_seq}{binder_sequence}\n")
        
        # Line 4: >101
        f.write(">101\n")
        
        # Line 5: TARGET_SEQ + [BINDER_LEN gaps]
        gaps = "-" * binder_len
        f.write(f"{target_seq}{gaps}\n")
        
        # MSA entries: each sequence + [BINDER_LEN gaps]
        for header, sequence in msa_entries:
            f.write(f"{header}\n")
            f.write(f"{sequence}{gaps}\n")
        
        # Second-to-last line: >102
        f.write(">102\n")
        
        # Last line: [TARGET_LEN gaps] + BINDER_SEQ
        target_gaps = "-" * target_len
        f.write(f"{target_gaps}{binder_sequence}\n")


def extract_target_name(full_name: str) -> str:
    """
    Extract target name from binder name.
    
    Examples:
        AQP_1 -> AQP
        SLP_23 -> SLP
        CLN5_5 -> CLN5
    """
    # Split by underscore and take all but last part
    parts = full_name.split('_')
    if len(parts) < 2:
        raise ValueError(f"Invalid binder name format: {full_name} (expected NAME_NUMBER)")
    
    return '_'.join(parts[:-1])


def process_csv(
    csv_path: Path,
    target_msa_dir: Path,
    output_dir: Path
) -> Dict[str, int]:
    """
    Process CSV and generate multimer a3m files.
    
    Args:
        csv_path: Input CSV with name,target,binder columns
        target_msa_dir: Directory containing target MSA files
        output_dir: Output directory for multimer MSAs
    
    Returns:
        Statistics dict
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    stats = {
        'total': 0,
        'success': 0,
        'failed': 0,
        'targets_used': set()
    }
    
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            stats['total'] += 1
            
            binder_name = row['name'].strip()
            target_seq = row['target'].strip()
            binder_seq = row['binder'].strip()
            
            # Extract target name
            try:
                target_name = extract_target_name(binder_name)
            except ValueError as e:
                print(f"ERROR: {e}", file=sys.stderr)
                stats['failed'] += 1
                continue
            
            # Find target MSA file
            target_a3m = target_msa_dir / f"{target_name}.a3m"
            if not target_a3m.exists():
                print(
                    f"ERROR: Target MSA not found: {target_a3m}",
                    file=sys.stderr
                )
                stats['failed'] += 1
                continue
            
            # Generate multimer MSA
            output_path = output_dir / f"{binder_name}.a3m"
            
            try:
                generate_multimer_a3m(
                    target_a3m,
                    binder_name,
                    binder_seq,
                    output_path
                )
                
                stats['success'] += 1
                stats['targets_used'].add(target_name)
                
                print(f"✓ Generated: {binder_name}.a3m", file=sys.stderr)
                
            except Exception as e:
                print(
                    f"ERROR generating {binder_name}: {e}",
                    file=sys.stderr
                )
                stats['failed'] += 1
    
    return stats


def validate_multimer_a3m(a3m_path: Path) -> bool:
    """
    Validate generated multimer a3m file.
    
    Checks:
    - Header format
    - Sequence lengths
    - Gap patterns
    """
    try:
        with open(a3m_path) as f:
            lines = f.readlines()
        
        if len(lines) < 6:
            print(f"ERROR: Too few lines in {a3m_path}", file=sys.stderr)
            return False
        
        # Parse header
        header = lines[0].strip()
        if not header.startswith('#'):
            print(f"ERROR: Invalid header: {header}", file=sys.stderr)
            return False
        
        lengths = header[1:].split('\t')[0]
        target_len, binder_len = map(int, lengths.split(','))
        expected_len = target_len + binder_len
        
        # Check line 3 (full sequence)
        full_seq = lines[2].strip()
        if len(full_seq) != expected_len:
            print(
                f"ERROR: Line 3 length mismatch: expected {expected_len}, "
                f"got {len(full_seq)}",
                file=sys.stderr
            )
            return False
        
        # Check line 5 (target + gaps)
        target_with_gaps = lines[4].strip()
        if len(target_with_gaps) != expected_len:
            print(
                f"ERROR: Line 5 length mismatch: expected {expected_len}, "
                f"got {len(target_with_gaps)}",
                file=sys.stderr
            )
            return False
        
        # Check last line (gaps + binder)
        last_line = lines[-1].strip()
        if len(last_line) != expected_len:
            print(
                f"ERROR: Last line length mismatch: expected {expected_len}, "
                f"got {len(last_line)}",
                file=sys.stderr
            )
            return False
        
        return True
        
    except Exception as e:
        print(f"ERROR validating {a3m_path}: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Generate multimer a3m files from target MSAs and binder sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python generate_multimer_a3m.py \\
        --input_csv input.csv \\
        --target_msa_dir targets/msa \\
        --output_dir project/colabfold/msa_multimer

Input CSV format:
    name,target,binder
    AQP_1,MSEQ_TARGET,MSEQ_BINDER1
    AQP_2,MSEQ_TARGET,MSEQ_BINDER2
    SLP_1,MSEQ_TARGET,MSEQ_BINDER3

Requirements:
    - Target MSA files must exist: targets/msa/AQP.a3m, targets/msa/SLP.a3m
    - Binder names must follow NAME_NUMBER format (e.g., AQP_1, SLP_23)
        """
    )
    
    parser.add_argument(
        '--input_csv',
        type=Path,
        required=True,
        help='Input CSV with name,target,binder columns'
    )
    
    parser.add_argument(
        '--target_msa_dir',
        type=Path,
        required=True,
        help='Directory containing target MSA files (e.g., targets/msa/)'
    )
    
    parser.add_argument(
        '--output_dir',
        type=Path,
        required=True,
        help='Output directory for multimer MSA files'
    )
    
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate generated files'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.input_csv.exists():
        print(f"ERROR: Input CSV not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)
    
    if not args.target_msa_dir.exists():
        print(f"ERROR: Target MSA directory not found: {args.target_msa_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Process CSV
    print("Generating multimer a3m files...", file=sys.stderr)
    stats = process_csv(args.input_csv, args.target_msa_dir, args.output_dir)
    
    # Validate if requested
    if args.validate and stats['success'] > 0:
        print("\nValidating generated files...", file=sys.stderr)
        for a3m_file in args.output_dir.glob("*.a3m"):
            if not validate_multimer_a3m(a3m_file):
                print(f"✗ Validation failed: {a3m_file.name}", file=sys.stderr)
            else:
                print(f"✓ Validated: {a3m_file.name}", file=sys.stderr)
    
    # Print summary
    print("\n" + "="*60, file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print("="*60, file=sys.stderr)
    print(f"Total entries: {stats['total']}", file=sys.stderr)
    print(f"Successful: {stats['success']}", file=sys.stderr)
    print(f"Failed: {stats['failed']}", file=sys.stderr)
    print(f"Unique targets used: {len(stats['targets_used'])}", file=sys.stderr)
    print(f"Targets: {', '.join(sorted(stats['targets_used']))}", file=sys.stderr)
    print(f"Output directory: {args.output_dir}", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    sys.exit(0 if stats['failed'] == 0 else 1)


if __name__ == '__main__':
    main()