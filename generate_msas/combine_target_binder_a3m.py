#!/usr/bin/env python3
"""
Combine individual target and binder A3M MSAs into ColabFold multimer A3M format.

For each row in the input CSV (name, target_name, binder_seq), this script:
  1. Reads the target MSA from target_msa_dir/{target_name}.a3m
  2. Reads the binder MSA from binder_msa_dir/{name}.a3m
  3. Produces a multimer A3M at output_dir/{name}.a3m

Output format (ColabFold unpaired multimer A3M):
  Line 1: #<len_target>,<len_binder>\t1,1
  Line 2: >101\t102
  Line 3: <target_query><binder_query>
  >101
  <target_query><'-' * len_binder>
  <target_hit_1_header>
  <target_hit_1_seq><'-' * len_binder>
  ...
  >102
  <'-' * len_target><binder_query>
  <binder_hit_1_header>
  <'-' * len_target><binder_hit_1_seq>
  ...
"""

import argparse
import csv
import os
import sys


def parse_a3m(filepath, chain_id=None):
    """Parse an A3M file, returning (query_seq, [(header, seq), ...]).

    The first entry is the query sequence. Subsequent entries are MSA hits.
    The A3M header line (starting with #) is skipped.

    ColabFold server a3m files contain two >101 blocks: the first for UniRef
    hits, the second for environmental DB hits. When chain_id is provided,
    any secondary >101 headers (environmental DB separators) are renamed to
    >chain_id so they carry the correct chain identifier in the combined
    multimer a3m output.
    """
    headers = []
    sequences = []
    current_header = None
    current_seq_parts = []
    seen_101 = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append("".join(current_seq_parts))
                current_header = line
                # Rename secondary >101 to >chain_id
                if current_header == ">101":
                    if seen_101:
                        if chain_id is not None:
                            current_header = f">{chain_id}"
                    else:
                        seen_101 = True
                headers.append(current_header)
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

    # Flush last entry
    if current_header is not None:
        sequences.append("".join(current_seq_parts))

    if len(headers) == 0:
        raise ValueError(f"No entries found in {filepath}")

    query_seq = sequences[0]
    entries = list(zip(headers[1:], sequences[1:]))

    return query_seq, entries


def compute_query_length(query_seq):
    """Compute the aligned length of the query.

    In A3M format, the query is the first sequence and defines the alignment
    columns. Its length is simply len(query_seq) since the query itself has
    no lowercase insertion characters.
    """
    return len(query_seq)


def combine_a3m(target_a3m_path, binder_a3m_path, output_path):
    """Combine a target A3M and binder A3M into a multimer A3M."""
    target_query, target_entries = parse_a3m(target_a3m_path, chain_id="101")
    binder_query, binder_entries = parse_a3m(binder_a3m_path, chain_id="102")

    len_target = compute_query_length(target_query)
    len_binder = compute_query_length(binder_query)

    target_gap = "-" * len_binder
    binder_gap = "-" * len_target

    with open(output_path, "w") as out:
        # Header
        out.write(f"#{len_target},{len_binder}\t1,1\n")
        out.write(f">101\t102\n")
        out.write(f"{target_query}{binder_query}\n")

        # >101 block: target unpaired MSA
        out.write(">101\n")
        out.write(f"{target_query}{target_gap}\n")
        for header, seq in target_entries:
            out.write(f"{header}\n")
            out.write(f"{seq}{target_gap}\n")

        # >102 block: binder unpaired MSA
        out.write(">102\n")
        out.write(f"{binder_gap}{binder_query}\n")
        for header, seq in binder_entries:
            out.write(f"{header}\n")
            out.write(f"{binder_gap}{seq}\n")

    return len_target, len_binder, len(target_entries), len(binder_entries)


def main():
    parser = argparse.ArgumentParser(
        description="Combine target and binder A3M MSAs into ColabFold multimer A3M format."
    )
    parser.add_argument(
        "--input_csv",
        required=True,
        help="CSV file with columns: name, target, binder (header row required)",
    )
    parser.add_argument(
        "--target_msa_dir",
        required=True,
        help="Directory containing target MSAs ({target_name}.a3m)",
    )
    parser.add_argument(
        "--binder_msa_dir",
        required=True,
        help="Directory containing binder MSAs ({name}.a3m)",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for multimer A3M files",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate output files after generation",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Read all available target MSAs and index by query sequence
    target_msa_files = [
        f for f in os.listdir(args.target_msa_dir) if f.endswith(".a3m")
    ]
    target_seq_to_name = {}
    for msa_file in target_msa_files:
        target_name = msa_file.replace(".a3m", "")
        filepath = os.path.join(args.target_msa_dir, msa_file)
        try:
            query_seq, _ = parse_a3m(filepath)
            target_seq_to_name[query_seq] = target_name
        except Exception as e:
            print(f"[WARN] Could not parse {filepath}: {e}", file=sys.stderr)

    if not target_seq_to_name:
        print("[ERROR] No valid target MSA files found", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loaded {len(target_seq_to_name)} target MSAs")

    # Process CSV
    n_success = 0
    n_fail = 0
    errors = []

    with open(args.input_csv, "r") as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # skip header

        for row in reader:
            if len(row) < 3 or not row[0].strip():
                continue

            name = row[0].strip()
            target_seq = row[1].strip()

            # Find matching target MSA
            target_name = target_seq_to_name.get(target_seq)
            if target_name is None:
                msg = f"No target MSA found for '{name}' (sequence not matched)"
                print(f"[ERROR] {msg}", file=sys.stderr)
                errors.append(msg)
                n_fail += 1
                continue

            target_a3m = os.path.join(args.target_msa_dir, f"{target_name}.a3m")
            binder_a3m = os.path.join(args.binder_msa_dir, f"{name}.a3m")
            output_a3m = os.path.join(args.output_dir, f"{name}.a3m")

            if not os.path.isfile(target_a3m):
                msg = f"Target MSA not found: {target_a3m}"
                print(f"[ERROR] {msg}", file=sys.stderr)
                errors.append(msg)
                n_fail += 1
                continue

            if not os.path.isfile(binder_a3m):
                msg = f"Binder MSA not found: {binder_a3m}"
                print(f"[ERROR] {msg}", file=sys.stderr)
                errors.append(msg)
                n_fail += 1
                continue

            try:
                lt, lb, nt, nb = combine_a3m(target_a3m, binder_a3m, output_a3m)
                print(
                    f"[INFO] {name}: target_len={lt}, binder_len={lb}, "
                    f"target_hits={nt}, binder_hits={nb}"
                )
                n_success += 1
            except Exception as e:
                msg = f"Failed to combine MSAs for '{name}': {e}"
                print(f"[ERROR] {msg}", file=sys.stderr)
                errors.append(msg)
                n_fail += 1

    print(f"\n[INFO] Results: {n_success} succeeded, {n_fail} failed")

    if args.validate and n_success > 0:
        print("[INFO] Validating output files...")
        n_valid = 0
        for f in os.listdir(args.output_dir):
            if not f.endswith(".a3m"):
                continue
            filepath = os.path.join(args.output_dir, f)
            try:
                with open(filepath, "r") as fh:
                    first_line = fh.readline().strip()
                    if not first_line.startswith("#") or "\t" not in first_line:
                        print(f"[WARN] Invalid header in {f}: {first_line}")
                        continue
                    second_line = fh.readline().strip()
                    if not second_line.startswith(">101\t102"):
                        print(f"[WARN] Invalid chain header in {f}: {second_line}")
                        continue
                n_valid += 1
            except Exception as e:
                print(f"[WARN] Validation error for {f}: {e}")
        print(f"[INFO] Validation: {n_valid}/{n_success} files valid")

    if n_fail > 0:
        print(f"\n[ERROR] {n_fail} failures occurred:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()