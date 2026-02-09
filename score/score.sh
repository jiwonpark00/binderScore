#!/usr/bin/env bash
#===============================================================================
# Stage 3: Metrics Extraction and Scoring for binder screening pipeline
#
# Inputs (named flags):
#   --project              - Project directory (from predict_structure.sh)
#
# Optional flags:
#   --pae-threshold        (default: 10.0)
#   --pde-threshold        (default: 2.0)
#   --interface-threshold  (default: 4.0)
#   --pdockq-threshold     (default: 8.0)
#
# Outputs:
#   - $project/metrics/{af2_multimer,af2_ptm_complex,boltz,chai,binder}_metrics.csv
#   - $project/metrics/binderScore.csv (final merged, one row per binder)
#   - $project/symlinks/{method}/*.pdb (symlinks to structure files)
#   - $project/metrics/{method}_rosettaScore.csv (Rosetta interface metrics)
#   - Completion marker: $project/.scoring_complete
#===============================================================================

set -euo pipefail
export PYTHONWARNINGS="ignore"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#===============================================================================
# Logging functions
#===============================================================================

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') $*"
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

log_stage() {
    echo ""
    echo "==================================================="
    echo " $*"
    echo "==================================================="
    echo ""
}

#===============================================================================
# Help message
#===============================================================================

show_help() {
    cat <<EOF
Usage: $0 [OPTIONS]

Stage 3: Metrics Extraction and Scoring for binder screening pipeline

Required flags:
  --project              Project directory (from predict_structure.sh)

Optional flags:
  --pae-threshold        PAE threshold (default: 10.0)
  --pde-threshold        PDE threshold for Boltz (default: 2.0)
  --interface-threshold  Interface distance threshold (default: 4.0)
  --pdockq-threshold     pDockQ threshold (default: 8.0)

Options:
  --help, -h             Show this help message and exit

Examples:
  $0 --project /path/to/project
  $0 --project /path/to/project --pae-threshold 8.0 --pdockq-threshold 6.0
EOF
}

#===============================================================================
# Argument parsing
#===============================================================================

project=""
pae_threshold="10.0"
pde_threshold="2.0"
interface_threshold="4.0"
pdockq_threshold="8.0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)             project="$2"; shift 2 ;;
        --pae-threshold)       pae_threshold="$2"; shift 2 ;;
        --pde-threshold)       pde_threshold="$2"; shift 2 ;;
        --interface-threshold) interface_threshold="$2"; shift 2 ;;
        --pdockq-threshold)    pdockq_threshold="$2"; shift 2 ;;
        --help|-h)             show_help; exit 0 ;;
        *)                     log_error "Unknown flag: $1"; show_help; exit 1 ;;
    esac
done

# Validate required flags
if [[ -z "$project" ]]; then
    log_error "Missing required flag: --project"
    show_help
    exit 1
fi

# Validate inputs
[[ -d "$project" ]] || { log_error "Project directory not found: $project"; exit 1; }
[[ -f "$project/.prediction_complete" ]] || { log_error "Structure prediction not complete (missing .prediction_complete)"; exit 1; }
[[ -f "$project/summary.csv" ]] || { log_error "summary.csv not found in project directory"; exit 1; }

DALPHABALL_PATH="/home/ubuntu/miniforge3/envs/pyrosetta/bin/dalphaball"

# Check if already complete
if [[ -f "$project/.scoring_complete" ]]; then
    log_warn "Scoring already complete (found .scoring_complete marker)"
    log_info "To regenerate, delete: $project/.scoring_complete"
    exit 0
fi

#===============================================================================
# Environment setup
#===============================================================================

log_stage "Stage 3: Metrics Extraction and Scoring"
log_info "Project: $project"
log_info "Thresholds: PAE=$pae_threshold, PDE=$pde_threshold, Interface=$interface_threshold, pDockQ=$pdockq_threshold"
log_info ""
source ~/miniforge3/etc/profile.d/conda.sh

# Check for required scripts
for script in colabfold_extract_metrics.py boltz_extract_metrics.py chai_extract_metrics.py binder_extract_metrics.py merge.py compute_rosetta_metrics.py compute_rmsd.py; do
    if [[ ! -f "$SCRIPT_DIR/$script" ]]; then
        log_error "Required script not found: $SCRIPT_DIR/$script"
        exit 1
    fi
done

#===============================================================================
# Directory setup
#===============================================================================

metrics_dir="$project/metrics"
output_dir="$project/filtered_structures"
symlinks_dir="$project/symlinks"
relaxed_dir="$project/relaxed"

rm -rf "$metrics_dir"
mkdir -p "$metrics_dir"

[[ -d "$output_dir" ]] || { log_error "Filtered structures directory not found: $output_dir"; exit 1; }

#===============================================================================
# Activate pyrosetta environment
#===============================================================================

log_info "Activating pyrosetta environment..."
conda activate pyrosetta

#===============================================================================
# Run metrics extraction scripts
#===============================================================================

log_stage "Extracting Metrics"

# --- AF2-multimer ---
if [[ -d "$output_dir/af2_multimer" ]] && [[ -n "$(ls -A "$output_dir/af2_multimer" 2>/dev/null)" ]]; then
    log_info ""
    log_info "Extracting AF2-multimer metrics..."
    python "$SCRIPT_DIR/colabfold_extract_metrics.py" \
        --parent_dir "$output_dir/af2_multimer" \
        --prefix "mul" \
        --pae_threshold "$pae_threshold" \
        --interface_threshold "$interface_threshold" \
        --pdockq_threshold "$pdockq_threshold" \
        --output_filename "$metrics_dir/af2_multimer_metrics.csv"
    log_info "AF2-multimer metrics saved to: $metrics_dir/af2_multimer_metrics.csv"
else
    log_warn "Skipping AF2-multimer (directory empty or missing)"
fi

# --- AF2-ptm complex ---
if [[ -d "$output_dir/af2_ptm_complex" ]] && [[ -n "$(ls -A "$output_dir/af2_ptm_complex" 2>/dev/null)" ]]; then
    log_info ""
    log_info "Extracting AF2-ptm complex metrics..."
    python "$SCRIPT_DIR/colabfold_extract_metrics.py" \
        --parent_dir "$output_dir/af2_ptm_complex" \
        --prefix "mon" \
        --pae_threshold "$pae_threshold" \
        --interface_threshold "$interface_threshold" \
        --pdockq_threshold "$pdockq_threshold" \
        --output_filename "$metrics_dir/af2_ptm_complex_metrics.csv"
    log_info "AF2-ptm complex metrics saved to: $metrics_dir/af2_ptm_complex_metrics.csv"
else
    log_warn "Skipping AF2-ptm complex (directory empty or missing)"
fi

# --- Boltz ---
if [[ -d "$output_dir/boltz" ]] && [[ -n "$(ls -A "$output_dir/boltz" 2>/dev/null)" ]]; then
    log_info ""
    log_info "Extracting Boltz metrics..."
    python "$SCRIPT_DIR/boltz_extract_metrics.py" \
        --parent_dir "$output_dir/boltz" \
        --prefix "bol" \
        --pae_threshold "$pae_threshold" \
        --pde_threshold "$pde_threshold" \
        --interface_threshold "$interface_threshold" \
        --pdockq_threshold "$pdockq_threshold" \
        --output_filename "$metrics_dir/boltz_metrics.csv"
    log_info "Boltz metrics saved to: $metrics_dir/boltz_metrics.csv"
else
    log_warn "Skipping Boltz (directory empty or missing)"
fi

# --- Chai ---
if [[ -d "$output_dir/chai" ]] && [[ -n "$(ls -A "$output_dir/chai" 2>/dev/null)" ]]; then
    log_info ""
    log_info "Extracting Chai metrics..."
    python "$SCRIPT_DIR/chai_extract_metrics.py" \
        --parent_dir "$output_dir/chai" \
        --prefix "cha" \
        --interface_threshold "$interface_threshold" \
        --pdockq_threshold "$pdockq_threshold" \
        --output_filename "$metrics_dir/chai_metrics.csv"
    log_info "Chai metrics saved to: $metrics_dir/chai_metrics.csv"
else
    log_warn "Skipping Chai (directory empty or missing)"
fi

# --- Binder ---
if [[ -d "$output_dir/binder" ]] && [[ -n "$(ls -A "$output_dir/binder" 2>/dev/null)" ]]; then
    log_info ""
    log_info "Extracting Binder metrics..."
    python "$SCRIPT_DIR/binder_extract_metrics.py" \
        --parent_dir "$output_dir/binder" \
        --prefix "bin" \
        --output_filename "$metrics_dir/binder_metrics.csv"
    log_info "Binder metrics saved to: $metrics_dir/binder_metrics.csv"
else
    log_warn "Skipping Binder (directory empty or missing)"
fi

#===============================================================================
# Merge into binderScore.csv
#===============================================================================

log_stage "Generating binderScore.csv"

log_info "Running merge script..."
python "$SCRIPT_DIR/merge.py" "$project" || {
    log_error "Merge script failed"
    conda deactivate
    exit 1
}

mv "$metrics_dir/binderScore.csv" "$project/binderScore.csv"

# Verify non-empty output
row_count=$(python -c "import pandas as pd; print(len(pd.read_csv('${project}/binderScore.csv')))")
if [[ "$row_count" -eq 0 ]]; then
    log_warn "binderScore.csv has 0 rows — no binders were scored"
fi

log_info "binderScore.csv generated successfully"

#===============================================================================
# Part 2: Create Symlinks
#===============================================================================

log_stage "Creating Structure Symlinks"

# Define methods to process
methods=(af2_multimer af2_ptm_complex boltz chai binder)

# Clean and create symlinks directory
rm -rf "$symlinks_dir"
mkdir -p "$symlinks_dir"

for method in "${methods[@]}"; do
    method_src="$output_dir/$method"
    method_dst="$symlinks_dir/$method"
    
    # Skip if source directory doesn't exist or is empty
    if [[ ! -d "$method_src" ]] || [[ -z "$(ls -A "$method_src" 2>/dev/null)" ]]; then
        log_warn "Skipping symlinks for $method (directory empty or missing)"
        continue
    fi
    
    mkdir -p "$method_dst"
    symlink_count=0
    
    # Iterate over binder subdirectories
    for binder_dir in "$method_src"/*/; do
        # Skip if not a directory
        [[ -d "$binder_dir" ]] || continue
        
        binder_name=$(basename "$binder_dir")
        
        # Find structure file: prefer .pdb over .cif
        structure_file=""
        
        # First, look for .pdb files
        pdb_file=$(find "$binder_dir" -maxdepth 1 -name "*.pdb" -type f | head -n 1)
        if [[ -n "$pdb_file" ]]; then
            structure_file="$pdb_file"
            ext="pdb"
        else
            # Fall back to .cif files
            cif_file=$(find "$binder_dir" -maxdepth 1 -name "*.cif" -type f | head -n 1)
            if [[ -n "$cif_file" ]]; then
                structure_file="$cif_file"
                ext="cif"
            fi
        fi
        
        # Create symlink if structure file found
        if [[ -n "$structure_file" ]]; then
            ln -sf "$(realpath "$structure_file")" "$method_dst/${binder_name}.${ext}"
            symlink_count=$((symlink_count+1))
            
        else
            log_warn "No PDB/CIF found for $method/$binder_name"
        fi
    done
    
    log_info "Created $symlink_count symlinks for $method"
done

log_info "Symlinks created in: $symlinks_dir"

#===============================================================================
# Part 3: Compute Rosetta Metrics
#===============================================================================

log_stage "Computing Rosetta Interface Metrics"

# Create relaxed structures directory
mkdir -p "$relaxed_dir"

# Define method to prefix mapping (matching existing metrics prefixes)
declare -A method_prefixes=(
    ["af2_multimer"]="mul"
    ["af2_ptm_complex"]="mon"
    ["boltz"]="bol"
    ["chai"]="cha"
    ["binder"]="bin"
)

for method in "${methods[@]}"; do
    # Skip binder method - single-chain structures have no interface to score
    if [[ "$method" == "binder" ]]; then
        log_info "Skipping Rosetta interface metrics for $method (single-chain, no interface)"
        continue
    fi

    method_symlinks="$symlinks_dir/$method"
    prefix="${method_prefixes[$method]}"
    
    # Skip if symlinks directory doesn't exist or is empty
    if [[ ! -d "$method_symlinks" ]] || [[ -z "$(ls -A "$method_symlinks" 2>/dev/null)" ]]; then
        log_warn "Skipping Rosetta metrics for $method (no symlinks found)"
        continue
    fi
    
    # Count structure files (both PDB and CIF are now supported)
    structure_count=$(find "$method_symlinks" \( -name "*.pdb" -o -name "*.cif" \) -type l | wc -l)
    if [[ "$structure_count" -eq 0 ]]; then
        log_warn "Skipping Rosetta metrics for $method (no structure files found)"
        continue
    fi
    
    log_info ""
    log_info "Computing Rosetta metrics for $method ($structure_count structures, prefix: $prefix)..."
    
    python "$SCRIPT_DIR/compute_rosetta_metrics.py" \
        --folder "$prefix:$method_symlinks" \
        --binder-chain B \
        --target-chain A \
        --out-csv "$metrics_dir/${method}_rosettaScore.csv" \
        --relaxed-dir "$relaxed_dir/$method" \
        --dalphaball-path "$DALPHABALL_PATH" \
    || {
        log_error "Rosetta metrics failed for $method"
        # Continue with other methods instead of exiting
        continue
    }
    
    log_info "Rosetta metrics saved to: $metrics_dir/${method}_rosettaScore.csv"
done

log_info "Relaxed structures cached in: $relaxed_dir"

#===============================================================================
# Part 4: Merge Rosetta Metrics into binderRosettaScore.csv
#===============================================================================

log_stage "Merging Rosetta Metrics"

# Collect all rosettaScore.csv files that exist
rosetta_csvs=()
for method in "${methods[@]}"; do
    csv_file="$metrics_dir/${method}_rosettaScore.csv"
    if [[ -f "$csv_file" ]]; then
        rosetta_csvs+=("$csv_file")
    fi
done

if [[ ${#rosetta_csvs[@]} -eq 0 ]]; then
    log_warn "No Rosetta metrics CSVs found to merge"
else
    log_info "Merging ${#rosetta_csvs[@]} Rosetta metrics files..."
    
    # Use Python to merge CSVs on binder_id (outer join)
    python - "${rosetta_csvs[@]}" "$metrics_dir/binderRosettaScore.csv" <<'PYTHON_SCRIPT'
import sys
import pandas as pd

csv_files = sys.argv[1:-1]
output_file = sys.argv[-1]

# Read and merge all CSVs
merged_df = None
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on="binder", how="outer")

# Write merged CSV
if merged_df is not None:
    merged_df.to_csv(output_file, index=False)
    print(f"[INFO] Merged {len(csv_files)} files into {output_file} ({len(merged_df)} rows, {len(merged_df.columns)} columns)")
PYTHON_SCRIPT

    log_info "Merged Rosetta metrics saved to: $metrics_dir/binderRosettaScore.csv"
fi

if [[ -f "$metrics_dir/binderRosettaScore.csv" ]]; then
    mv "$metrics_dir/binderRosettaScore.csv" "$project/binderRosettaScore.csv"
fi

if [[ -f "$project/binderScore.csv" ]] && [[ -f "$project/binderRosettaScore.csv" ]]; then
    python -c "
import pandas as pd
binder = pd.read_csv('${project}/binderScore.csv')
rosetta = pd.read_csv('${project}/binderRosettaScore.csv')
merged = pd.merge(binder, rosetta, on='binder', how='outer')
merged.to_csv('${project}/binderFinalScore.csv', index=False)
print(f'[INFO] Merged final score: {len(merged)} rows, {len(merged.columns)} columns')
"
else
    log_warn "Skipping final merge (binderScore.csv or binderRosettaScore.csv missing)"
fi

#===============================================================================
# Part 5: Compute Cross-Method RMSD
#===============================================================================

log_stage "Computing Cross-Method Binder RMSD"

if [[ -d "$relaxed_dir" ]] && [[ -n "$(ls -A "$relaxed_dir" 2>/dev/null)" ]]; then
    python "$SCRIPT_DIR/compute_rmsd.py" \
        --project "$project" \
    || {
        log_error "RMSD computation failed"
    }

    # Merge RMSD.csv into binderFinalScore.csv
    if [[ -f "$project/RMSD.csv" ]] && [[ -f "$project/binderFinalScore.csv" ]]; then
        python -c "
import pandas as pd
final = pd.read_csv('${project}/binderFinalScore.csv')
rmsd = pd.read_csv('${project}/RMSD.csv')
merged = pd.merge(final, rmsd, on='binder', how='outer')
merged.to_csv('${project}/binderFinalScore.csv', index=False)
print(f'[INFO] Merged RMSD into binderFinalScore: {len(merged)} rows, {len(merged.columns)} columns')
"
    elif [[ -f "$project/RMSD.csv" ]] && [[ ! -f "$project/binderFinalScore.csv" ]]; then
        log_warn "binderFinalScore.csv not found, RMSD.csv saved standalone"
    fi
else
    log_warn "Skipping RMSD computation (relaxed directory empty or missing)"
fi

#===============================================================================
# Deactivate environment
#===============================================================================

conda deactivate

#===============================================================================
# Completion
#===============================================================================

log_stage "Scoring Complete"

# Create completion marker
cat > "$project/.scoring_complete" <<EOF
# Scoring Complete
# Generated: $(date)
# Thresholds: PAE=$pae_threshold, PDE=$pde_threshold, Interface=$interface_threshold, pDockQ=$pdockq_threshold
# Output: $metrics_dir/binderScore.csv
# Rosetta metrics: $metrics_dir/*_rosettaScore.csv
# Symlinks: $symlinks_dir
# Relaxed structures: $relaxed_dir
EOF

log_info "Results:"
log_info "  - Individual metrics: $metrics_dir/*_metrics.csv"
log_info "  - Final scores: ./binderScore.csv"
log_info "  - Rosetta metrics (per method): $metrics_dir/*_rosettaScore.csv"
log_info "  - Rosetta metrics (merged): $metrics_dir/binderRosettaScore.csv"
log_info "  - Structure symlinks: $symlinks_dir"
log_info "  - Relaxed structures: $relaxed_dir"
log_info "Completion marker: $project/.scoring_complete"
log_info ""