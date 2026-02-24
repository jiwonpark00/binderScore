#!/bin/bash
#===============================================================================
# rerun_filter.sh - Rerun filter_and_select.py with per-target distance cutoffs
#
# This script allows you to rerun the filtering step without redoing structure
# prediction. Useful when you want to adjust per-target distance cutoffs.
#
# The residue CSV (with target_name, residue_number, distance_cutoff columns)
# is read from the .prediction_complete marker by default, or can be overridden
# with --residue-csv.
#
# Usage:
#   ./rerun_filter.sh <project_dir> [--residue-csv /path/to/residues.csv]
#
# Example:
#   # Re-filter using the original residue CSV from the pipeline run
#   ./rerun_filter.sh ./my_project
#
#   # Re-filter with updated per-target cutoffs in a new residue CSV
#   ./rerun_filter.sh ./my_project --residue-csv /path/to/updated_residues.csv
#
# Requirements:
#   - distances.csv must exist in project_dir (from previous run)
#   - predictions/ directory must exist
#   - Python environment with required packages
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#===============================================================================
# Logging functions
#===============================================================================

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') $*"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2
}

#===============================================================================
# Argument parsing and validation
#===============================================================================

if [[ $# -lt 1 ]]; then
    log_error "Usage: $0 <project_dir> [--residue-csv /path/to/residues.csv]"
    log_error ""
    log_error "Arguments:"
    log_error "  project_dir    - Project directory (must contain distances.csv)"
    log_error ""
    log_error "Options:"
    log_error "  --residue-csv  - Path to residue CSV with per-target cutoffs"
    log_error "                   (default: reads path from .prediction_complete marker)"
    log_error ""
    log_error "The residue CSV must have columns: target_name, residue_number, distance_cutoff"
    log_error ""
    log_error "Example:"
    log_error "  $0 ./my_project"
    log_error "  $0 ./my_project --residue-csv updated_residues.csv"
    exit 1
fi

project="$1"
shift

# Parse optional --residue-csv flag
residue_csv=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --residue-csv) residue_csv="$2"; shift 2 ;;
        *) log_error "Unknown flag: $1"; exit 1 ;;
    esac
done

# Validate project directory
if [[ ! -d "$project" ]]; then
    log_error "Project directory not found: $project"
    exit 1
fi

# Validate distances.csv exists
if [[ ! -f "$project/distances.csv" ]]; then
    log_error "distances.csv not found in $project"
    log_error "Run the full pipeline first to generate distances.csv"
    exit 1
fi

# Validate predictions directory exists
if [[ ! -d "$project/predictions" ]]; then
    log_error "predictions/ directory not found in $project"
    exit 1
fi

# If no --residue-csv provided, read from .prediction_complete marker
if [[ -z "$residue_csv" ]]; then
    if [[ -f "$project/.prediction_complete" ]]; then
        residue_csv=$(grep "^# Residue CSV:" "$project/.prediction_complete" | sed 's/^# Residue CSV: //')
    fi
    if [[ -z "$residue_csv" ]]; then
        log_error "Could not determine residue CSV path from .prediction_complete marker"
        log_error "Please provide --residue-csv explicitly"
        exit 1
    fi
fi

# Validate residue CSV exists
if [[ ! -f "$residue_csv" ]]; then
    log_error "Residue CSV not found: $residue_csv"
    exit 1
fi

#===============================================================================
# Main
#===============================================================================

log_info "Rerunning filter with per-target distance cutoffs"
log_info "  Project: $project"
log_info "  Residue CSV: $residue_csv"
log_info "  Distances CSV: $project/distances.csv"

# Read binder_msa_mode from .prediction_complete marker
binder_msa_mode=""
if [[ -f "$project/.prediction_complete" ]]; then
    binder_msa_mode=$(grep "^# Binder MSA mode:" "$project/.prediction_complete" | sed 's/^# Binder MSA mode: //')
fi
if [[ -z "$binder_msa_mode" ]]; then
    log_warn "Could not determine binder_msa_mode from .prediction_complete, defaulting to 'false'"
    binder_msa_mode="false"
fi
log_info "  Binder MSA mode: $binder_msa_mode"

# Clear old filtered_structures
if [[ -d "$project/filtered_structures" ]]; then
    log_info "Removing old filtered_structures directory..."
    rm -rf "$project/filtered_structures"
fi

# Backup old summary if it exists
if [[ -f "$project/summary.csv" ]]; then
    backup_name="$project/summary_$(date '+%Y%m%d_%H%M%S').csv.bak"
    log_info "Backing up old summary.csv to $backup_name"
    mv "$project/summary.csv" "$backup_name"
fi

# Activate conda environment if available
if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    # Try to activate chai environment (has biopython)
    if conda activate chai 2>/dev/null; then
        log_info "Activated conda environment: chai"
    fi
fi

# Run filter_and_select.py
log_info "Running filter_and_select.py..."

# Check if retry distances exist
retry_args=""
if [[ -f "$project/distances_retry.csv" && -d "$project/predictions_retry" ]]; then
    log_info "Found retry data, including in filter..."
    retry_args="--distances_csv_retry $project/distances_retry.csv --predictions_dir_retry $project/predictions_retry"
fi

python "$SCRIPT_DIR/filter_and_select.py" \
    --distances_csv "$project/distances.csv" \
    --predictions_dir "$project/predictions" \
    --residue_csv "$residue_csv" \
    --output_dir "$project/filtered_structures" \
    --summary_csv "$project/summary.csv" \
    --binder_msa "$binder_msa_mode" \
    $retry_args

log_info "Done!"
log_info "Results:"
log_info "  - Filtered structures: $project/filtered_structures/"
log_info "  - Summary: $project/summary.csv"