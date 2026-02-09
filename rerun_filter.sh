#!/bin/bash
#===============================================================================
# rerun_filter.sh - Rerun filter_and_select.py with a different distance cutoff
#
# This script allows you to rerun the filtering step without redoing structure
# prediction. Useful when too many structures failed the initial cutoff.
#
# Usage:
#   ./rerun_filter.sh <project_dir> <new_cutoff>
#
# Example:
#   ./rerun_filter.sh ./my_project 15.0
#   ./rerun_filter.sh /path/to/project 25.0
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

if [[ $# -lt 2 ]]; then
    log_error "Usage: $0 <project_dir> <new_cutoff>"
    log_error ""
    log_error "Arguments:"
    log_error "  project_dir  - Project directory (must contain distances.csv)"
    log_error "  new_cutoff   - New distance cutoff in Angstroms (e.g., 15.0, 25.0)"
    log_error ""
    log_error "Example:"
    log_error "  $0 ./my_project 15.0"
    exit 1
fi

project="$1"
cutoff="$2"

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

# Validate cutoff is a number
if ! [[ "$cutoff" =~ ^[0-9]+\.?[0-9]*$ ]]; then
    log_error "Invalid cutoff value: $cutoff (must be a positive number)"
    exit 1
fi

#===============================================================================
# Main
#===============================================================================

log_info "Rerunning filter with new cutoff"
log_info "  Project: $project"
log_info "  New cutoff: $cutoff Å"
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
    --cutoff "$cutoff" \
    --output_dir "$project/filtered_structures" \
    --summary_csv "$project/summary.csv" \
    --binder_msa "$binder_msa_mode" \
    $retry_args

log_info "Done!"
log_info "Results:"
log_info "  - Filtered structures: $project/filtered_structures/"
log_info "  - Summary: $project/summary.csv"