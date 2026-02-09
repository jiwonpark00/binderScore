#!/usr/bin/env bash
#===============================================================================
# Stage 2: Structure Prediction for binder screening pipeline
#
# Inputs (named flags):
#   --csv              - CSV file (name,target,binder)
#   --project          - Project directory (contains MSA outputs from generate_msas.sh)
#   --target-dir       - Target directory (contains /msa and /pdb subdirectories)
#   --residue-csv      - Target residue CSV file (target_name,residue_number)
#   --distance-cutoff  - Distance cutoff in Angstroms (e.g., 20)
#   --binder-msa       - binder_msa_mode flag (true/false)
#
# Optional flags - ColabFold (AF2):
#   --num-recycle             - Number of recycles for AF2 (default: 2)
#   --recycle-early-stop-tol  - Early stop tolerance for AF2 (default: 0.5)
#   --num-seeds               - Number of seeds to try per model (default: 1)
#   --num-models              - Number of AF2 models to run (default: 5)
#
# Optional flags - Chai:
#   --chai-num-trunk-recycles  - Number of trunk recycles (default: 3)
#   --chai-num-diffn-timesteps - Number of diffusion timesteps (default: 200)
#   --chai-diffusion-samples   - Number of diffusion samples (default: 5)
#
# Optional flags - Boltz:
#   --boltz-recycling-steps    - Number of recycling steps (default: 3)
#   --boltz-sampling-steps     - Number of sampling steps (default: 200)
#   --boltz-diffusion-samples  - Number of diffusion samples (default: 5)
#   --use-potentials           - Use inference-time potentials (default: true)
#
# Outputs:
#   - predictions/{af2_multimer,af2_ptm_complex,chai,boltz,binder}/
#   - distances.csv
#   - filtered_structures/
#   - summary.csv
#===============================================================================

set -euo pipefail
shopt -s nullglob
export DISABLE_PANDERA_IMPORT_WARNING=True
export PYTHONWARNINGS="ignore::UserWarning"

# NOTE: set -e is active during setup/validation. It is disabled before
# running prediction methods so that individual method failures don't
# kill the entire pipeline (fail-loud-but-continue pattern).

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
# Helper function: Check if predictions exist in output directory
#===============================================================================

predictions_exist() {
    local output_dir="$1"
    [[ -d "$output_dir" ]] || return 1
    # Check if find returns any results
    [ -n "$(find "$output_dir" -type f \( -name "*.pdb" -o -name "*.cif" \) -print -quit)" ]
}

#===============================================================================
# Help message
#===============================================================================

show_help() {
    cat <<EOF
Usage: $0 [OPTIONS]

Stage 2: Structure Prediction for binder screening pipeline

Required flags:
  --csv                        CSV file with name,target,binder columns
  --project                    Project directory (from generate_msas.sh)
  --target-dir                 Target directory (contains /msa and /pdb)
  --residue-csv                CSV with target_name,residue_number columns
  --distance-cutoff            Distance cutoff in Angstroms
  --binder-msa                 Whether binder MSAs exist (true/false)

ColabFold (AF2) optional flags:
  --num-recycle                Number of recycles for AF2 (default: 2)
  --recycle-early-stop-tol     Early stop tolerance for AF2 recycles (default: 0.5)
  --num-seeds                  Number of seeds per model (default: 1)
  --num-models                 Number of AF2 models to run, 1-5 (default: 5)

Chai optional flags:
  --chai-num-trunk-recycles    Number of trunk recycles (default: 3)
  --chai-num-diffn-timesteps   Number of diffusion timesteps (default: 200)
  --chai-diffusion-samples     Number of diffusion samples (default: 5)

Boltz optional flags:
  --boltz-recycling-steps      Number of recycling steps (default: 3)
  --boltz-sampling-steps       Number of sampling steps (default: 200)
  --boltz-diffusion-samples    Number of diffusion samples (default: 5)
  --use-potentials             Use inference-time potentials (default: true)

Options:
  --resume                     Resume from existing outputs (skip completed steps)
  --help, -h                   Show this help message and exit

Example:
  $0 --csv input.csv --project project/ --target-dir targets/ \\
     --residue-csv residues.csv --distance-cutoff 20 --binder-msa false

  $0 --csv input.csv --project project/ --target-dir targets/ \\
     --residue-csv residues.csv --distance-cutoff 20 --binder-msa true \\
     --num-recycle 3 --boltz-diffusion-samples 10 --use-potentials false
EOF
}

#===============================================================================
# Argument parsing
#===============================================================================

csv=""
project=""
target_dir=""
residue_csv=""
distance_cutoff=""
binder_msa_mode=""
# ColabFold defaults
num_recycle="2"
recycle_early_stop_tol="0.5"
num_seeds="1"
num_models="5"
# Chai defaults
chai_num_trunk_recycles="3"
chai_num_diffn_timesteps="200"
chai_diffusion_samples="5"
# Boltz defaults
boltz_recycling_steps="3"
boltz_sampling_steps="200"
boltz_diffusion_samples="5"
use_potentials="true"
# Resume mode
resume="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --csv)                         csv="$2"; shift 2 ;;
        --project)                     project="$2"; shift 2 ;;
        --target-dir)                  target_dir="$2"; shift 2 ;;
        --residue-csv)                 residue_csv="$2"; shift 2 ;;
        --distance-cutoff)             distance_cutoff="$2"; shift 2 ;;
        --binder-msa)                  binder_msa_mode="$2"; shift 2 ;;
        # ColabFold flags
        --num-recycle)                 num_recycle="$2"; shift 2 ;;
        --recycle-early-stop-tol)      recycle_early_stop_tol="$2"; shift 2 ;;
        --num-seeds)                   num_seeds="$2"; shift 2 ;;
        --num-models)                  num_models="$2"; shift 2 ;;
        # Chai flags
        --chai-num-trunk-recycles)     chai_num_trunk_recycles="$2"; shift 2 ;;
        --chai-num-diffn-timesteps)    chai_num_diffn_timesteps="$2"; shift 2 ;;
        --chai-diffusion-samples)      chai_diffusion_samples="$2"; shift 2 ;;
        # Boltz flags
        --boltz-recycling-steps)       boltz_recycling_steps="$2"; shift 2 ;;
        --boltz-sampling-steps)        boltz_sampling_steps="$2"; shift 2 ;;
        --boltz-diffusion-samples)     boltz_diffusion_samples="$2"; shift 2 ;;
        --use-potentials)              use_potentials="$2"; shift 2 ;;
        --resume)                      resume="true"; shift ;;
        --help|-h)                     show_help; exit 0 ;;
        *)                             log_error "Unknown flag: $1"; show_help; exit 1 ;;
    esac
done

# Validate required flags
missing=()
[[ -z "$csv" ]]             && missing+=("--csv")
[[ -z "$project" ]]         && missing+=("--project")
[[ -z "$target_dir" ]]      && missing+=("--target-dir")
[[ -z "$residue_csv" ]]     && missing+=("--residue-csv")
[[ -z "$distance_cutoff" ]] && missing+=("--distance-cutoff")
[[ -z "$binder_msa_mode" ]] && missing+=("--binder-msa")

if [[ ${#missing[@]} -gt 0 ]]; then
    log_error "Missing required flags: ${missing[*]}"
    show_help
    exit 1
fi

# Validate inputs
[[ -f "$csv" ]] || { log_error "CSV file not found: $csv"; exit 1; }
[[ -d "$project" ]] || { log_error "Project directory not found: $project"; exit 1; }
[[ -f "$project/.msa_complete" ]] || { log_error "MSA generation not complete (missing .msa_complete)"; exit 1; }
[[ -d "$target_dir" ]] || { log_error "Target directory not found: $target_dir"; exit 1; }
[[ -d "$target_dir/msa" ]] || { log_error "Target MSA directory not found: $target_dir/msa"; exit 1; }
[[ -d "$target_dir/pdb" ]] || { log_error "Target PDB directory not found: $target_dir/pdb"; exit 1; }
[[ -f "$residue_csv" ]] || { log_error "Residue CSV file not found: $residue_csv"; exit 1; }

if [[ "$binder_msa_mode" != "true" && "$binder_msa_mode" != "false" ]]; then
    log_error "binder_msa_mode must be 'true' or 'false', got: $binder_msa_mode"
    exit 1
fi

if [[ "$use_potentials" != "true" && "$use_potentials" != "false" ]]; then
    log_error "use_potentials must be 'true' or 'false', got: $use_potentials"
    exit 1
fi

# Ensure newline at end of CSV
lastchar=$(tail -c 1 "$csv")
if [[ "$lastchar" != '' ]]; then
    printf '\n' >> "$csv"
fi

#===============================================================================
# Environment setup
#===============================================================================

log_stage "Stage 2: Structure Prediction"
log_info "CSV: $csv"
log_info "Project: $project"
log_info "Target dir: $target_dir"
log_info "Residue CSV: $residue_csv"
log_info "Distance cutoff: $distance_cutoff Å"
log_info "Binder MSA mode: $binder_msa_mode"
log_info "--- ColabFold (AF2) ---"
log_info "  Num recycle: $num_recycle"
log_info "  Recycle early stop tolerance: $recycle_early_stop_tol"
log_info "  Num seeds: $num_seeds"
log_info "  Num models: $num_models"
log_info "--- Chai ---"
log_info "  Num trunk recycles: $chai_num_trunk_recycles"
log_info "  Num diffusion timesteps: $chai_num_diffn_timesteps"
log_info "  Diffusion samples: $chai_diffusion_samples"
log_info "--- Boltz ---"
log_info "  Recycling steps: $boltz_recycling_steps"
log_info "  Sampling steps: $boltz_sampling_steps"
log_info "  Diffusion samples: $boltz_diffusion_samples"
log_info "  Use potentials: $use_potentials"

source $HOME/miniforge3/etc/profile.d/conda.sh

# Track errors across prediction methods (fail-loud-but-continue)
pipeline_errors=()

# Prepare Boltz --use_potentials flag
if [[ "$use_potentials" == "true" ]]; then
    use_potentials_flag="--use_potentials"
else
    use_potentials_flag=""
fi

# Generate random seeds for ColabFold (0-10 range for reproducibility across small set)
cf_random_seed=$(shuf -i 0-10 -n 1)
cf_retry_seed=$(shuf -i 0-10 -n 1)
log_info "ColabFold random seed: $cf_random_seed"
log_info "ColabFold retry random seed: $cf_retry_seed"

# Detect GPUs
NGPUS=$(nvidia-smi -L 2>/dev/null | wc -l)
if (( NGPUS == 0 )); then
    log_error "No GPUs detected"
    exit 1
fi

TOTAL_CPUS=$(nproc)
WORKERS_PER_GPU=$(( (TOTAL_CPUS * 7 / 10) / NGPUS ))
[[ $WORKERS_PER_GPU -lt 1 ]] && WORKERS_PER_GPU=1

CUDA_VISIBLE_DEVICES=$(printf '%s,' $(seq 0 $((NGPUS-1))))
CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES%,}
export CUDA_VISIBLE_DEVICES

log_info "GPUs: $NGPUS, Total CPUs: $TOTAL_CPUS, Workers/GPU: $WORKERS_PER_GPU"

# Check for required scripts
for script in generate_boltz_yaml.py calculate_distances.py filter_and_select.py prepare_retry_inputs.py chai_run.py; do
    if [[ ! -f "$SCRIPT_DIR/$script" ]]; then
        log_error "Required script not found: $SCRIPT_DIR/$script"
        exit 1
    fi
done

#===============================================================================
# Directory setup
#===============================================================================

pred_dir="$project/predictions"
mkdir -p "$pred_dir"/{af2_multimer,af2_ptm_complex,chai,boltz,binder}
mkdir -p "$project/filtered_structures"

# MSA directories
multimer_msa_dir="$project/colabfold/msa_multimer"
binder_msa_dir="$project/colabfold/msa_binder"
chai_msa_dir="$project/chai/msas"

#===============================================================================
# Helper function: Parallel ColabFold batch
#===============================================================================

parallel_colabfold_batch() {
    local job_name="$1"; shift
    local input_dir="$1"; shift
    local output_dir="$1"; shift
    local -a extra_args=( "$@" )

    local -a a3m_files=( "$input_dir"/*.a3m )
    if (( ${#a3m_files[@]} == 0 )); then
        log_warn "[$job_name] No .a3m files in $input_dir"
        return 1
    fi

    log_info "[$job_name] Found ${#a3m_files[@]} a3m files"
    mkdir -p "$output_dir"

    # Create per-GPU input directories
    for g in $(seq 0 $((NGPUS-1))); do
        mkdir -p "$input_dir/in_$g" "$output_dir/gpu_$g"
    done

    # Distribute a3m files across GPUs
    local i=0
    for f in "${a3m_files[@]}"; do
        local g=$(( i % NGPUS ))
        cp "$f" "$input_dir/in_$g/"
        i=$((i+1))
    done

    # Launch parallel jobs
    local -a pids=()
    for g in $(seq 0 $((NGPUS-1))); do
        local in_dir="$input_dir/in_$g"
        local out_dir="$output_dir/gpu_$g"
        local -a files=( "$in_dir"/*.a3m )
        (( ${#files[@]} )) || continue

        log_info "[$job_name] GPU $g: ${#files[@]} a3m files"

        (
            export CUDA_VISIBLE_DEVICES="$g"
            colabfold_batch \
                "${extra_args[@]}" \
                "$in_dir" "$out_dir/"
        ) >"$out_dir/run.log" 2>&1 &

        pids+=( "$!" )
    done

    # Wait for completion
    local fail=0
    for pid in "${pids[@]}"; do
        wait "$pid" || fail=1
    done

    if (( fail )); then
        log_error "[$job_name] One or more workers failed. Check $output_dir/gpu_*/run.log"
        return 1
    fi

    log_info "[$job_name] Complete. Outputs in $output_dir"
    return 0
}

#===============================================================================
# Helper function: Run Chai-1 predictions
#===============================================================================

run_chai_predictions() {
    local fasta_dir="$1"
    local msa_dir="$2"
    local output_dir="$3"
    local num_samples="$4"
    local num_trunk_recycles="$5"
    local num_diffn_timesteps="$6"

    local -a fasta_files=( "$fasta_dir"/*.fasta )
    if (( ${#fasta_files[@]} == 0 )); then
        log_warn "[chai] No FASTA files in $fasta_dir"
        return 1
    fi

    log_info "[chai] Found ${#fasta_files[@]} FASTA files"
    mkdir -p "$output_dir/chai_log"

    local -a pids=()
    local i=0
    local fail=0

    for f in "${fasta_files[@]}"; do
        local gpu=$(( i % NGPUS ))
        local name="$(basename "$f" .fasta)"
        local out="$output_dir/$name"
        mkdir -p "$out"

        log_info "[chai] Running $name on GPU $gpu"

        (
            local msa_args=""
            if [[ -n "$msa_dir" && -d "$msa_dir" ]]; then
                msa_args="--msa_directory $msa_dir"
            fi

            python $SCRIPT_DIR/chai_run.py \
                --fasta "$f" \
                --output_dir "$out" \
                --gpu "$gpu" \
                --num_diffn_samples "$num_samples" \
                --num_trunk_recycles "$num_trunk_recycles" \
                --num_diffn_timesteps "$num_diffn_timesteps" \
                $msa_args
        ) >"$output_dir/chai_log/${name}_run.log" 2>&1 &

        pids+=( "$!" )
        i=$((i+1))

        # Throttle to NGPUS concurrent jobs
        if (( ${#pids[@]} >= NGPUS )); then
            wait "${pids[0]}" || fail=1
            pids=( "${pids[@]:1}" )
        fi
    done

    # Wait for remaining
    for pid in "${pids[@]}"; do
        wait "$pid" || fail=1
    done

    if (( fail )); then
        log_error "[chai] One or more workers failed. Check $output_dir/chai_log/*_run.log"
        return 1
    fi

    log_info "[chai] Complete. Outputs in $output_dir"
    return 0
}

#===============================================================================
# Helper function: Flatten gpu_* directories and cleanup in_* directories
#===============================================================================

flatten_gpu_dirs() {
    local output_dir="$1"
    local input_dir="$2"

    # Consolidate gpu_* outputs into parent
    for sd in "$output_dir"/gpu_*; do
        [[ -d "$sd" ]] || continue
        mv "$sd"/* "$output_dir/" 2>/dev/null || true
        rmdir "$sd" 2>/dev/null || true
    done

    # Cleanup temp per-GPU input directories
    rm -rf "$input_dir"/in_*/
}

#===============================================================================
# Helper functions: Determine completed binders per method
#===============================================================================

# AF2 (multimer & ptm_complex & binder-AF2): done when <name>.done.txt exists
get_done_binders_af2() {
    local output_dir="$1"
    [[ -d "$output_dir" ]] || return 0
    for f in "$output_dir"/*.done.txt; do
        [[ -f "$f" ]] || continue
        basename "$f" .done.txt
    done
}

# Chai: done when subdirectory <name>/ exists with at least one .cif or .pdb
get_done_binders_chai() {
    local output_dir="$1"
    [[ -d "$output_dir" ]] || return 0
    for d in "$output_dir"/*/; do
        [[ -d "$d" ]] || continue
        local name
        name="$(basename "$d")"
        # Skip log directories
        [[ "$name" == "chai_log" ]] && continue
        if [ -n "$(find "$d" -maxdepth 2 -type f \( -name '*.pdb' -o -name '*.cif' \) -print -quit)" ]; then
            echo "$name"
        fi
    done
}

# Boltz: done when predictions/<name>/ exists under boltz_results_yaml/
get_done_binders_boltz() {
    local output_dir="$1"
    local boltz_pred_dir="$output_dir/boltz_results_yaml/predictions"
    [[ -d "$boltz_pred_dir" ]] || return 0
    for d in "$boltz_pred_dir"/*/; do
        [[ -d "$d" ]] || continue
        local name
        name="$(basename "$d")"
        if [ -n "$(find "$d" -maxdepth 2 -type f \( -name '*.pdb' -o -name '*.cif' \) -print -quit)" ]; then
            echo "$name"
        fi
    done
}

#===============================================================================
# Helper function: Prepare a subset input directory for resume
#===============================================================================

# Creates a temporary subset directory containing only input files for
# binders that have NOT completed. Copies files (not symlinks) for safety.
# Args: source_dir, file_extension (.a3m/.fasta/.yaml), done_file (one name per line)
# Prints: path to the temporary subset directory
prepare_resume_subset() {
    local source_dir="$1"
    local ext="$2"
    local done_file="$3"

    local subset_dir
    subset_dir="$(mktemp -d "${source_dir}/resume_subset_XXXX")"

    for f in "$source_dir"/*"$ext"; do
        [[ -f "$f" ]] || continue
        local name
        name="$(basename "$f" "$ext")"
        # If this binder is NOT in the done list, copy it to subset
        if ! grep -qxF "$name" "$done_file"; then
            cp "$f" "$subset_dir/"
        fi
    done

    echo "$subset_dir"
}

#===============================================================================
# Prepare inputs
#===============================================================================

# Check if pipeline already complete
if [[ "$resume" == "true" ]] && [[ -f "$project/.prediction_complete" ]]; then
    log_info "Pipeline already complete (found .prediction_complete marker)"
    log_info "Results available at:"
    log_info "  - Predictions: $project/predictions"
    log_info "  - Filtered structures: $project/filtered_structures/"
    log_info "  - Summary: $project/summary.csv"
    exit 0
fi

log_stage "Preparing Inputs"

# --- Chai-1 FASTA files (multimer) ---
chai_fasta_dir="$project/chai/fasta_multimer"
mkdir -p "$chai_fasta_dir"

log_info "Generating Chai-1 FASTA files..."
awk -F',' '
NR==1 { next }
{
    gsub(/\r/, "", $0)
    name=$1; target=$2; binder=$3
    gsub(/^[ \t]+|[ \t]+$/, "", name)
    gsub(/^[ \t]+|[ \t]+$/, "", target)
    gsub(/^[ \t]+|[ \t]+$/, "", binder)
    if (name == "" || target == "" || binder == "") next

    out = outdir "/" name ".fasta"
    print ">protein|name=target" > out
    print target >> out
    print ">protein|name=binder" >> out
    print binder >> out
}
' outdir="$chai_fasta_dir" "$csv"

n_chai_fasta=$(find "$chai_fasta_dir" -name "*.fasta" -type f | wc -l)
log_info "Created $n_chai_fasta Chai-1 FASTA files"

# Ensure we're in a conda environment with Python for the next steps
conda activate chai

# --- Boltz-2 YAML files ---
# WARNING: Do not rename this directory. Boltz derives its output subdirectory
# name from the input directory name (e.g., "yaml" -> "boltz_results_yaml").
# Changing this breaks find_structures_boltz() in calculate_distances.py.
boltz_yaml_dir="$project/boltz/yaml"
mkdir -p "$boltz_yaml_dir"

log_info "Generating Boltz-2 YAML files..."
python $SCRIPT_DIR/generate_boltz_yaml.py \
    --csv "$csv" \
    --target_msa_dir "$target_dir/msa" \
    --target_pdb_dir "$target_dir/pdb" \
    --binder_msa_dir "$binder_msa_dir" \
    --binder_msa_mode "$binder_msa_mode" \
    --output_dir "$boltz_yaml_dir"

n_boltz_yaml=$(find "$boltz_yaml_dir" -name "*.yaml" -type f | wc -l)
log_info "Created $n_boltz_yaml Boltz-2 YAML files"

conda deactivate

# --- Binder-only FASTA files (for Chai no-MSA mode) ---
if [[ "$binder_msa_mode" == "false" ]]; then
    binder_fasta_dir="$project/chai/fasta_binder"
    mkdir -p "$binder_fasta_dir"

    log_info "Generating binder-only FASTA files for Chai-1..."
    awk -F',' '
    NR==1 { next }
    {
        gsub(/\r/, "", $0)
        name=$1; binder=$3
        gsub(/^[ \t]+|[ \t]+$/, "", name)
        gsub(/^[ \t]+|[ \t]+$/, "", binder)
        if (name == "" || binder == "") next

        out = outdir "/" name ".fasta"
        print ">protein|name=binder" > out
        print binder >> out
    }
    ' outdir="$binder_fasta_dir" "$csv"
fi

#===============================================================================
# Run Structure Predictions
#===============================================================================

log_stage "Running Structure Predictions"

# Disable set -e: individual method failures are tracked in pipeline_errors
set +e

# --- 1. ColabFold AF2-multimer ---
run_af2_multimer=true
af2_multimer_input_dir="$multimer_msa_dir"

if [[ "$resume" == "true" ]]; then
    # Always flatten first in case a previous run left gpu_* dirs
    flatten_gpu_dirs "$pred_dir/af2_multimer" "$multimer_msa_dir"

    done_file=$(mktemp)
    get_done_binders_af2 "$pred_dir/af2_multimer" > "$done_file"
    n_done=$(wc -l < "$done_file")
    n_total=$(find "$multimer_msa_dir" -maxdepth 1 -name '*.a3m' -type f 2>/dev/null | wc -l)
    n_remaining=$((n_total - n_done))

    if (( n_remaining <= 0 )); then
        log_info "Skipping AF2-multimer - all $n_total binders complete (--resume)"
        run_af2_multimer=false
    elif (( n_done > 0 )); then
        log_info "Resuming AF2-multimer: $n_remaining of $n_total binders remaining"
        af2_multimer_input_dir=$(prepare_resume_subset "$multimer_msa_dir" ".a3m" "$done_file")
    fi
    rm -f "$done_file"
fi

if [[ "$run_af2_multimer" == "true" ]]; then
    log_info "Starting AF2-multimer predictions..."
    conda activate colabfold

    if ! parallel_colabfold_batch "af2_multimer" \
        "$af2_multimer_input_dir" \
        "$pred_dir/af2_multimer" \
        --num-recycle "$num_recycle" \
        --recycle-early-stop-tolerance "$recycle_early_stop_tol" \
        --num-seeds "$num_seeds" \
        --num-models "$num_models" \
        --random-seed "$cf_random_seed" \
        --model-type alphafold2_multimer_v3 \
        --calc-extra-ptm \
        --rank iptm; then
        log_error "AF2-multimer predictions failed"
        pipeline_errors+=("af2_multimer")
    fi

    flatten_gpu_dirs "$pred_dir/af2_multimer" "$af2_multimer_input_dir"

    # Cleanup resume subset if one was created
    [[ "$af2_multimer_input_dir" != "$multimer_msa_dir" ]] && rm -rf "$af2_multimer_input_dir"

    conda deactivate
fi

# --- 2. ColabFold AF2-ptm on complex ---
run_af2_ptm=true
af2_ptm_input_dir="$multimer_msa_dir"

if [[ "$resume" == "true" ]]; then
    flatten_gpu_dirs "$pred_dir/af2_ptm_complex" "$multimer_msa_dir"

    done_file=$(mktemp)
    get_done_binders_af2 "$pred_dir/af2_ptm_complex" > "$done_file"
    n_done=$(wc -l < "$done_file")
    n_total=$(find "$multimer_msa_dir" -maxdepth 1 -name '*.a3m' -type f 2>/dev/null | wc -l)
    n_remaining=$((n_total - n_done))

    if (( n_remaining <= 0 )); then
        log_info "Skipping AF2-ptm complex - all $n_total binders complete (--resume)"
        run_af2_ptm=false
    elif (( n_done > 0 )); then
        log_info "Resuming AF2-ptm complex: $n_remaining of $n_total binders remaining"
        af2_ptm_input_dir=$(prepare_resume_subset "$multimer_msa_dir" ".a3m" "$done_file")
    fi
    rm -f "$done_file"
fi

if [[ "$run_af2_ptm" == "true" ]]; then
    log_info "Starting AF2-ptm complex predictions..."
    conda activate colabfold

    if ! parallel_colabfold_batch "af2_ptm_complex" \
        "$af2_ptm_input_dir" \
        "$pred_dir/af2_ptm_complex" \
        --num-recycle "$num_recycle" \
        --recycle-early-stop-tolerance "$recycle_early_stop_tol" \
        --num-seeds "$num_seeds" \
        --num-models "$num_models" \
        --random-seed "$cf_random_seed" \
        --model-type alphafold2_ptm \
        --calc-extra-ptm \
        --rank iptm; then
        log_error "AF2-ptm complex predictions failed"
        pipeline_errors+=("af2_ptm_complex")
    fi

    flatten_gpu_dirs "$pred_dir/af2_ptm_complex" "$af2_ptm_input_dir"

    [[ "$af2_ptm_input_dir" != "$multimer_msa_dir" ]] && rm -rf "$af2_ptm_input_dir"

    conda deactivate
fi

# --- 3. Chai-1 with MSAs ---
run_chai=true
chai_input_dir="$chai_fasta_dir"

if [[ "$resume" == "true" ]]; then
    done_file=$(mktemp)
    get_done_binders_chai "$pred_dir/chai" > "$done_file"
    n_done=$(wc -l < "$done_file")
    n_total=$(find "$chai_fasta_dir" -maxdepth 1 -name '*.fasta' -type f 2>/dev/null | wc -l)
    n_remaining=$((n_total - n_done))

    if (( n_remaining <= 0 )); then
        log_info "Skipping Chai-1 - all $n_total binders complete (--resume)"
        run_chai=false
    elif (( n_done > 0 )); then
        log_info "Resuming Chai-1: $n_remaining of $n_total binders remaining"
        chai_input_dir=$(prepare_resume_subset "$chai_fasta_dir" ".fasta" "$done_file")
    fi
    rm -f "$done_file"
fi

if [[ "$run_chai" == "true" ]]; then
    log_info "Starting Chai-1 predictions..."
    conda activate chai

    if ! run_chai_predictions \
        "$chai_input_dir" \
        "$chai_msa_dir" \
        "$pred_dir/chai" \
        "$chai_diffusion_samples" \
        "$chai_num_trunk_recycles" \
        "$chai_num_diffn_timesteps"; then
        log_error "Chai-1 predictions failed"
        pipeline_errors+=("chai")
    fi

    [[ "$chai_input_dir" != "$chai_fasta_dir" ]] && rm -rf "$chai_input_dir"

    conda deactivate
fi

# --- 4. Boltz-2 with per-chain MSAs ---
run_boltz=true
boltz_input_dir="$boltz_yaml_dir"

if [[ "$resume" == "true" ]]; then
    done_file=$(mktemp)
    get_done_binders_boltz "$pred_dir/boltz" > "$done_file"
    n_done=$(wc -l < "$done_file")
    n_total=$(find "$boltz_yaml_dir" -maxdepth 1 -name '*.yaml' -type f 2>/dev/null | wc -l)
    n_remaining=$((n_total - n_done))

    if (( n_remaining <= 0 )); then
        log_info "Skipping Boltz-2 - all $n_total binders complete (--resume)"
        run_boltz=false
    elif (( n_done > 0 )); then
        log_info "Resuming Boltz-2: $n_remaining of $n_total binders remaining"
        boltz_input_dir=$(prepare_resume_subset "$boltz_yaml_dir" ".yaml" "$done_file")
    fi
    rm -f "$done_file"
fi

if [[ "$run_boltz" == "true" ]]; then
    log_info "Starting Boltz-2 predictions..."
    conda activate boltz

    if ! boltz predict "$boltz_input_dir" \
        --out_dir "$pred_dir/boltz" \
        --output_format pdb \
        --recycling_steps "$boltz_recycling_steps" \
        --sampling_steps "$boltz_sampling_steps" \
        --diffusion_samples "$boltz_diffusion_samples" \
        --num_workers "$WORKERS_PER_GPU" \
        --devices "$NGPUS" \
        $use_potentials_flag \
        --override 2>"$pred_dir/boltz/boltz_stderr.log"; then
        log_error "Boltz-2 predictions failed (see $pred_dir/boltz/boltz_stderr.log)"
        pipeline_errors+=("boltz")
    fi

    [[ "$boltz_input_dir" != "$boltz_yaml_dir" ]] && rm -rf "$boltz_input_dir"

    conda deactivate
fi

# --- 5. Binder-only predictions ---
run_binder=true

if [[ "$binder_msa_mode" == "true" ]]; then
    # Binder uses AF2-ptm → check .done.txt
    binder_input_dir="$binder_msa_dir"

    if [[ "$resume" == "true" ]]; then
        flatten_gpu_dirs "$pred_dir/binder" "$binder_msa_dir"

        done_file=$(mktemp)
        get_done_binders_af2 "$pred_dir/binder" > "$done_file"
        n_done=$(wc -l < "$done_file")
        n_total=$(find "$binder_msa_dir" -maxdepth 1 -name '*.a3m' -type f 2>/dev/null | wc -l)
        n_remaining=$((n_total - n_done))

        if (( n_remaining <= 0 )); then
            log_info "Skipping binder-only (AF2-ptm) - all $n_total binders complete (--resume)"
            run_binder=false
        elif (( n_done > 0 )); then
            log_info "Resuming binder-only (AF2-ptm): $n_remaining of $n_total binders remaining"
            binder_input_dir=$(prepare_resume_subset "$binder_msa_dir" ".a3m" "$done_file")
        fi
        rm -f "$done_file"
    fi

    if [[ "$run_binder" == "true" ]]; then
        log_info "Starting binder-only (AF2-ptm) predictions..."
        conda activate colabfold

        if ! parallel_colabfold_batch "binder" \
            "$binder_input_dir" \
            "$pred_dir/binder" \
            --num-recycle "$num_recycle" \
            --recycle-early-stop-tolerance "$recycle_early_stop_tol" \
            --num-seeds "$num_seeds" \
            --num-models "$num_models" \
            --random-seed "$cf_random_seed" \
            --model-type alphafold2_ptm \
            --rank ptm; then
            log_error "Binder-only (AF2-ptm) predictions failed"
            pipeline_errors+=("binder")
        fi

        flatten_gpu_dirs "$pred_dir/binder" "$binder_input_dir"

        [[ "$binder_input_dir" != "$binder_msa_dir" ]] && rm -rf "$binder_input_dir"

        conda deactivate
    fi
else
    # Binder uses Chai → check subdirectories
    binder_chai_input_dir="$binder_fasta_dir"

    if [[ "$resume" == "true" ]]; then
        done_file=$(mktemp)
        get_done_binders_chai "$pred_dir/binder" > "$done_file"
        n_done=$(wc -l < "$done_file")
        n_total=$(find "$binder_fasta_dir" -maxdepth 1 -name '*.fasta' -type f 2>/dev/null | wc -l)
        n_remaining=$((n_total - n_done))

        if (( n_remaining <= 0 )); then
            log_info "Skipping binder-only (Chai) - all $n_total binders complete (--resume)"
            run_binder=false
        elif (( n_done > 0 )); then
            log_info "Resuming binder-only (Chai): $n_remaining of $n_total binders remaining"
            binder_chai_input_dir=$(prepare_resume_subset "$binder_fasta_dir" ".fasta" "$done_file")
        fi
        rm -f "$done_file"
    fi

    if [[ "$run_binder" == "true" ]]; then
        log_info "Starting binder-only (Chai) predictions..."
        conda activate chai

        if ! run_chai_predictions \
            "$binder_chai_input_dir" \
            "" \
            "$pred_dir/binder" \
            "$chai_diffusion_samples" \
            "$chai_num_trunk_recycles" \
            "$chai_num_diffn_timesteps"; then
            log_error "Binder-only (Chai) predictions failed"
            pipeline_errors+=("binder")
        fi

        [[ "$binder_chai_input_dir" != "$binder_fasta_dir" ]] && rm -rf "$binder_chai_input_dir"

        conda deactivate
    fi
fi

#===============================================================================
# Distance Calculation
#===============================================================================

if [[ "$resume" == "true" ]] && [[ -f "$project/distances.csv" ]]; then
    log_info "Skipping distance calculation - distances.csv already exists (--resume mode)"
else
    log_stage "Calculating Distances"

    conda activate chai  # Using chai env which has biopython

    python $SCRIPT_DIR/calculate_distances.py \
        --predictions_dir "$pred_dir" \
        --residue_csv "$residue_csv" \
        --csv "$csv" \
        --output_csv "$project/distances.csv" \
        --attempt 1

    conda deactivate

    log_info "Distance calculations saved to: $project/distances.csv"
fi

#===============================================================================
# Filter and Select Structures
#===============================================================================

#===============================================================================
# Initial Filter (Attempt 1)
#===============================================================================

if [[ "$resume" == "true" ]] && [[ -f "$project/summary_attempt1.csv" ]]; then
    log_info "Skipping initial filtering - summary_attempt1.csv already exists (--resume mode)"
else
    log_stage "Initial Filtering (Attempt 1)"

    conda activate chai  # Using chai env for Python

    # First pass filtering to identify failures
    python $SCRIPT_DIR/filter_and_select.py \
        --distances_csv "$project/distances.csv" \
        --predictions_dir "$pred_dir" \
        --cutoff "$distance_cutoff" \
        --output_dir "$project/filtered_structures_attempt1" \
        --summary_csv "$project/summary_attempt1.csv" \
        --binder_msa "$binder_msa_mode"

    conda deactivate
fi

#===============================================================================
# Check for Failures and Prepare Retry
#===============================================================================

# Count failures (excluding binder method)
n_failures=$(awk -F',' 'NR>1 && $8=="failed" && $2!="binder" {count++} END {print count+0}' "$project/summary_attempt1.csv")

log_info "Attempt 1 complete. Found $n_failures failures requiring retry."

if (( n_failures > 0 )); then
    log_stage "Preparing Retry Inputs (Attempt 2)"
    
    conda activate chai
    
    # Prepare retry inputs
    python $SCRIPT_DIR/prepare_retry_inputs.py \
        --summary_csv "$project/summary_attempt1.csv" \
        --project_dir "$project" \
        --retry_inputs_dir "$project/retry_inputs"
    
    conda deactivate
    
    #===============================================================================
    # Run Retry Predictions (Attempt 2)
    #===============================================================================
    
    log_stage "Running Retry Predictions (Attempt 2)"
    
    pred_dir_retry="$project/predictions_retry"
    mkdir -p "$pred_dir_retry"/{af2_multimer,af2_ptm_complex,chai,boltz}
    
    retry_msa_dir="$project/retry_inputs/colabfold/msa_multimer"
    retry_chai_fasta_dir="$project/retry_inputs/chai/fasta_multimer"
    retry_boltz_yaml_dir="$project/retry_inputs/boltz/yaml"
    
    # --- Retry AF2-multimer ---
    if [[ -d "$retry_msa_dir" ]] && ls "$retry_msa_dir"/*.a3m 1>/dev/null 2>&1; then
        run_retry_af2=true
        retry_af2_input_dir="$retry_msa_dir"

        if [[ "$resume" == "true" ]]; then
            flatten_gpu_dirs "$pred_dir_retry/af2_multimer" "$retry_msa_dir"

            done_file=$(mktemp)
            get_done_binders_af2 "$pred_dir_retry/af2_multimer" > "$done_file"
            n_done=$(wc -l < "$done_file")
            n_total=$(find "$retry_msa_dir" -maxdepth 1 -name '*.a3m' -type f 2>/dev/null | wc -l)
            n_remaining=$((n_total - n_done))

            if (( n_remaining <= 0 )); then
                log_info "Skipping AF2-multimer retry - all $n_total binders complete (--resume)"
                run_retry_af2=false
            elif (( n_done > 0 )); then
                log_info "Resuming AF2-multimer retry: $n_remaining of $n_total binders remaining"
                retry_af2_input_dir=$(prepare_resume_subset "$retry_msa_dir" ".a3m" "$done_file")
            fi
            rm -f "$done_file"
        fi

        if [[ "$run_retry_af2" == "true" ]]; then
            log_info "Retrying AF2-multimer predictions..."
            conda activate colabfold
            
            log_info "ColabFold retry random seed: $cf_retry_seed"
            
            if ! parallel_colabfold_batch "af2_multimer_retry" \
                "$retry_af2_input_dir" \
                "$pred_dir_retry/af2_multimer" \
                --num-recycle "$num_recycle" \
                --recycle-early-stop-tolerance "$recycle_early_stop_tol" \
                --num-seeds "$num_seeds" \
                --num-models "$num_models" \
                --random-seed "$cf_retry_seed" \
                --model-type alphafold2_multimer_v3 \
                --calc-extra-ptm \
                --rank iptm; then
                log_error "AF2-multimer retry failed"
                pipeline_errors+=("af2_multimer_retry")
            fi

            flatten_gpu_dirs "$pred_dir_retry/af2_multimer" "$retry_af2_input_dir"

            [[ "$retry_af2_input_dir" != "$retry_msa_dir" ]] && rm -rf "$retry_af2_input_dir"
            
            conda deactivate
        fi
    fi
    
    # --- Retry AF2-ptm complex ---
    if [[ -d "$retry_msa_dir" ]] && ls "$retry_msa_dir"/*.a3m 1>/dev/null 2>&1; then
        run_retry_ptm=true
        retry_ptm_input_dir="$retry_msa_dir"

        if [[ "$resume" == "true" ]]; then
            flatten_gpu_dirs "$pred_dir_retry/af2_ptm_complex" "$retry_msa_dir"

            done_file=$(mktemp)
            get_done_binders_af2 "$pred_dir_retry/af2_ptm_complex" > "$done_file"
            n_done=$(wc -l < "$done_file")
            n_total=$(find "$retry_msa_dir" -maxdepth 1 -name '*.a3m' -type f 2>/dev/null | wc -l)
            n_remaining=$((n_total - n_done))

            if (( n_remaining <= 0 )); then
                log_info "Skipping AF2-ptm complex retry - all $n_total binders complete (--resume)"
                run_retry_ptm=false
            elif (( n_done > 0 )); then
                log_info "Resuming AF2-ptm complex retry: $n_remaining of $n_total binders remaining"
                retry_ptm_input_dir=$(prepare_resume_subset "$retry_msa_dir" ".a3m" "$done_file")
            fi
            rm -f "$done_file"
        fi

        if [[ "$run_retry_ptm" == "true" ]]; then
            log_info "Retrying AF2-ptm complex predictions..."
            conda activate colabfold
            
            if ! parallel_colabfold_batch "af2_ptm_complex_retry" \
                "$retry_ptm_input_dir" \
                "$pred_dir_retry/af2_ptm_complex" \
                --num-recycle "$num_recycle" \
                --recycle-early-stop-tolerance "$recycle_early_stop_tol" \
                --num-seeds "$num_seeds" \
                --num-models "$num_models" \
                --random-seed "$cf_retry_seed" \
                --model-type alphafold2_ptm \
                --calc-extra-ptm \
                --rank iptm; then
                log_error "AF2-ptm complex retry failed"
                pipeline_errors+=("af2_ptm_complex_retry")
            fi

            flatten_gpu_dirs "$pred_dir_retry/af2_ptm_complex" "$retry_ptm_input_dir"

            [[ "$retry_ptm_input_dir" != "$retry_msa_dir" ]] && rm -rf "$retry_ptm_input_dir"
            
            conda deactivate
        fi
    fi
    
    # --- Retry Chai-1 ---
    if [[ -d "$retry_chai_fasta_dir" ]] && ls "$retry_chai_fasta_dir"/*.fasta 1>/dev/null 2>&1; then
        run_retry_chai=true
        retry_chai_input_dir="$retry_chai_fasta_dir"

        if [[ "$resume" == "true" ]]; then
            done_file=$(mktemp)
            get_done_binders_chai "$pred_dir_retry/chai" > "$done_file"
            n_done=$(wc -l < "$done_file")
            n_total=$(find "$retry_chai_fasta_dir" -maxdepth 1 -name '*.fasta' -type f 2>/dev/null | wc -l)
            n_remaining=$((n_total - n_done))

            if (( n_remaining <= 0 )); then
                log_info "Skipping Chai-1 retry - all $n_total binders complete (--resume)"
                run_retry_chai=false
            elif (( n_done > 0 )); then
                log_info "Resuming Chai-1 retry: $n_remaining of $n_total binders remaining"
                retry_chai_input_dir=$(prepare_resume_subset "$retry_chai_fasta_dir" ".fasta" "$done_file")
            fi
            rm -f "$done_file"
        fi

        if [[ "$run_retry_chai" == "true" ]]; then
            log_info "Retrying Chai-1 predictions..."
            conda activate chai
            
            # Note: Chai MSA directory is reused (hash-based lookup)
            if ! run_chai_predictions \
                "$retry_chai_input_dir" \
                "$chai_msa_dir" \
                "$pred_dir_retry/chai" \
                "$chai_diffusion_samples" \
                "$chai_num_trunk_recycles" \
                "$chai_num_diffn_timesteps"; then
                log_error "Chai-1 retry failed"
                pipeline_errors+=("chai_retry")
            fi

            [[ "$retry_chai_input_dir" != "$retry_chai_fasta_dir" ]] && rm -rf "$retry_chai_input_dir"
            
            conda deactivate
        fi
    fi
    
    # --- Retry Boltz-2 ---
    if [[ -d "$retry_boltz_yaml_dir" ]] && ls "$retry_boltz_yaml_dir"/*.yaml 1>/dev/null 2>&1; then
        run_retry_boltz=true
        retry_boltz_input_dir="$retry_boltz_yaml_dir"

        if [[ "$resume" == "true" ]]; then
            done_file=$(mktemp)
            get_done_binders_boltz "$pred_dir_retry/boltz" > "$done_file"
            n_done=$(wc -l < "$done_file")
            n_total=$(find "$retry_boltz_yaml_dir" -maxdepth 1 -name '*.yaml' -type f 2>/dev/null | wc -l)
            n_remaining=$((n_total - n_done))

            if (( n_remaining <= 0 )); then
                log_info "Skipping Boltz-2 retry - all $n_total binders complete (--resume)"
                run_retry_boltz=false
            elif (( n_done > 0 )); then
                log_info "Resuming Boltz-2 retry: $n_remaining of $n_total binders remaining"
                retry_boltz_input_dir=$(prepare_resume_subset "$retry_boltz_yaml_dir" ".yaml" "$done_file")
            fi
            rm -f "$done_file"
        fi

        if [[ "$run_retry_boltz" == "true" ]]; then
            log_info "Retrying Boltz-2 predictions..."
            conda activate boltz
            
            if ! boltz predict "$retry_boltz_input_dir" \
                --out_dir "$pred_dir_retry/boltz" \
                --output_format pdb \
                --recycling_steps "$boltz_recycling_steps" \
                --sampling_steps "$boltz_sampling_steps" \
                --diffusion_samples "$boltz_diffusion_samples" \
                --num_workers "$WORKERS_PER_GPU" \
                --devices "$NGPUS" \
                $use_potentials_flag \
                --override 2>"$pred_dir_retry/boltz/boltz_stderr.log"; then
                log_error "Boltz-2 retry failed (see $pred_dir_retry/boltz/boltz_stderr.log)"
                pipeline_errors+=("boltz_retry")
            fi

            [[ "$retry_boltz_input_dir" != "$retry_boltz_yaml_dir" ]] && rm -rf "$retry_boltz_input_dir"
            
            conda deactivate
        fi
    fi
    
    #===============================================================================
    # Calculate Distances for Retry (Attempt 2)
    #===============================================================================
    
    if [[ "$resume" == "true" ]] && [[ -f "$project/distances_retry.csv" ]]; then
        log_info "Skipping retry distance calculation - distances_retry.csv already exists (--resume mode)"
    else
        log_stage "Calculating Distances for Retry (Attempt 2)"
        
        conda activate chai
        
        python $SCRIPT_DIR/calculate_distances.py \
            --predictions_dir "$pred_dir_retry" \
            --residue_csv "$residue_csv" \
            --csv "$csv" \
            --output_csv "$project/distances_retry.csv" \
            --attempt 2
        
        conda deactivate
        
        log_info "Retry distance calculations saved to: $project/distances_retry.csv"
    fi
    
    #===============================================================================
    # Final Filter (Merge Attempt 1 + Attempt 2)
    #===============================================================================
    
    if [[ "$resume" == "true" ]] && [[ -d "$project/filtered_structures" ]] && [[ -f "$project/summary.csv" ]]; then
        log_info "Skipping final filtering - filtered_structures and summary.csv already exist (--resume mode)"
    else
        log_stage "Final Filtering (Merging Attempts)"
        
        conda activate chai
        
        # Remove attempt 1 filtered structures (will be replaced by final)
        rm -rf "$project/filtered_structures_attempt1"
        
        python $SCRIPT_DIR/filter_and_select.py \
            --distances_csv "$project/distances.csv" \
            --predictions_dir "$pred_dir" \
            --cutoff "$distance_cutoff" \
            --output_dir "$project/filtered_structures" \
            --summary_csv "$project/summary.csv" \
            --binder_msa "$binder_msa_mode" \
            --distances_csv_retry "$project/distances_retry.csv" \
            --predictions_dir_retry "$pred_dir_retry"
        
        conda deactivate
    fi

else
    # No failures - just rename attempt 1 outputs to final
    if [[ "$resume" == "true" ]] && [[ -d "$project/filtered_structures" ]] && [[ -f "$project/summary.csv" ]]; then
        log_info "Skipping final filtering - filtered_structures and summary.csv already exist (--resume mode)"
    else
        log_info "No failures detected. Using attempt 1 results as final."
        mv "$project/filtered_structures_attempt1" "$project/filtered_structures"
        mv "$project/summary_attempt1.csv" "$project/summary.csv"
    fi
fi

#===============================================================================
# Completion
#===============================================================================

log_stage "Pipeline Complete"

log_info "Results:"
log_info "  - Predictions: $pred_dir"
log_info "  - Distances: $project/distances.csv"
log_info "  - Filtered structures: $project/filtered_structures/"
log_info "  - Summary: $project/summary.csv"

# Create completion marker
cat > "$project/.prediction_complete" <<EOF
# Structure Prediction Complete
# Generated: $(date)
# Residue CSV: $residue_csv
# Distance cutoff: $distance_cutoff
# Binder MSA mode: $binder_msa_mode
# --- ColabFold ---
# Num recycle: $num_recycle
# Recycle early stop tolerance: $recycle_early_stop_tol
# Num seeds: $num_seeds
# Num models: $num_models
# Random seed: $cf_random_seed
# --- Chai ---
# Num trunk recycles: $chai_num_trunk_recycles
# Num diffusion timesteps: $chai_num_diffn_timesteps
# Diffusion samples: $chai_diffusion_samples
# --- Boltz ---
# Recycling steps: $boltz_recycling_steps
# Sampling steps: $boltz_sampling_steps
# Diffusion samples: $boltz_diffusion_samples
# Use potentials: $use_potentials
EOF

log_info "Completion marker: $project/.prediction_complete"

# Report pipeline errors
if (( ${#pipeline_errors[@]} > 0 )); then
    log_warn "Pipeline completed with errors in: ${pipeline_errors[*]}"
    exit 1
fi