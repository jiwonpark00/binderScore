#!/usr/bin/env bash
#===============================================================================
# binderScore.sh - Full binder screening pipeline
#
# Runs all three stages sequentially:
#   Stage 1: MSA Generation          (generate_msas.sh)
#   Stage 2: Structure Prediction    (predict_structure.sh)
#   Stage 3: Metrics & Scoring       (score.sh)
#
# Each stage checks completion markers and skips if already done.
# To re-run a stage, delete the corresponding marker file:
#   Stage 1: $project/.msa_complete
#   Stage 2: $project/.prediction_complete
#   Stage 3: $project/.scoring_complete
#===============================================================================

set -euo pipefail

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

Full binder screening pipeline: MSA generation → structure prediction → scoring.

Required flags:
  --csv                        CSV file with name,target,binder columns
  --project                    Output project directory
  --target-dir                 Target directory (must contain targets.fasta and /pdb)
  --binder-msa                 Generate binder MSAs (true/false)
  --residue-csv                CSV with target_name,residue_number columns
  --distance-cutoff            Distance cutoff in Angstroms

Conditional flags (Stage 1 — required only if --binder-msa true):
  --msa-method                 MSA generation method: 'server' or 'mem'
  --db                         ColabFold database directory (required only if --msa-method mem)

ColabFold (AF2) optional flags (Stage 2):
  --num-recycle                Number of recycles for AF2 (default: 2)
  --recycle-early-stop-tol     Early stop tolerance for AF2 recycles (default: 0.5)
  --num-seeds                  Number of seeds per model (default: 1)
  --num-models                 Number of AF2 models to run, 1-5 (default: 5)

Chai optional flags (Stage 2):
  --chai-num-trunk-recycles    Number of trunk recycles (default: 3)
  --chai-num-diffn-timesteps   Number of diffusion timesteps (default: 200)
  --chai-diffusion-samples     Number of diffusion samples (default: 5)

Boltz optional flags (Stage 2):
  --boltz-recycling-steps      Number of recycling steps (default: 3)
  --boltz-sampling-steps       Number of sampling steps (default: 200)
  --boltz-diffusion-samples    Number of diffusion samples (default: 5)
  --use-potentials             Use inference-time potentials (default: true)

Scoring optional flags (Stage 3):
  --pae-threshold              PAE threshold (default: 10.0)
  --pde-threshold              PDE threshold for Boltz (default: 2.0)
  --interface-threshold        Interface distance threshold (default: 4.0)
  --pdockq-threshold           pDockQ threshold (default: 8.0)

Options:
  --help, -h                   Show this help message and exit

Examples:
  # Minimal (target-only MSAs, all defaults)
  $0 --csv input.csv --project my_project/ --target-dir targets/ \\
     --binder-msa false --residue-csv residues.csv --distance-cutoff 20

  # Binder MSAs via ColabFold server
  $0 --csv input.csv --project my_project/ --target-dir targets/ \\
     --binder-msa true --msa-method server \\
     --residue-csv residues.csv --distance-cutoff 20

  # Binder MSAs via local DB + custom prediction/scoring params
  $0 --csv input.csv --project my_project/ --target-dir targets/ \\
     --binder-msa true --msa-method mem --db /data/colabfold_db \\
     --residue-csv residues.csv --distance-cutoff 20 \\
     --num-recycle 3 --boltz-diffusion-samples 10 --pae-threshold 8.0
EOF
}

#===============================================================================
# Argument parsing
#===============================================================================

# --- Always required ---
csv=""
project=""
target_dir=""
binder_msa=""
residue_csv=""
distance_cutoff=""

# --- Conditional (Stage 1) ---
msa_method=""
db=""

# --- ColabFold / AF2 defaults (Stage 2) ---
num_recycle="2"
recycle_early_stop_tol="0.5"
num_seeds="1"
num_models="5"

# --- Chai defaults (Stage 2) ---
chai_num_trunk_recycles="3"
chai_num_diffn_timesteps="200"
chai_diffusion_samples="5"

# --- Boltz defaults (Stage 2) ---
boltz_recycling_steps="3"
boltz_sampling_steps="200"
boltz_diffusion_samples="5"
use_potentials="true"

# --- Scoring defaults (Stage 3) ---
pae_threshold="10.0"
pde_threshold="2.0"
interface_threshold="4.0"
pdockq_threshold="8.0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        # Required
        --csv)                         csv="$2"; shift 2 ;;
        --project)                     project="$2"; shift 2 ;;
        --target-dir)                  target_dir="$2"; shift 2 ;;
        --binder-msa)                  binder_msa="$2"; shift 2 ;;
        --residue-csv)                 residue_csv="$2"; shift 2 ;;
        --distance-cutoff)             distance_cutoff="$2"; shift 2 ;;
        # Conditional (Stage 1)
        --msa-method)                  msa_method="$2"; shift 2 ;;
        --db)                          db="$2"; shift 2 ;;
        # ColabFold / AF2 (Stage 2)
        --num-recycle)                 num_recycle="$2"; shift 2 ;;
        --recycle-early-stop-tol)      recycle_early_stop_tol="$2"; shift 2 ;;
        --num-seeds)                   num_seeds="$2"; shift 2 ;;
        --num-models)                  num_models="$2"; shift 2 ;;
        # Chai (Stage 2)
        --chai-num-trunk-recycles)     chai_num_trunk_recycles="$2"; shift 2 ;;
        --chai-num-diffn-timesteps)    chai_num_diffn_timesteps="$2"; shift 2 ;;
        --chai-diffusion-samples)      chai_diffusion_samples="$2"; shift 2 ;;
        # Boltz (Stage 2)
        --boltz-recycling-steps)       boltz_recycling_steps="$2"; shift 2 ;;
        --boltz-sampling-steps)        boltz_sampling_steps="$2"; shift 2 ;;
        --boltz-diffusion-samples)     boltz_diffusion_samples="$2"; shift 2 ;;
        --use-potentials)              use_potentials="$2"; shift 2 ;;
        # Scoring (Stage 3)
        --pae-threshold)               pae_threshold="$2"; shift 2 ;;
        --pde-threshold)               pde_threshold="$2"; shift 2 ;;
        --interface-threshold)         interface_threshold="$2"; shift 2 ;;
        --pdockq-threshold)            pdockq_threshold="$2"; shift 2 ;;
        # Help
        --help|-h)                     show_help; exit 0 ;;
        *)                             log_error "Unknown flag: $1"; show_help; exit 1 ;;
    esac
done

#===============================================================================
# Validate required flags
#===============================================================================

missing=()
[[ -z "$csv" ]]             && missing+=("--csv")
[[ -z "$project" ]]         && missing+=("--project")
[[ -z "$target_dir" ]]      && missing+=("--target-dir")
[[ -z "$binder_msa" ]]      && missing+=("--binder-msa")
[[ -z "$residue_csv" ]]     && missing+=("--residue-csv")
[[ -z "$distance_cutoff" ]] && missing+=("--distance-cutoff")

if [[ ${#missing[@]} -gt 0 ]]; then
    log_error "Missing required flags: ${missing[*]}"
    show_help
    exit 1
fi

# Validate binder_msa value
if [[ "$binder_msa" != "true" && "$binder_msa" != "false" ]]; then
    log_error "--binder-msa must be 'true' or 'false', got: $binder_msa"
    exit 1
fi

# Conditional validation: --msa-method required when --binder-msa true
if [[ "$binder_msa" == "true" ]]; then
    if [[ -z "$msa_method" ]]; then
        log_error "--msa-method is required when --binder-msa true"
        log_error "Use --msa-method server or --msa-method mem"
        show_help
        exit 1
    fi
    if [[ "$msa_method" != "server" && "$msa_method" != "mem" ]]; then
        log_error "--msa-method must be 'server' or 'mem', got: $msa_method"
        exit 1
    fi
    if [[ "$msa_method" == "mem" && -z "$db" ]]; then
        log_error "--db is required when --msa-method mem"
        exit 1
    fi
fi

# Validate use_potentials value
if [[ "$use_potentials" != "true" && "$use_potentials" != "false" ]]; then
    log_error "--use-potentials must be 'true' or 'false', got: $use_potentials"
    exit 1
fi

# Resolve script directory (assumes all scripts are co-located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#===============================================================================
# Pipeline configuration summary
#===============================================================================

log_stage "Binder Screening Pipeline"
log_info "Configuration:"
log_info "  CSV:                    $csv"
log_info "  Project:                $project"
log_info "  Target dir:             $target_dir"
log_info "  Binder MSA:             $binder_msa"
if [[ "$binder_msa" == "true" ]]; then
    log_info "  MSA method:             $msa_method"
    [[ -n "$db" ]] && log_info "  Database:               $db"
fi
log_info "  Residue CSV:            $residue_csv"
log_info "  Distance cutoff:        $distance_cutoff Å"
log_info "--- ColabFold (AF2) ---"
log_info "  Num recycle:            $num_recycle"
log_info "  Recycle early stop tol: $recycle_early_stop_tol"
log_info "  Num seeds:              $num_seeds"
log_info "  Num models:             $num_models"
log_info "--- Chai ---"
log_info "  Num trunk recycles:     $chai_num_trunk_recycles"
log_info "  Num diffn timesteps:    $chai_num_diffn_timesteps"
log_info "  Diffusion samples:      $chai_diffusion_samples"
log_info "--- Boltz ---"
log_info "  Recycling steps:        $boltz_recycling_steps"
log_info "  Sampling steps:         $boltz_sampling_steps"
log_info "  Diffusion samples:      $boltz_diffusion_samples"
log_info "  Use potentials:         $use_potentials"
log_info "--- Scoring ---"
log_info "  PAE threshold:          $pae_threshold"
log_info "  PDE threshold:          $pde_threshold"
log_info "  Interface threshold:    $interface_threshold"
log_info "  pDockQ threshold:       $pdockq_threshold"

#===============================================================================
# Stage 1: MSA Generation
#===============================================================================

log_stage "Starting Stage 1: MSA Generation"

# Build Stage 1 command — only pass conditional flags when relevant
stage1_cmd=(
    bash "$SCRIPT_DIR/binderScore/generate_msas/generate_msas.sh"
    --csv "$csv"
    --project "$project"
    --target-dir "$target_dir"
    --binder-msa "$binder_msa"
)

if [[ "$binder_msa" == "true" ]]; then
    stage1_cmd+=(--msa-method "$msa_method")
    if [[ -n "$db" ]]; then
        stage1_cmd+=(--db "$db")
    fi
fi

"${stage1_cmd[@]}"

log_info "Stage 1 complete."

#===============================================================================
# Stage 2: Structure Prediction
#===============================================================================

log_stage "Starting Stage 2: Structure Prediction"

bash "$SCRIPT_DIR/binderScore/predict_structure/predict_structure.sh" \
    --csv "$csv" \
    --project "$project" \
    --target-dir "$target_dir" \
    --residue-csv "$residue_csv" \
    --distance-cutoff "$distance_cutoff" \
    --binder-msa "$binder_msa" \
    --num-recycle "$num_recycle" \
    --recycle-early-stop-tol "$recycle_early_stop_tol" \
    --num-seeds "$num_seeds" \
    --num-models "$num_models" \
    --chai-num-trunk-recycles "$chai_num_trunk_recycles" \
    --chai-num-diffn-timesteps "$chai_num_diffn_timesteps" \
    --chai-diffusion-samples "$chai_diffusion_samples" \
    --boltz-recycling-steps "$boltz_recycling_steps" \
    --boltz-sampling-steps "$boltz_sampling_steps" \
    --boltz-diffusion-samples "$boltz_diffusion_samples" \
    --use-potentials "$use_potentials"

log_info "Stage 2 complete."

#===============================================================================
# Stage 3: Metrics & Scoring
#===============================================================================

log_stage "Starting Stage 3: Metrics & Scoring"

bash "$SCRIPT_DIR/binderScore/score/score.sh" \
    --project "$project" \
    --pae-threshold "$pae_threshold" \
    --pde-threshold "$pde_threshold" \
    --interface-threshold "$interface_threshold" \
    --pdockq-threshold "$pdockq_threshold"

log_info "Stage 3 complete."

#===============================================================================
# Pipeline complete
#===============================================================================

log_stage "Pipeline Complete"
log_info "All three stages finished successfully."
log_info "Project directory: $project"
log_info ""
log_info "Key outputs:"
if [[ -f "$project/binderFinalScore.csv" ]]; then
    log_info "  - Final scores:        $project/binderFinalScore.csv"
else
    log_info "  - Final scores:        $project/binderScore.csv"
fi
log_info "  - Rosetta metrics:     $project/binderRosettaScore.csv"
log_info "  - RMSD metrics:        $project/RMSD.csv"
log_info "  - Filtered structures: $project/filtered_structures/"
log_info "  - Relaxed structures:  $project/relaxed/"
log_info "  - Structure symlinks:  $project/symlinks/"
log_info "  - Per-method metrics:  $project/metrics/"