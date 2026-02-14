#!/usr/bin/env bash
#===============================================================================
# Stage 1: MSA Generation for binder screening pipeline
#
# Inputs (named flags):
#   --csv          - CSV file (name,target,binder)
#   --project      - Project directory
#   --target-dir   - Target directory (contains /pdb and /msa subdirectories)
#   --binder-msa   - Generate binder MSAs (true/false)
#   --msa-method   - MSA generation method: server or mem (required if binder-msa=true)
#   --db           - ColabFold database path (required if msa-method=mem)
#
# Outputs:
#   If binder_msa=false:
#     - Target MSAs in $target_dir/msa/
#     - Multimer MSAs in $project/colabfold/msa_multimer/ (template-based)
#     - Chai-1 parquet MSAs in $project/chai/msas/ (from target MSAs only)
#   If binder_msa=true:
#     - Target MSAs in $target_dir/msa/
#     - Binder MSAs in $project/colabfold/msa_binder/
#     - Multimer MSAs in $project/colabfold/msa_multimer/ (combined from target + binder)
#     - Chai-1 parquet MSAs in $project/chai/msas/ (from target + binder MSAs)
#   - Completion marker: $project/.msa_complete
#===============================================================================

set -euo pipefail
export DISABLE_PANDERA_IMPORT_WARNING=True
export PYTHONWARNINGS="ignore"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#===============================================================================
# Simple logging functions using echo
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

Stage 1: MSA Generation for binder screening pipeline

Required flags:
  --csv          CSV file with name,target,binder columns
  --project      Output project directory
  --target-dir   Target directory (must contain targets.fasta and /pdb)
  --binder-msa   Generate binder MSAs (true/false)

Conditional flags (required only if --binder-msa true):
  --msa-method   MSA generation method: 'server' or 'mem'
                   server = use ColabFold MSA server (no local DB needed)
                   mem    = use local DB loaded into memory via vmtouch
                            (requires --db and >=700GB RAM; falls back to server if not enough RAM)
  --db           ColabFold database directory (required only if --msa-method mem)

Options:
  --help, -h     Show this help message and exit

Examples:
  # Generate MSA for targets only and use as template to generate other MSAs
  $0 --csv input.csv --project project/ --target-dir targets/ --binder-msa false

  # Generate binder MSAs via ColabFold server, then combine with target MSAs
  $0 --csv input.csv --project project/ --target-dir targets/ --binder-msa true --msa-method server

  # Generate binder MSAs via local DB in memory
  $0 --csv input.csv --db /data/db --project project/ --target-dir targets/ --binder-msa true --msa-method mem
EOF
}

#===============================================================================
# Argument parsing
#===============================================================================

csv=""
db=""
project=""
target_dir=""
binder_msa=""
msa_method=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --csv)          csv="$2"; shift 2 ;;
        --db)           db="$2"; shift 2 ;;
        --project)      project="$2"; shift 2 ;;
        --target-dir)   target_dir="$2"; shift 2 ;;
        --binder-msa)   binder_msa="$2"; shift 2 ;;
        --msa-method)   msa_method="$2"; shift 2 ;;
        --help|-h)      show_help; exit 0 ;;
        *)              log_error "Unknown flag: $1"; show_help; exit 1 ;;
    esac
done

# Validate required flags (always required)
missing=()
[[ -z "$csv" ]]        && missing+=("--csv")
[[ -z "$project" ]]    && missing+=("--project")
[[ -z "$target_dir" ]] && missing+=("--target-dir")
[[ -z "$binder_msa" ]] && missing+=("--binder-msa")

if [[ ${#missing[@]} -gt 0 ]]; then
    log_error "Missing required flags: ${missing[*]}"
    show_help
    exit 1
fi

# Validate inputs
[[ -f "$csv" ]] || { log_error "CSV file not found: $csv"; exit 1; }
lastchar=$(tail -c 1 "$csv")
if [[ "$lastchar" != '' ]]; then
    printf '\n' >> "$csv"
fi

[[ -d "$target_dir" ]] || { log_error "Target directory not found: $target_dir"; exit 1; }

# Validate binder_msa flag
if [[ "$binder_msa" != "true" && "$binder_msa" != "false" ]]; then
    log_error "binder_msa must be 'true' or 'false', got: $binder_msa"
    exit 1
fi

# Conditional validation for binder_msa=true
if [[ "$binder_msa" == "true" ]]; then
    if [[ -z "$msa_method" ]]; then
        log_error "--msa-method is required when --binder-msa true"
        log_error "Use --msa-method server or --msa-method mem"
        show_help
        exit 1
    fi

    if [[ "$msa_method" != "server" && "$msa_method" != "mem" ]]; then
        log_error "msa_method must be 'server' or 'mem', got: $msa_method"
        exit 1
    fi

    if [[ "$msa_method" == "mem" ]]; then
        if [[ -z "$db" ]]; then
            log_error "--db is required when --msa-method mem"
            exit 1
        fi
        [[ -d "$db" ]] || { log_error "Database path not found: $db"; exit 1; }

        # Check available RAM (require >=700GB)
        total_ram_kb=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
        total_ram_gb=$(( total_ram_kb / 1024 / 1024 ))
        if (( total_ram_gb < 700 )); then
            log_warn "System has ${total_ram_gb}GB RAM, but msa-method=mem requires >=700GB"
            log_warn "Forcing msa-method=server"
            msa_method="server"
        else
            log_info "System has ${total_ram_gb}GB RAM - sufficient for msa-method=mem"
        fi
    fi
fi

# Create project directory
mkdir -p "$project"

# Check if already complete (only for binder_msa=true mode)
if [[ "$binder_msa" == "true" && -f "$project/.msa_complete" ]]; then
    log_warn "MSAs already generated (found .msa_complete marker)"
    log_info "To regenerate, delete: $project/.msa_complete"
    exit 0
fi

#===============================================================================
# Environment setup
#===============================================================================

log_stage "Stage 1: MSA Generation"
log_info "CSV: $csv"
log_info "Project: $project"
log_info "Target dir: $target_dir"
if [[ "$binder_msa" == "true" ]]; then
    log_info "Mode: Binder MSAs (full)"
    log_info "MSA method: $msa_method"
    [[ -n "$db" ]] && log_info "Database: $db"
else
    log_info "Mode: Target MSAs only (template-based)"
fi

source ~/miniforge3/etc/profile.d/conda.sh
shopt -s nullglob

# Detect GPUs
NGPUS=$(nvidia-smi -L 2>/dev/null | wc -l)
if (( NGPUS == 0 )); then
    log_error "No GPUs detected (nvidia-smi -L failed)"
    exit 1
fi
log_info "Detected $NGPUS GPUs"

CUDA_VISIBLE_DEVICES=$(printf '%s,' $(seq 0 $((NGPUS-1))))
CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES%,}
export CUDA_VISIBLE_DEVICES

#===============================================================================
# Dependency checks
#===============================================================================

check_command() {
    command -v "$1" >/dev/null 2>&1 || { 
        log_error "Missing dependency: $1"
        return 1
    }
}

check_command vmtouch || log_warn "vmtouch not found - DB loading will be skipped"
check_command jq || { log_error "jq required for JSON parsing"; exit 1; }

# Verify Python scripts based on mode
if [[ "$binder_msa" == "false" ]]; then
    # Target-only mode needs generate_multimer_a3m.py
    if [[ ! -f "$SCRIPT_DIR/generate_multimer_a3m.py" ]]; then
        log_error "generate_multimer_a3m.py is required for target-only mode"
        log_error "This script generates multimer MSAs from target templates"
        exit 1
    fi
elif [[ "$binder_msa" == "true" ]]; then
    # Binder mode needs combine_target_binder_a3m.py
    if [[ ! -f "$SCRIPT_DIR/combine_target_binder_a3m.py" ]]; then
        log_error "combine_target_binder_a3m.py is required for binder MSA mode"
        log_error "This script combines target + binder MSAs into multimer A3Ms"
        exit 1
    fi
fi

# Both modes need the new parallel conversion script
if [[ ! -f "$SCRIPT_DIR/a3m_to_pqt_parallel.py" ]]; then
    log_error "a3m_to_pqt_parallel.py is required for MSA conversion"
    log_error "Expected at: $SCRIPT_DIR/a3m_to_pqt_parallel.py"
    exit 1
fi

#===============================================================================
# Activate ColabFold environment
#===============================================================================

log_info "Activating ColabFold environment..."
conda activate colabfold

# Verify ColabFold tools
check_command mmseqs || { log_error "mmseqs not found in environment"; exit 1; }
check_command colabfold_search || { log_error "colabfold_search not found"; exit 1; }
check_command colabfold_batch || { log_error "colabfold_batch not found"; exit 1; }

#===============================================================================
# Convert CSV to FASTA
#===============================================================================

log_info "Converting CSV to FASTA inputs..."

colabfold_dir="$project/colabfold"
mkdir -p "$colabfold_dir"

msa_dir_multimer="$colabfold_dir/msa_multimer"
msa_dir_binder="$colabfold_dir/msa_binder"
mkdir -p "$msa_dir_multimer" "$msa_dir_binder"

# Multimer FASTA (target:binder)
colabfold_input_multimer="$colabfold_dir/colabfold_input_multimer.fasta"
awk -F',' 'NR>1 {print ">"$1"\n"$2":"$3}' "$csv" > "$colabfold_input_multimer"
n_multimer=$(grep -c "^>" "$colabfold_input_multimer" || echo 0)
log_info "Created multimer FASTA with $n_multimer sequences"

# Binder-only FASTA
colabfold_input_binder="$colabfold_dir/colabfold_input_binder.fasta"
awk -F',' 'NR>1 {print ">"$1"\n"$3}' "$csv" > "$colabfold_input_binder"
n_binder=$(grep -c "^>" "$colabfold_input_binder" || echo 0)
log_info "Created binder FASTA with $n_binder sequences"

if (( n_multimer == 0 )); then
    log_error "No sequences found in CSV. Check CSV format: name,target,binder"
    exit 1
fi

#===============================================================================
# SHARED FUNCTION: Generate Target MSAs
# This runs in BOTH modes to ensure target MSAs exist for Chai-1
#===============================================================================

generate_target_msas() {
    log_stage "Generating Target MSAs (Shared Step)"
    
    # Validate target directory structure
    targets_fasta="$target_dir/targets.fasta"
    if [[ ! -f "$targets_fasta" ]]; then
        log_error "targets.fasta not found: $targets_fasta"
        log_error "Expected: $target_dir/targets.fasta"
        log_error ""
        log_error "Create this file with target sequences:"
        log_error "  >AQP"
        log_error "  MSEQUENCE..."
        log_error "  >SLP"
        log_error "  MSEQUENCE..."
        exit 1
    fi
    
    # Create target MSA directory
    target_msa_dir="$target_dir/msa"
    mkdir -p "$target_msa_dir"
    
    #===========================================================================
    # Extract target names from targets.fasta
    #===========================================================================
    
    log_info "Extracting target names from targets.fasta..."
    target_names=()
    
    while IFS= read -r line; do
        # Check if line is a FASTA header
        if [[ "$line" =~ ^\>(.+)$ ]]; then
            target_name="${BASH_REMATCH[1]}"
            # Remove any whitespace or special characters
            target_name=$(echo "$target_name" | awk '{print $1}')
            target_names+=("$target_name")
            log_info "  Found target: $target_name"
        fi
    done < "$targets_fasta"
    
    if [[ ${#target_names[@]} -eq 0 ]]; then
        log_error "No target sequences found in $targets_fasta"
        log_error "Expected format:"
        log_error "  >TARGET_NAME"
        log_error "  SEQUENCE..."
        exit 1
    fi
    
    log_info "Total targets found: ${#target_names[@]}"
    
    #===========================================================================
    # Check which targets need MSAs generated
    #===========================================================================
    
    targets_to_generate=()
    targets_existing=()
    
    for target_name in "${target_names[@]}"; do
        if [[ -f "$target_msa_dir/${target_name}.a3m" ]]; then
            targets_existing+=("$target_name")
            log_info "  ✓ MSA exists: $target_name"
        else
            targets_to_generate+=("$target_name")
            log_info "  ✗ MSA missing: $target_name"
        fi
    done
    
    log_info "MSA status: ${#targets_existing[@]} existing, ${#targets_to_generate[@]} to generate"
    
    #===========================================================================
    # Generate missing target MSAs
    #===========================================================================
    
    if [[ ${#targets_to_generate[@]} -gt 0 ]]; then
        log_info "Generating MSAs for ${#targets_to_generate[@]} targets..."
        log_info "Command: colabfold_batch --msa-only $targets_fasta $target_msa_dir"
        
        colabfold_batch \
            --msa-only \
            "$targets_fasta" \
            "$target_msa_dir" \
            >> "$target_msa_dir/colabfold_batch.log" 2>&1
        
        # Validate MSA generation
        log_info "Validating MSA generation..."
        all_generated=true
        for target_name in "${targets_to_generate[@]}"; do
            if [[ -f "$target_msa_dir/${target_name}.a3m" ]]; then
                log_info "  ✓ Generated: $target_name"
            else
                log_error "  ✗ Failed: $target_name"
                all_generated=false
            fi
        done
        
        if [[ "$all_generated" != "true" ]]; then
            log_error "Some target MSAs failed to generate"
            log_error "Check log: $target_msa_dir/colabfold_batch.log"
            exit 1
        fi
    else
        log_info "All target MSAs already exist - skipping generation"
    fi
    
    # Export for use in calling scope
    export TARGET_MSA_DIR="$target_msa_dir"
    export TARGET_NAMES="${target_names[*]}"
    export N_TARGETS="${#target_names[@]}"
}

#===============================================================================
# SHARED FUNCTION: Convert MSAs to Chai-1 Parquet Format
#===============================================================================

convert_to_chai_parquet() {
    local source_dir="$1"
    local output_dir="$2"
    local label="$3"
    local log_suffix="$4"
    
    log_info "Converting $label MSAs to Chai-1 parquet format..."
    log_info "Source directory: $source_dir"
    log_info "Output directory: $output_dir"
    
    local conversion_script="$SCRIPT_DIR/a3m_to_pqt_parallel.py"
    
    python "$conversion_script" \
        "$source_dir" \
        "$output_dir" \
        --workers auto \
        --source uniref90 \
        2>&1 | tee "$output_dir/conversion_${log_suffix}.log"
    
    return ${PIPESTATUS[0]}
}

#===============================================================================
# TARGET-ONLY MSA MODE (binder_msa=false)
# Generate MSAs for targets only, then create multimer MSAs using templates
#===============================================================================

if [[ "$binder_msa" == "false" ]]; then
    log_stage "TARGET-ONLY MSA MODE"
    log_info "Generating target MSAs to use as templates for multimer MSAs"
    
    # Generate target MSAs (shared step)
    generate_target_msas
    target_msa_dir="$TARGET_MSA_DIR"
    
    #===========================================================================
    # Generate multimer MSAs using template
    #===========================================================================
    
    log_stage "Generating Multimer MSAs using target MSAs"
    
    log_info "Generating multimer MSAs"
    
    python "$SCRIPT_DIR/generate_multimer_a3m.py" \
        --input_csv "$csv" \
        --target_msa_dir "$target_msa_dir" \
        --output_dir "$msa_dir_multimer" \
        --validate \
        2>&1 | tee "$msa_dir_multimer/generation.log"
    
    # Check if successful
    if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
        n_multimer=$(find "$msa_dir_multimer" -name "*.a3m" -type f | wc -l)
        log_info "Successfully generated $n_multimer multimer MSAs using templates"
    else
        log_error "Multimer MSA generation failed"
        log_error "Check log: $msa_dir_multimer/generation.log"
        exit 1
    fi

    #===========================================================================
    # Convert to Chai-1 parquet (target MSAs only)
    #===========================================================================
    
    log_stage "Converting MSAs to Chai-1 Parquet Format"
    
    chai_msa_dir="$project/chai/msas"
    mkdir -p "$chai_msa_dir"
    
    # Deactivate ColabFold environment
    log_info "Deactivating ColabFold environment..."
    conda deactivate
    
    # Activate Chai environment
    log_info "Activating Chai environment..."
    conda activate chai || {
        log_error "Failed to activate chai environment"
        log_error "Make sure 'conda activate chai' works in your shell"
        exit 1
    }
    
    # Convert target MSAs only (Chai-1 needs individual chain MSAs, not multimer)
    if ! convert_to_chai_parquet "$target_msa_dir" "$chai_msa_dir" "target" "targets"; then
        log_error "Chai-1 MSA conversion failed"
        log_error "Check log: $chai_msa_dir/conversion_targets.log"
        exit 1
    fi

    n_parquet=$(find "$chai_msa_dir" -name "*.aligned.pqt" -type f 2>/dev/null | wc -l)
    if (( n_parquet == 0 )); then
        log_error "No parquet files created - Chai-1 MSA conversion failed"
        exit 1
    fi
    log_info "Successfully converted $n_parquet target MSAs to Chai-1 format"
    touch "$chai_msa_dir/.conversion_complete"
    
    # Return to ColabFold environment for consistency
    log_info "Returning to ColabFold environment..."
    conda activate colabfold
    
    #===========================================================================
    # Validation and completion
    #===========================================================================
    
    log_info "Validating outputs..."
    
    n_expected=$(awk -F',' 'NR>1 && $0 ~ /[^[:space:]]/' "$csv" | wc -l)
    n_multimer=$(find "$msa_dir_multimer" -name "*.a3m" -type f | wc -l)
    n_target_a3m=$(find "$target_msa_dir" -maxdepth 1 -name "*.a3m" -type f | wc -l)
    n_chai=$(find "$chai_msa_dir" -name "*.aligned.pqt" -type f 2>/dev/null | wc -l)
    if (( n_chai == 0 )); then
        n_chai=$(find "$chai_msa_dir" -name "*.parquet" -type f 2>/dev/null | wc -l)
    fi
    
    log_info "Validation results:"
    log_info "  Expected binders: $n_expected"
    log_info "  Target MSAs: $n_target_a3m"
    log_info "  Multimer MSAs: $n_multimer"
    log_info "  Chai-1 parquet (targets): $n_chai"
    
    if (( n_multimer != n_expected )); then
        log_error "Too few multimer MSAs generated: $n_multimer < $n_expected"
        exit 1
    fi
    
    # Mark completion
    cat > "$project/.msa_complete" <<EOF
# MSA Generation Complete (Target-Only Mode)
# Generated: $(date)
# Method: Template-based (fast!)
# Targets: $N_TARGETS ($TARGET_NAMES)
# Binders processed: $n_expected
# Multimer MSAs: $n_multimer
# Chai-1 MSAs (targets): $n_chai
EOF
    
    conda deactivate
    
    log_info "Target-only MSA generation complete!"
    log_info "Completion marker: $project/.msa_complete"
    log_info ""
    log_info "Summary:"
    log_info "  - Target MSAs: $target_msa_dir/*.a3m (reusable)"
    log_info "  - Multimer MSAs: $msa_dir_multimer/*.a3m (template-based, for ColabFold)"
    log_info "  - Chai-1 MSAs: $chai_msa_dir/*.aligned.pqt (target sequences only)"
    log_info ""
    log_info "Next step: predict_structure.sh"
    
    exit 0
fi

#===============================================================================
# BINDER MSA MODE (binder_msa=true)
# 1. Generate target MSAs (shared step)
# 2. Generate binder MSAs (via --msa-method server or mem)
# 3. Combine target + binder MSAs into multimer A3Ms
# 4. Convert target + binder MSAs to Chai-1 parquet
#===============================================================================
if [[ "$binder_msa" == "true" ]]; then

    log_stage "BINDER MSA MODE"
    log_info "Generating MSAs for target:binder complexes"

    # Generate target MSAs (shared step)
    generate_target_msas
    target_msa_dir="$TARGET_MSA_DIR"

    #===============================================================================
    # Generate Binder MSAs
    #===============================================================================

    log_stage "Generating Binder MSAs (method: $msa_method)"

    binder_msa_outdir="$msa_dir_binder"

    # Stale MSA check: if a3m files exist but count doesn't match expected, wipe and regenerate
    n_expected=$(awk -F',' 'NR>1 && $0 ~ /[^[:space:]]/' "$csv" | wc -l)
    existing_binder_a3m=($(find "$binder_msa_outdir" -maxdepth 1 -name "*.a3m" -type f 2>/dev/null))
    n_existing_binder=${#existing_binder_a3m[@]}

    if (( n_existing_binder > 0 )); then
        if (( n_existing_binder == n_expected )); then
            log_info "Binder MSAs already exist ($n_existing_binder files match expected $n_expected) - skipping generation"
        else
            log_warn "Binder MSA count mismatch: found $n_existing_binder, expected $n_expected"
            log_warn "Wiping $binder_msa_outdir and regenerating all binder MSAs"
            rm -f "$binder_msa_outdir"/*.a3m
            n_existing_binder=0
        fi
    fi

    if (( n_existing_binder == 0 )); then
        if [[ "$msa_method" == "server" ]]; then
            #===================================================================
            # Server method: colabfold_batch --msa-only
            #===================================================================
            log_info "Generating binder MSAs using ColabFold server..."
            log_info "Input: $colabfold_input_binder"
            log_info "Output: $binder_msa_outdir"

            colabfold_batch \
                --msa-only \
                "$colabfold_input_binder" \
                "$binder_msa_outdir" \
                >> "$binder_msa_outdir/colabfold_batch.log" 2>&1

            n_binder_a3m=$(find "$binder_msa_outdir" -maxdepth 1 -name "*.a3m" -type f | wc -l)
            log_info "Generated $n_binder_a3m binder MSA files"

            if (( n_binder_a3m == 0 )); then
                log_error "No binder MSA files generated"
                log_error "Check log: $binder_msa_outdir/colabfold_batch.log"
                exit 1
            fi

        elif [[ "$msa_method" == "mem" ]]; then
            #===================================================================
            # Memory method: vmtouch + colabfold_search
            #===================================================================

            # Load database into memory
            load_database_to_memory() {
                local db_path="$1"
                
                if ! command -v vmtouch >/dev/null 2>&1; then
                    log_warn "vmtouch not available - skipping database preload"
                    log_warn "Performance may be degraded on first run"
                    return 0
                fi
                
                local idx_files=("$db_path"/*.idx*)
                if [[ ${#idx_files[@]} -eq 0 ]]; then
                    log_warn "No .idx files found in $db_path"
                    return 0
                fi
                
                log_info "Checking database memory status..."
                local current_resident=$(vmtouch -v "$db_path"/*.idx* 2>/dev/null | \
                    awk -F': *' '/Resident Pages:/ {print $2; exit}' | \
                    awk '{gsub(/[()]/, "", $NF); print $NF}')
                
                if [[ -n "$current_resident" ]]; then
                    local pct="${current_resident%\%}"
                    if (( $(echo "$pct > 80" | bc -l 2>/dev/null || echo 0) )); then
                        log_info "Database already ${current_resident} in memory - skipping load"
                        return 0
                    fi
                fi
                
                log_info "Loading database index files into memory..."
                log_info "This requires ~600GB RAM and may take up to 1 hour"
                log_info "Progress will be shown every 10 seconds"
                
                vmtouch -t "$db_path"/*.idx* >/tmp/vmtouch_load.log 2>&1 &
                local pid=$!
                
                local last_pct="0%"
                while kill -0 "$pid" 2>/dev/null; do
                    local line="$(vmtouch -v "$db_path"/*.idx* 2>/dev/null | \
                        awk -F': *' '/Resident Pages:/ {print $2; exit}')"
                    
                    if [[ -n "$line" ]]; then
                        local pct="$(awk '{print $NF}' <<<"$line")"
                        local pages="$(awk '{print $1}' <<<"$line")"
                        local bytes="$(awk '{print $2}' <<<"$line")"
                        
                        if [[ "$pct" != "$last_pct" ]]; then
                            printf "\rProgress: %6s  (pages %s, bytes %s)" "$pct" "$pages" "$bytes"
                            last_pct="$pct"
                        fi
                    fi
                    
                    sleep 10
                done
                
                wait "$pid"
                echo ""
                log_info "Database index files loaded into memory"
            }

            # Load database
            if [[ -d "$db" ]]; then
                load_database_to_memory "$db"
            else
                log_error "Database path not found or not a directory: $db"
                exit 1
            fi

            log_info "Generating binder MSAs using local database..."
            log_info "Input: $colabfold_input_binder"
            log_info "Output: $binder_msa_outdir"

            CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES" colabfold_search \
                "$colabfold_input_binder" \
                "$db" \
                "$binder_msa_outdir" \
                --gpu 1 \
                --db-load-mode 2 \
                >> "$binder_msa_outdir/colabfold_search.log" 2>&1

            n_binder_a3m=$(find "$binder_msa_outdir" -maxdepth 1 -name "*.a3m" -type f | wc -l)
            log_info "Generated $n_binder_a3m binder MSA files"

            if (( n_binder_a3m == 0 )); then
                log_error "No binder MSA files generated"
                log_error "Check log: $binder_msa_outdir/colabfold_search.log"
                exit 1
            fi
        fi
    fi

    #===============================================================================
    # Combine Target + Binder MSAs into Multimer MSAs
    #===============================================================================

    log_stage "Combining Target + Binder MSAs into Multimer MSAs"

    # Stale multimer MSA check
    existing_multimer_a3m=($(find "$msa_dir_multimer" -name "*.a3m" -type f 2>/dev/null))
    n_existing_multimer=${#existing_multimer_a3m[@]}

    if (( n_existing_multimer > 0 )); then
        if (( n_existing_multimer == n_expected )); then
            log_info "Multimer MSAs already exist ($n_existing_multimer files match expected $n_expected) - skipping"
        else
            log_warn "Multimer MSA count mismatch: found $n_existing_multimer, expected $n_expected"
            log_warn "Wiping $msa_dir_multimer and regenerating all multimer MSAs"
            rm -f "$msa_dir_multimer"/*.a3m
            n_existing_multimer=0
        fi
    fi

    if (( n_existing_multimer == 0 )); then
        log_info "Combining MSAs using combine_target_binder_a3m.py"

        python "$SCRIPT_DIR/combine_target_binder_a3m.py" \
            --input_csv "$csv" \
            --target_msa_dir "$target_msa_dir" \
            --binder_msa_dir "$msa_dir_binder" \
            --output_dir "$msa_dir_multimer" \
            --validate \
            2>&1 | tee "$msa_dir_multimer/generation.log"

        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            log_error "Multimer MSA combination failed"
            log_error "Check log: $msa_dir_multimer/generation.log"
            exit 1
        fi

        n_multimer_a3m=$(find "$msa_dir_multimer" -name "*.a3m" -type f | wc -l)
        log_info "Successfully generated $n_multimer_a3m multimer MSAs"
    fi

    #===============================================================================
    # Convert MSAs to Chai-1 Parquet Format (target + binder MSAs)
    #===============================================================================

    log_stage "Converting MSAs to Chai-1 Parquet Format"

    chai_msa_dir="$project/chai/msas"
    mkdir -p "$chai_msa_dir"

    # Deactivate ColabFold environment
    log_info "Deactivating ColabFold environment..."
    conda deactivate

    # Activate Chai environment
    log_info "Activating Chai environment..."
    conda activate chai || {
        log_error "Failed to activate chai environment"
        log_error "Make sure 'conda activate chai' works in your shell"
        exit 1
    }

    # Convert target MSAs
    log_info "Converting target MSAs..."
    if ! convert_to_chai_parquet "$target_msa_dir" "$chai_msa_dir" "target" "targets"; then
        log_error "Target MSA conversion to Chai-1 parquet failed"
        log_error "Check log: $chai_msa_dir/conversion_targets.log"
        exit 1
    fi
    n_target_parquet=$(find "$chai_msa_dir" -name "*.aligned.pqt" -type f 2>/dev/null | wc -l)
    log_info "Converted $n_target_parquet target MSAs"

    # Convert binder MSAs
    log_info "Converting binder MSAs..."
    if ! convert_to_chai_parquet "$msa_dir_binder" "$chai_msa_dir" "binder" "binders"; then
        log_error "Binder MSA conversion to Chai-1 parquet failed"
        log_error "Check log: $chai_msa_dir/conversion_binders.log"
        exit 1
    fi
    n_total_parquet=$(find "$chai_msa_dir" -name "*.aligned.pqt" -type f 2>/dev/null | wc -l)
    log_info "Total Chai-1 parquet files: $n_total_parquet"

    # Validate parquet count
    if (( n_total_parquet == 0 )); then
        log_error "No parquet files created - Chai-1 MSA conversion failed"
        exit 1
    fi
    touch "$chai_msa_dir/.conversion_complete"

    # Return to ColabFold environment for consistency
    log_info "Returning to ColabFold environment..."
    conda activate colabfold

    #===============================================================================
    # Validation and Completion
    #===============================================================================

    log_info "Validating MSA generation..."

    n_multimer_a3m=$(find "$msa_dir_multimer" -name "*.a3m" -type f | wc -l)
    n_binder_a3m=$(find "$msa_dir_binder" -maxdepth 1 -name "*.a3m" -type f | wc -l)
    n_target_a3m=$(find "$target_msa_dir" -maxdepth 1 -name "*.a3m" -type f | wc -l)
    n_chai_parquet=$(find "$chai_msa_dir" -name "*.aligned.pqt" -type f 2>/dev/null | wc -l)

    log_info "Validation results:"
    log_info "  Expected sequences: $n_expected"
    log_info "  Target MSAs: $n_target_a3m"
    log_info "  Binder MSAs: $n_binder_a3m"
    log_info "  Multimer MSAs: $n_multimer_a3m"
    log_info "  Chai-1 parquet (targets + binders): $n_chai_parquet"

    if (( n_multimer_a3m != n_expected )); then
        log_error "Too few multimer MSAs generated: $n_multimer_a3m < $n_expected"
        exit 1
    fi

    if (( n_binder_a3m != n_expected )); then
        log_error "Too few binder MSAs generated: $n_binder_a3m < $n_expected"
        exit 1
    fi

    #===============================================================================
    # Mark completion
    #===============================================================================

    cat > "$project/.msa_complete" <<EOF
# MSA Generation Complete (Binder Mode)
# Generated: $(date)
# Method: $msa_method
# Targets: $N_TARGETS ($TARGET_NAMES)
# Sequences: $n_expected
# Target MSAs: $n_target_a3m
# Binder MSAs: $n_binder_a3m
# Multimer MSAs: $n_multimer_a3m
# Chai-1 MSAs (targets + binders): $n_chai_parquet
EOF

    conda deactivate

    log_info "MSA generation complete!"
    log_info "Completion marker: $project/.msa_complete"
    log_info ""
    log_info "Summary:"
    log_info "  - Target MSAs: $target_msa_dir/*.a3m (reusable)"
    log_info "  - Binder MSAs: $msa_dir_binder/*.a3m"
    log_info "  - Multimer MSAs: $msa_dir_multimer/*.a3m (combined target + binder, for ColabFold)"
    log_info "  - Chai-1 MSAs: $chai_msa_dir/*.aligned.pqt (targets + binders)"
    log_info ""
    log_info "Next step: predict_structure.sh"

    exit 0
fi
