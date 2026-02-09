#!/usr/bin/env bash
#===============================================================================
# check.sh - Dependency checker for the binderScore pipeline
#
# Checks all required:
#   - System binaries
#   - Conda installation & environments
#   - Per-environment Python packages
#   - Per-environment CLI tools
#   - GPU availability
#   - Required file paths
#   - Pipeline scripts
#
# Usage:
#   bash check.sh [--script-dir /path/to/binderScore]
#
# Exit code: 0 if all critical dependencies met, 1 otherwise
#===============================================================================

set -uo pipefail

#===============================================================================
# Configuration
#===============================================================================

SCRIPT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --script-dir) SCRIPT_DIR="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: $0 [--script-dir /path/to/binderScore]"
            echo ""
            echo "Checks all dependencies for the binderScore pipeline."
            echo "If --script-dir is not provided, checks for scripts in ~/bin/binderScore"
            exit 0
            ;;
        *) echo "Unknown flag: $1"; exit 1 ;;
    esac
done

[[ -z "$SCRIPT_DIR" ]] && SCRIPT_DIR="$HOME/bin/binderScore"

#===============================================================================
# Counters & formatting
#===============================================================================

PASS=0
FAIL=0
WARN=0

MISSING_CRITICAL=()
MISSING_OPTIONAL=()

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

pass() {
    echo -e "  ${GREEN}✓${NC} $*"
    PASS=$((PASS+1))
}

fail() {
    echo -e "  ${RED}✗${NC} $*"
    FAIL=$((FAIL+1))
    MISSING_CRITICAL+=("$1")
}

warn() {
    echo -e "  ${YELLOW}⚠${NC} $*"
    WARN=$((WARN+1))
    MISSING_OPTIONAL+=("$1")
}

section() {
    echo ""
    echo -e "${BOLD}${CYAN}━━━ $* ━━━${NC}"
}

#===============================================================================
# 1. System binaries
#===============================================================================

section "System Binaries"

check_binary() {
    local name="$1"
    local critical="${2:-true}"

    if command -v "$name" &>/dev/null; then
        local location
        location=$(command -v "$name")
        pass "$name  ($location)"
    else
        if [[ "$critical" == "true" ]]; then
            fail "$name" "MISSING (critical)"
        else
            warn "$name" "MISSING (optional, but recommended)"
        fi
    fi
}

check_binary bash
check_binary awk
check_binary grep
check_binary find
check_binary nproc
check_binary bc
check_binary jq
check_binary nvidia-smi
check_binary vmtouch false

#===============================================================================
# 2. GPU availability
#===============================================================================

section "GPU Availability"

if command -v nvidia-smi &>/dev/null; then
    NGPUS=$(nvidia-smi -L 2>/dev/null | wc -l)
    if (( NGPUS > 0 )); then
        pass "NVIDIA GPUs detected: $NGPUS"
        nvidia-smi -L 2>/dev/null | while read -r line; do
            echo "       $line"
        done
    else
        fail "NVIDIA GPUs" "nvidia-smi found but no GPUs detected"
    fi
else
    fail "NVIDIA GPUs" "nvidia-smi not found — no GPU support"
fi

#===============================================================================
# 3. Conda installation
#===============================================================================

section "Conda Installation"

CONDA_SH="$HOME/miniforge3/etc/profile.d/conda.sh"
if [[ -f "$CONDA_SH" ]]; then
    pass "conda.sh found ($CONDA_SH)"
    # shellcheck disable=SC1090
    source "$CONDA_SH"
else
    fail "conda.sh" "Not found at $CONDA_SH"
    echo -e "  ${RED}   Cannot proceed with environment checks without conda${NC}"
    # Print summary and exit early
    section "Summary"
    echo -e "  ${GREEN}Passed:${NC}  $PASS"
    echo -e "  ${RED}Failed:${NC}  $FAIL"
    echo -e "  ${YELLOW}Warned:${NC} $WARN"
    exit 1
fi

#===============================================================================
# 4. Conda environments
#===============================================================================

section "Conda Environments"

COLABFOLD_ENV="colabfold"
CHAI_ENV="chai"
BOLTZ_ENV="boltz"
PYROSETTA_ENV="pyrosetta"

check_conda_env() {
    local env_name="$1"
    local env_path="$2"

    if conda activate "$env_path" 2>/dev/null; then
        pass "Environment: $env_name  ($env_path)"
        conda deactivate 2>/dev/null
    else
        fail "Environment: $env_name" "Cannot activate: $env_path"
    fi
}

check_conda_env "colabfold" "$COLABFOLD_ENV"
check_conda_env "chai"      "$CHAI_ENV"
check_conda_env "boltz"     "$BOLTZ_ENV"
check_conda_env "pyrosetta" "$PYROSETTA_ENV"

#===============================================================================
# 5. ColabFold environment — CLI tools
#===============================================================================

section "ColabFold Environment (CLI tools)"

if conda activate "$COLABFOLD_ENV" 2>/dev/null; then
    for tool in mmseqs colabfold_search colabfold_batch; do
        if command -v "$tool" &>/dev/null; then
            pass "$tool"
        else
            fail "$tool" "Not found in colabfold environment"
        fi
    done
    conda deactivate 2>/dev/null
else
    fail "colabfold tools" "Cannot activate colabfold environment — skipping"
fi

#===============================================================================
# 6. Chai environment — Python packages
#===============================================================================

section "Chai Environment (Python packages)"

if conda activate "$CHAI_ENV" 2>/dev/null; then
    for pkg in chai_lab pandas numpy biopython; do
        # Map display names to import names
        import_name="$pkg"
        [[ "$pkg" == "biopython" ]] && import_name="Bio"
        [[ "$pkg" == "chai_lab" ]] && import_name="chai_lab"

        if python -c "import $import_name" 2>/dev/null; then
            version=$(python -c "import $import_name; print(getattr($import_name, '__version__', 'unknown'))" 2>/dev/null)
            pass "$pkg  (v$version)"
        else
            fail "$pkg" "Not importable in chai environment"
        fi
    done

    # Check specific chai_lab submodules used in scripts
    for submod in "chai_lab.chai1" "chai_lab.data.parsing.msas.aligned_pqt" "chai_lab.data.parsing.msas.data_source"; do
        if python -c "import $submod" 2>/dev/null; then
            pass "$submod"
        else
            fail "$submod" "Submodule not importable"
        fi
    done

    conda deactivate 2>/dev/null
else
    fail "chai packages" "Cannot activate chai environment — skipping"
fi

#===============================================================================
# 7. Boltz environment — Python packages & CLI
#===============================================================================

section "Boltz Environment (Python packages & CLI)"

if conda activate "$BOLTZ_ENV" 2>/dev/null; then
    # Check boltz CLI
    if command -v boltz &>/dev/null; then
        pass "boltz CLI"
    else
        fail "boltz CLI" "Not found in boltz environment"
    fi

    # Check Python packages
    for pkg in boltz numpy; do
        if python -c "import $pkg" 2>/dev/null; then
            version=$(python -c "import $pkg; print(getattr($pkg, '__version__', 'unknown'))" 2>/dev/null)
            pass "$pkg  (v$version)"
        else
            fail "$pkg" "Not importable in boltz environment"
        fi
    done

    conda deactivate 2>/dev/null
else
    fail "boltz packages" "Cannot activate boltz environment — skipping"
fi

#===============================================================================
# 8. PyRosetta environment — Python packages & binaries
#===============================================================================

section "PyRosetta Environment (Python packages & binaries)"

if conda activate "$PYROSETTA_ENV" 2>/dev/null; then

    # Python packages
    for pkg in pyrosetta numpy scipy pandas biopython; do
        import_name="$pkg"
        [[ "$pkg" == "biopython" ]] && import_name="Bio"

        if python -c "import $import_name" 2>/dev/null; then
            version=$(python -c "import $import_name; print(getattr($import_name, '__version__', 'unknown'))" 2>/dev/null)
            pass "$pkg  (v$version)"
        else
            fail "$pkg" "Not importable in pyrosetta environment"
        fi
    done

    # Check pyrosetta submodules
    for submod in \
        "pyrosetta.rosetta.protocols.relax" \
        "pyrosetta.rosetta.protocols.analysis" \
        "pyrosetta.rosetta.core.select.residue_selector" \
        "pyrosetta.rosetta.core.simple_metrics.metrics"; do
        if python -c "import $submod" 2>/dev/null; then
            pass "$submod"
        else
            fail "$submod" "PyRosetta submodule not importable"
        fi
    done

    # Check SAP score (optional)
    if python -c "from pyrosetta.rosetta.core.pack.guidance_scoreterms.sap import SapScoreMetric" 2>/dev/null; then
        pass "SapScoreMetric (SAP scoring)"
    else
        warn "SapScoreMetric" "Not available — SAP scores will be NaN"
    fi

    # Check mkdssp / dssp binary (used by DSSP in biopython)
    if command -v mkdssp &>/dev/null; then
        pass "mkdssp"
    elif command -v dssp &>/dev/null; then
        pass "dssp (as mkdssp alternative)"
    else
        fail "mkdssp/dssp" "Not found — secondary structure calculation will fail"
    fi

    # Check DAlphaBall binary
    DALPHABALL_PATH="/home/ubuntu/miniforge3/envs/pyrosetta/bin/dalphaball"
    if [[ -f "$DALPHABALL_PATH" ]] && [[ -x "$DALPHABALL_PATH" ]]; then
        pass "DAlphaBall  ($DALPHABALL_PATH)"
    else
        warn "DAlphaBall" "Not found at $DALPHABALL_PATH — hole scoring disabled"
    fi

    conda deactivate 2>/dev/null
else
    fail "pyrosetta packages" "Cannot activate pyrosetta environment — skipping"
fi

#===============================================================================
# 9. Pipeline scripts
#===============================================================================

section "Pipeline Scripts (in $SCRIPT_DIR)"

if [[ -d "$SCRIPT_DIR" ]]; then
    pass "Script directory exists"

    EXPECTED_SCRIPTS=(
        "binderScore.sh"
        "generate_msas/generate_msas.sh"
        "generate_msas/generate_multimer_a3m.py"
        "generate_msas/a3m_to_pqt_parallel.py"
        "predict_structure/predict_structure.sh"
        "predict_structure/generate_boltz_yaml.py"
        "predict_structure/calculate_distances.py"
        "predict_structure/filter_and_select.py"
        "predict_structure/prepare_retry_inputs.py"
        "predict_structure/chai_run.py"
        "predict_structure/rerun_filter.sh"
        "score/score.sh"
        "score/colabfold_extract_metrics.py"
        "score/boltz_extract_metrics.py"
        "score/chai_extract_metrics.py"
        "score/binder_extract_metrics.py"
        "score/compute_rosetta_metrics.py"
        "score/merge.py"
    )

    for script in "${EXPECTED_SCRIPTS[@]}"; do
        full_path="$SCRIPT_DIR/$script"
        if [[ -f "$full_path" ]]; then
            if [[ -x "$full_path" ]] || [[ "$script" == *.py ]]; then
                pass "$script"
            else
                warn "$script" "Exists but not executable"
            fi
        else
            fail "$script" "MISSING"
        fi
    done

else
    fail "Script directory" "Not found: $SCRIPT_DIR"
fi

#===============================================================================
# 10. Key filesystem paths
#===============================================================================

section "Key Filesystem Paths"

check_path() {
    local label="$1"
    local path="$2"
    local critical="${3:-true}"

    if [[ -e "$path" ]]; then
        pass "$label  ($path)"
    else
        if [[ "$critical" == "true" ]]; then
            fail "$label" "Not found: $path"
        else
            warn "$label" "Not found: $path"
        fi
    fi
}

check_path "Miniforge3"        "$HOME/miniforge3"
check_path "ColabFold install" "$HOME/localcolabfold" false
check_path "~/bin directory"   "$HOME/bin" false

#===============================================================================
# Summary
#===============================================================================

section "Summary"

TOTAL=$((PASS + FAIL + WARN))

echo -e "  ${GREEN}Passed:${NC}   $PASS / $TOTAL"
echo -e "  ${RED}Failed:${NC}   $FAIL / $TOTAL"
echo -e "  ${YELLOW}Warnings:${NC} $WARN / $TOTAL"
echo ""

if (( FAIL > 0 )); then
    echo -e "${RED}${BOLD}Missing critical dependencies:${NC}"
    for item in "${MISSING_CRITICAL[@]}"; do
        echo -e "  ${RED}•${NC} $item"
    done
    echo ""
fi

if (( WARN > 0 )); then
    echo -e "${YELLOW}${BOLD}Missing optional dependencies:${NC}"
    for item in "${MISSING_OPTIONAL[@]}"; do
        echo -e "  ${YELLOW}•${NC} $item"
    done
    echo ""
fi

if (( FAIL == 0 )); then
    echo -e "${GREEN}${BOLD}All critical dependencies satisfied! Pipeline is ready to run.${NC}"
    exit 0
else
    echo -e "${RED}${BOLD}$FAIL critical dependencies missing. Please install them before running the pipeline.${NC}"
    echo ""
    echo -e "${BOLD}Quick install reference:${NC}"
    echo "  jq:          sudo apt install -y jq"
    echo "  vmtouch:     git clone https://github.com/hoytech/vmtouch.git && cd vmtouch && make && sudo make install"
    echo "  bc:          sudo apt install -y bc"
    echo "  Conda envs:  See setup.sh for full installation instructions"
    exit 1
fi