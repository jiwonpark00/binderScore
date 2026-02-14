# BinderScore

This is a practice project – most definitely not bug-free and coded with lots of help from Claude. 

A comprehensive, multi-method computational pipeline for screening *de novo* protein binder designs. BinderScore orchestrates MSA generation, structure prediction across four complementary methods (AlphaFold2-Multimer, AlphaFold2-pTM, Chai-1, and Boltz-2), distance-based filtering, confidence metric extraction, Rosetta interface scoring, and cross-method RMSD calculation — all in a single end-to-end workflow.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Setup](#setup)
- [Input Preparation](#input-preparation)
- [Usage](#usage)
- [Arguments and Options](#arguments-and-options)
- [Pipeline Stages in Detail](#pipeline-stages-in-detail)
- [Output Files](#output-files)
- [Metrics Reference](#metrics-reference)
- [Utility Scripts](#utility-scripts)
- [Caveats and Important Notes](#caveats-and-important-notes)
- [Acknowledgements](#acknowledgements)

---

## Overview

BinderScore takes a CSV of binder designs (each paired with a target protein) and runs them through a three-stage pipeline:

1. **Stage 1 — MSA Generation**: Generates multiple sequence alignments for targets and (optionally) binders, then assembles multimer MSAs for ColabFold and parquet MSAs for Chai-1.
2. **Stage 2 — Structure Prediction**: Runs four structure prediction methods in parallel (AF2-Multimer, AF2-pTM, Chai-1, Boltz-2) plus a binder-only prediction. Calculates distances from predicted structures to specified target residues, filters by a distance cutoff, and automatically retries failed predictions with a different random seed.
3. **Stage 3 — Scoring**: Extracts confidence metrics from each method's output, computes Rosetta interface energetics (with FastRelax), calculates cross-method binder RMSD, and merges everything into a single wide-format CSV.

The final output is a per-binder score table consolidating metrics from all methods, enabling informed ranking and downstream selection.

---

## Pipeline Architecture

```
binderScore.sh                      ← Top-level orchestrator
│
├── Stage 1: generate_msas/
│   ├── generate_msas.sh            ← MSA generation controller
│   ├── generate_multimer_a3m.py    ← Template-based multimer MSA (target-only mode)
│   ├── combine_target_binder_a3m.py← Combine target+binder MSAs (binder-MSA mode)
│   └── a3m_to_pqt_parallel.py     ← Convert A3M → Chai-1 parquet format
│
├── Stage 2: predict_structure/
│   ├── predict_structure.sh        ← Structure prediction controller
│   ├── generate_boltz_yaml.py      ← Create Boltz-2 YAML inputs
│   ├── chai_run.py                 ← Chai-1 inference wrapper
│   ├── calculate_distances.py      ← Measure target-residue to binder distance
│   ├── filter_and_select.py        ← Distance-based filtering & best-model selection
│   ├── prepare_retry_inputs.py     ← Set up retry inputs for failed predictions
│   └── rerun_filter.sh             ← Standalone re-filter with a new cutoff
│
├── Stage 3: score/
│   ├── score.sh                    ← Scoring controller
│   ├── colabfold_extract_metrics.py← AF2-multimer & AF2-pTM metric extraction
│   ├── boltz_extract_metrics.py    ← Boltz-2 metric extraction (includes PDE)
│   ├── chai_extract_metrics.py     ← Chai-1 metric extraction
│   ├── binder_extract_metrics.py   ← Binder-only metric extraction
│   ├── compute_rosetta_metrics.py  ← PyRosetta interface scoring + FastRelax
│   ├── compute_rmsd.py             ← Cross-method pairwise binder RMSD
│   └── merge.py                    ← Merge all metrics into final CSV
│
├── check.sh                        ← Dependency checker
└── setup.sh                        ← Installation reference (not executable)
```

---

## Requirements

### Hardware

- **GPU**: One or more NVIDIA GPUs (recommend at least 40GB VRAM for efficient pipeline)  (required for all structure prediction methods). The pipeline auto-detects GPU count and parallelises ColabFold and Chai-1 jobs across them.
- **RAM**: At least 64 GB for standard use. If using `--msa-method mem` (in-memory ColabFold database search), **≥ 800 GB RAM** is required. The pipeline will automatically fall back to `--msa-method server` if insufficient RAM is detected.
- **Disk**: Structure predictions and relaxed PDBs can be large. Plan for at least 5–10 GB per 100 binders.

### Software

- **OS**: Ubuntu 24 (tested); other Linux distributions should work.
- **Conda**: Miniforge3 (installed via `setup.sh`).
- **CUDA**: Compatible NVIDIA drivers and CUDA toolkit (12.4 recommended).
- **System utilities**: `bash`, `awk`, `grep`, `find`, `nproc`, `bc`, `jq`, `nvidia-smi`.
- **Optional**: `vmtouch` (for preloading ColabFold database into memory).

### Conda Environments

The pipeline uses four separate conda environments to isolate incompatible dependencies:

| Environment | Purpose | Key Packages |
|---|---|---|
| `colabfold` | MSA generation and AF2 predictions | `localcolabfold`, `mmseqs`, `colabfold_batch`, `colabfold_search`, CUDA 12.4 |
| `chai` | Chai-1 structure prediction | `chai_lab==0.6.1`, `pandas`, `numpy`, `biopython` |
| `boltz` | Boltz-2 structure prediction | `boltz[cuda]` |
| `pyrosetta` | Rosetta scoring and relaxation | `pyrosetta`, `pyrosetta-distributed`, `biopython`, `dssp`, `DAlphaBall` |

---

## Setup

> **`setup.sh` is NOT an executable script.** It is a reference document. Each block should be run independently, step by step, to catch and address errors as they arise.

### Step-by-step installation

**1. Check your  (if applicable):**


**2. Create the scripts directory and transfer BinderScore:**
```bash
mkdir -p ~/bin
# From your local machine:
rsync -avzP binderScore ubuntu@<YOUR_HOST>:~/bin/

# From this repo:
git clone https://github.com/jiwonpark00/binderScore ~/bin/binderScore
```

**3. Make scripts executable and add to PATH:**
```bash
ln -sf ~/bin/binderScore/binderScore.sh ~/bin/binderScore.sh
chmod +x ~/bin/*
chmod +x ~/bin/*/*
chmod +x ~/bin/*/*/*
echo 'export PATH="~/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**4. Install Miniforge3 (conda):**
```bash
wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-$(uname)-$(uname -m).sh
$HOME/miniforge3/bin/conda init bash
source ~/.bashrc
conda config --set auto_activate false
source ~/.bashrc
```

**5. Install LocalColabFold:**
```bash
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc
git clone https://github.com/yoshitakamo/localcolabfold.git
cd localcolabfold
pixi install && pixi run setup

conda create -n colabfold -c nvidia cuda-nvcc=12.4 -y
conda activate colabfold
mkdir -p ~/miniforge3/envs/colabfold/etc/conda/activate.d
mkdir -p ~/miniforge3/envs/colabfold/etc/conda/deactivate.d

echo 'export OLD_PATH="$PATH"' > ~/miniforge3/envs/colabfold/etc/conda/activate.d/env_vars.sh
echo 'export PATH="/home/ubuntu/localcolabfold/.pixi/envs/default/bin:$PATH"' >> ~/miniforge3/envs/colabfold/etc/conda/activate.d/env_vars.sh

echo 'export PATH="$OLD_PATH"' > ~/miniforge3/envs/colabfold/etc/conda/deactivate.d/env_vars.sh
echo 'unset OLD_PATH' >> ~/miniforge3/envs/colabfold/etc/conda/deactivate.d/env_vars.sh

conda deactivate
conda activate colabfold

# Verify
which colabfold_batch

conda deactivate
cd ~

rm -f ~/*.sh
```

**6. Install vmtouch (optional but recommended):**
```bash
git clone https://github.com/hoytech/vmtouch.git
cd vmtouch && make
sudo make install
cd ~
```

**7. Install jq:**
```bash
sudo apt-get update
sudo apt install -y jq
```

**8. Install Chai-1:**
```bash
conda create -n chai python==3.10 -y
conda activate chai
pip install --upgrade pip
pip install chai_lab==0.6.1
conda deactivate
```

**9. Install Boltz-2:**
```bash
conda create -n boltz python=3.11 -y
conda activate boltz
pip install boltz[cuda] -U
conda deactivate
```

**10. Install PyRosetta + DAlphaBall:**
```bash
conda create -n pyrosetta python=3.11 -y
conda activate pyrosetta
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(distributed=True,serialization=True)'
pip install pyrosetta-distributed
conda install dssp -y
pip install biopython

git clone https://github.com/outpace-bio/DAlphaBall.git
mkdir DAlphaBall_build
./DAlphaBall/dalphaball_docker_build.sh ~/DAlphaBall_build
conda install dalphaball -c ~/DAlphaBall_build -y
conda deactivate
```

### Verifying the Installation

After completing all steps, run the dependency checker:

```bash
bash ~/bin/binderScore/check.sh
```

This validates all system binaries, conda environments, Python packages, GPU availability, and pipeline scripts. It will report any critical or optional missing dependencies.

---

## Input Preparation

### 1. Binder CSV (`--csv`)

A CSV file with three columns (header required):

```
name,target,binder
AQP_01,MTARGETSEQUENCE...,MBINDERSEQUENCE...
AQP_02,MTARGETSEQUENCE...,MBINDERSEQUENCE...
SLP_01,MDIFFERENTTARGET...,MANOTHERBINDER...
```

**Naming convention is critical.** Binder names **must** follow the format `TARGET_NUMBER` (e.g., `AQP_01`, `IL6R_123`, `my_target_5`). The target name is extracted as everything before the last underscore; the suffix must be purely numeric. This convention is used to map binders to their targets throughout the pipeline.

### 2. Target Directory (`--target-dir`)

A directory containing:
```
targets/
├── targets.fasta     ← FASTA file with all target sequences
├── pdb/              ← One PDB per target (for Boltz-2 templates)
│   ├── AQP.pdb
│   └── SLP.pdb
└── msa/              ← Created automatically by Stage 1
    ├── AQP.a3m
    └── SLP.a3m
```

The `targets.fasta` file:
```
>AQP
MTARGETSEQUENCE...
>SLP
MDIFFERENTTARGET...
```

Target names in `targets.fasta` headers must match the target prefix extracted from binder names in the CSV.

### 3. Residue CSV (`--residue-csv`)

Specifies the target residue of interest for distance-based filtering:

```
target_name,residue_number
AQP,122
SLP,85
```

The pipeline calculates the minimum heavy-atom distance from this residue (in the target chain) to any residue in the binder chain. Structures where this distance exceeds `--distance-cutoff` are marked as "failed".

### 4. ColabFold Database (`--db`, only for `--msa-method mem`)

The ColabFold local database directory (typically 2+ TB). Only required when using in-memory MSA search. The pipeline loads `.idx` files into RAM using `vmtouch` for fast search.

---

## Usage

### Minimal example (target-only MSAs, all defaults)

```bash
binderScore.sh \
  --csv input.csv \
  --project my_project/ \
  --target-dir targets/ \
  --binder-msa false \
  --residue-csv residues.csv \
  --distance-cutoff 20
```

### With binder MSAs via ColabFold server

```bash
binderScore.sh \
  --csv input.csv \
  --project my_project/ \
  --target-dir targets/ \
  --binder-msa true \
  --msa-method server \
  --residue-csv residues.csv \
  --distance-cutoff 20
```

### With binder MSAs via local database + custom parameters

```bash
binderScore.sh \
  --csv input.csv \
  --project my_project/ \
  --target-dir targets/ \
  --binder-msa true \
  --msa-method mem \
  --db /data/colabfold_db \
  --residue-csv residues.csv \
  --distance-cutoff 20 \
  --num-recycle 3 \
  --boltz-diffusion-samples 10 \
  --pae-threshold 8.0
```

---

## Arguments and Options

### Required Flags

| Flag | Description |
|---|---|
| `--csv` | CSV file with `name,target,binder` columns |
| `--project` | Output project directory (created if it doesn't exist) |
| `--target-dir` | Target directory containing `targets.fasta` and `pdb/` subdirectory |
| `--binder-msa` | Whether to generate binder MSAs: `true` or `false` |
| `--residue-csv` | CSV with `target_name,residue_number` columns |
| `--distance-cutoff` | Distance cutoff in Ångströms for filtering |

### Conditional Flags (only when `--binder-msa true`)

| Flag | Description |
|---|---|
| `--msa-method` | `server` (ColabFold MSA server) or `mem` (local DB in memory) |
| `--db` | Path to ColabFold database directory (required if `--msa-method mem`) |

### ColabFold / AlphaFold2 Options (Stage 2)

| Flag | Default | Description |
|---|---|---|
| `--num-recycle` | 2 | Number of recycling iterations |
| `--recycle-early-stop-tol` | 0.5 | Early stop tolerance for recycling |
| `--num-seeds` | 1 | Number of random seeds per model |
| `--num-models` | 5 | Number of AF2 models to run (1–5) |

### Chai-1 Options (Stage 2)

| Flag | Default | Description |
|---|---|---|
| `--chai-num-trunk-recycles` | 3 | Number of trunk recycles |
| `--chai-num-diffn-timesteps` | 200 | Number of diffusion timesteps |
| `--chai-diffusion-samples` | 5 | Number of diffusion samples generated |

### Boltz-2 Options (Stage 2)

| Flag | Default | Description |
|---|---|---|
| `--boltz-recycling-steps` | 3 | Number of recycling steps |
| `--boltz-sampling-steps` | 200 | Number of sampling steps |
| `--boltz-diffusion-samples` | 5 | Number of diffusion samples generated |
| `--use-potentials` | true | Enable inference-time potentials |

### Scoring Thresholds (Stage 3)

| Flag | Default | Description |
|---|---|---|
| `--pae-threshold` | 10.0 | PAE threshold for interface contact quality (AF2 methods) |
| `--pde-threshold` | 2.0 | PDE threshold for Boltz-2 |
| `--interface-threshold` | 4.0 | Heavy-atom distance (Å) to define interface contacts |
| `--pdockq-threshold` | 8.0 | pDockQ distance threshold |

---

## Pipeline Stages in Detail

### Stage 1: MSA Generation (`generate_msas.sh`)

Two modes of operation:

**Target-only mode (`--binder-msa false`):**
1. Generates target MSAs via `colabfold_batch --msa-only` (using the ColabFold server).
2. Creates multimer A3M files by appending binder sequences with gap padding to the target MSA (template-based, via `generate_multimer_a3m.py`). This is fast since no binder MSA search is performed.
3. Converts target MSAs to Chai-1 parquet format (`a3m_to_pqt_parallel.py`).

**Binder MSA mode (`--binder-msa true`):**
1. Generates target MSAs (same as above).
2. Generates binder MSAs via either `colabfold_batch --msa-only` (server) or `colabfold_search` with the local database loaded into memory (mem).
3. Combines target and binder MSAs into multimer A3M files (`combine_target_binder_a3m.py`).
4. Converts both target and binder MSAs to Chai-1 parquet format.

**Completion marker:** `$project/.msa_complete`

### Stage 2: Structure Prediction (`predict_structure.sh`)

Runs five prediction jobs, distributing work across all available GPUs:

1. **AF2-Multimer** — ColabFold with `alphafold2_multimer_v3`, ranked by iPTM.
2. **AF2-pTM Complex** — ColabFold with `alphafold2_ptm` on the concatenated complex, ranked by iPTM.
3. **Chai-1** — Diffusion-based prediction with optional MSAs, generates CIF outputs.
4. **Boltz-2** — Diffusion-based prediction with per-chain MSAs and target PDB templates, generates PDB outputs.
5. **Binder-only** — Either AF2-pTM (if binder MSAs exist) or Chai-1 without MSAs (if not), for single-chain binder structure assessment.

After predictions complete, the pipeline:
- Calculates minimum heavy-atom distances from the specified target residue to the binder chain for every predicted structure.
- Filters structures by the distance cutoff, selecting the best-ranked passing model per (binder, method) pair. If no models pass, the best-ranked model is selected as a fallback.
- **Automatic retry**: Any (binder, method) pair that failed the distance filter is re-predicted with a different random seed, and the results are merged (attempt 1 preferred if it passed, attempt 2 used otherwise).

Individual method failures do not kill the pipeline; errors are tracked and reported at the end.

**Completion marker:** `$project/.prediction_complete`

### Stage 3: Scoring (`score.sh`)

1. **Confidence metric extraction** — Method-specific scripts parse prediction outputs:
   - AF2 methods: pTM, iPTM, pLDDT, PAE-based interface metrics, pDockQ, interface contacts, secondary structure composition.
   - Boltz-2: Similar, plus PDE-based metrics.
   - Chai-1: pLDDT, pTM, aggregate score, interface contacts, pDockQ (no PAE matrix available).
   - Binder-only: pTM, pLDDT, secondary structure composition.

2. **Merge** — All per-method CSVs are merged into a single `binderScore.csv` with prefixed columns (`mul_` for AF2-multimer, `mon_` for AF2-pTM, `bol_` for Boltz, `cha_` for Chai, `bin_` for binder-only).

3. **Symlink creation** — Convenient symlinks to selected structure files are created in `$project/symlinks/`.

4. **Rosetta interface scoring** — For each method, structures are relaxed with PyRosetta `FastRelax` (constrained to start coordinates), then scored using `InterfaceAnalyzerMover`. Metrics include binding energy (dG), interface SASA, shape complementarity, hydrogen bonds, buried unsatisfied H-bonds, packstat, SAP score, and more.

5. **Cross-method RMSD** — For each binder, relaxed structures are aligned on the target chain (chain A) CA atoms, and pairwise binder (chain B) CA RMSD is computed across all method pairs.

6. **Final merge** — `binderScore.csv` + Rosetta metrics + RMSD → `binderFinalScore.csv`.

**Completion marker:** `$project/.scoring_complete`

---

## Output Files

After a successful run, the project directory contains:

| Path | Description |
|---|---|
| `binderFinalScore.csv` | **Primary output.** One row per binder with all metrics from all methods, Rosetta scores, and cross-method RMSD. |
| `binderScore.csv` | Confidence metrics only (before Rosetta/RMSD merge). |
| `binderRosettaScore.csv` | Rosetta interface metrics merged across methods. |
| `RMSD.csv` | Pairwise cross-method binder RMSD. |
| `summary.csv` | Per-(binder, method) selection summary with pass/fail status and distance. |
| `distances.csv` | Raw distances for all predicted structures. |
| `filtered_structures/` | Best-ranked structures per method, organised by method then binder name. |
| `symlinks/` | Flat symlinks to filtered structures (`symlinks/{method}/{binder}.pdb`). |
| `relaxed/` | FastRelax-ed PDB structures (always PDB, even if input was CIF). |
| `metrics/` | Per-method metric CSVs and Rosetta score CSVs. |
| `predictions/` | Raw prediction outputs from all methods. |
| `predictions_retry/` | Retry prediction outputs (if any failures occurred). |
| `.msa_complete` | Stage 1 completion marker. |
| `.prediction_complete` | Stage 2 completion marker. |
| `.scoring_complete` | Stage 3 completion marker. |

---

## Metrics Reference

### Column Prefixes

All metrics columns are prefixed by method:

| Prefix | Method |
|---|---|
| `mul_` | AF2-Multimer |
| `mon_` | AF2-pTM (monomer model on complex) |
| `bol_` | Boltz-2 |
| `cha_` | Chai-1 |
| `bin_` | Binder-only |

### Confidence Metrics (per method)

Exact columns vary by method; common ones include: `ptm`, `iptm`, `plddt_mean`, `plddt_binder`, `plddt_target`, `pae_interface_mean`, `pdockq`, `interface_ncontacts`, `n_interface_residues_binder`, `n_interface_residues_target`, `helix_pct`, `sheet_pct`, `loop_pct`.

Boltz-2 additionally reports `pde_interface_mean` and related PDE metrics.

### Rosetta Interface Metrics (per method)

Each prefixed with the method prefix followed by `rosetta_*`:

| Metric | Description |
|---|---|
| `rosetta_interface_dG` | Binding free energy (REU) |
| `rosetta_interface_dSASA` | Change in solvent-accessible surface area upon binding |
| `rosetta_interface_dG_dSASA_ratio` | Efficiency: dG normalised by dSASA (×100) |
| `rosetta_interface_sc` | Shape complementarity score |
| `rosetta_interface_packstat` | Interface packing quality |
| `rosetta_interface_interface_hbonds` | Number of interface hydrogen bonds |
| `rosetta_interface_hbond_percentage` | H-bonds per interface residue (%) |
| `rosetta_interface_unsat_hbonds` | Buried unsatisfied hydrogen bond donors/acceptors |
| `rosetta_interface_delta_unsat_hbonds_percentage` | Unsatisfied H-bonds per interface residue (%) |
| `rosetta_interface_nres_binder` | Number of binder interface residues |
| `rosetta_interface_nres_target` | Number of target interface residues |
| `rosetta_interface_percent_of_binder_in_intf` | Fraction of binder surface in the interface (%) |
| `rosetta_interface_binder_hydrophobicity_percent` | Hydrophobic fraction of binder interface residues (%) |
| `rosetta_binder_score` | Total Rosetta energy of the binder chain |
| `rosetta_surface_hydrophobicity` | Fraction of surface-exposed apolar residues on binder |
| `sap_binder` | Spatial Aggregation Propensity of binder |
| `sap_target` | Spatial Aggregation Propensity of target |
| `sap_delta` | Δ SAP (binder SAP in complex minus binder SAP alone) |

### Cross-Method RMSD

Columns follow the pattern `{prefix1}_{prefix2}_binder_rmsd`, for example `mul_bol_binder_rmsd` is the binder CA RMSD between AF2-Multimer and Boltz-2 after alignment on the target chain.

---

## Utility Scripts

### `check.sh` — Dependency Checker

```bash
bash check.sh [--script-dir /path/to/binderScore]
```

Validates all system binaries, conda environments, Python packages, GPU availability, file paths, and pipeline scripts. Run this after installation to confirm readiness.

### `rerun_filter.sh` — Re-filter with a Different Distance Cutoff

```bash
bash predict_structure/rerun_filter.sh <project_dir> <new_cutoff>
```

Re-runs the distance filtering step without re-doing structure prediction. Useful when too many or too few structures passed the initial cutoff. Requires `distances.csv` from a previous run.

### Completion Markers

Each stage writes a hidden marker file on success. To **re-run a stage**, delete the corresponding marker:

| Marker | Stage |
|---|---|
| `$project/.msa_complete` | Stage 1: MSA Generation |
| `$project/.prediction_complete` | Stage 2: Structure Prediction |
| `$project/.scoring_complete` | Stage 3: Scoring |

---

## Caveats and Important Notes

### Binder Naming Convention
Binder names **must** follow the `TARGET_NUMBER` format (e.g., `AQP_01`). The target name is everything before the last underscore, and the suffix must be purely numeric. The pipeline uses this convention to map binders to their target sequences, MSAs, PDB templates, and residue definitions. Incorrectly named binders will cause lookup failures.

### `setup.sh` Is Not Executable
Do **not** run `setup.sh` as a single script. Each installation block should be executed independently and verified before proceeding. Errors in one step (e.g., conda environment creation) can cascade silently.

### Conda Environment Isolation
The pipeline switches between four conda environments during execution. It sources `~/miniforge3/etc/profile.d/conda.sh` at runtime. If your conda installation is elsewhere, you will need to update the path in the shell scripts.

### RAM for In-Memory MSA Search
Using `--msa-method mem` requires ≥ 800 GB of RAM to load ColabFold database index files. The pipeline checks available memory and automatically falls back to `--msa-method server` if insufficient. Loading the database can take up to one hour.

### GPU Parallelism
ColabFold predictions are distributed across GPUs by splitting A3M files into per-GPU input directories. Chai-1 jobs are queued across GPUs with throttling to prevent over-subscription. Boltz-2 uses its native multi-GPU support via the `--devices` flag.

### Fail-Loud-But-Continue
In Stage 2, `set -e` is disabled before running prediction methods. If one method fails (e.g., Chai-1 out of memory), the pipeline continues with the remaining methods and reports errors at the end. Check per-GPU and per-method logs in the `predictions/` directory if a method fails.

### Automatic Retry
Predictions that fail the distance filter on the first attempt are automatically retried with a different random seed. The retry only re-runs the specific (binder, method) pairs that failed, not the entire set. Attempt 1 results are preferred when they pass.

### Boltz YAML Directory Naming
Do **not** rename the `boltz/yaml/` input directory. Boltz derives its output subdirectory name from the input directory name. Changing it will break the distance calculation step's file discovery logic.

### DAlphaBall Path
The Rosetta scoring step expects DAlphaBall at `/home/ubuntu/miniforge3/envs/pyrosetta/bin/dalphaball`. If your installation is in a different location, update the `DALPHABALL_PATH` variable in `score.sh`. Without DAlphaBall, packstat scoring will not function.

### CIF vs PDB
Chai-1 outputs CIF files; all other methods output PDB. The pipeline handles both formats transparently in metric extraction and Rosetta scoring. Relaxed structures are always saved as PDB regardless of input format.

### SAP Score Availability
The SAP (Spatial Aggregation Propensity) score requires `SapScoreMetric` from PyRosetta. If this class is not available in your PyRosetta version, SAP columns will contain NaN values. The `check.sh` script warns about this.

### CSV Trailing Newline
The pipeline automatically appends a trailing newline to the input CSV if one is missing. This prevents the last row from being silently dropped by `awk`.

### PyRosetta Licence
PyRosetta requires a licence. Ensure you have obtained one from the RosettaCommons before installing.

---

## Acknowledgements

The `compute_rosetta_metrics.py` script was adapted from the methodology described in:

> **Computational pipeline for scoring protein binder designs.** bioRxiv, 2025. DOI: [10.1101/2025.08.14.670059v1](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1)

The pipeline integrates the following open-source tools and models:

- [ColabFold](https://github.com/sokrypton/ColabFold) / [LocalColabFold](https://github.com/yoshitakamo/localcolabfold) — MSA generation and AlphaFold2 inference
- [Chai-1](https://github.com/chaidiscovery/chai-lab) — Diffusion-based structure prediction
- [Boltz-2](https://github.com/jwohlwend/boltz) — Diffusion-based structure prediction with potentials
- [PyRosetta](https://www.pyrosetta.org/) — Interface scoring, relaxation, and energy calculations
- [DAlphaBall](https://github.com/outpace-bio/DAlphaBall) — Packing quality assessment
- [BioPython](https://biopython.org/) — Structure parsing and analysis
