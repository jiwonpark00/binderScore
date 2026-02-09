#!/usr/bin/env python3
import os
import re
import csv
import glob
import time
import argparse
import multiprocessing as mp
from typing import Dict, List, Tuple

# ---------------------------
# PyRosetta
import pyrosetta as pr
from pyrosetta import rosetta
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, OrResidueSelector
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.io import pose_from_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
# NEW: Import for CIF support
from pyrosetta.rosetta.core.import_pose import pose_from_file, FileType
try:
    from pyrosetta.rosetta.core.pack.guidance_scoreterms.sap import SapScoreMetric
    HAS_SAP = True
except Exception:
    HAS_SAP = False

# ---------------------------
# BioPython / SciPy
from Bio.PDB import PDBParser, MMCIFParser, Selection, Polypeptide
from Bio.PDB.Selection import unfold_entities
import numpy as np
from scipy.spatial import cKDTree

# ---------------------------
# I/O helpers
import pandas as pd

# ============================================================
# Utils
# ============================================================

def get_file_extension(filepath):
    """Get the file extension (pdb or cif)."""
    return os.path.splitext(filepath)[1].lower()

def load_pose(filepath):
    """
    Universal pose loader that handles both PDB and CIF formats.
    
    Args:
        filepath: Path to structure file (.pdb or .cif)
    
    Returns:
        pyrosetta.Pose object
    """
    ext = get_file_extension(filepath)
    
    if ext == '.pdb':
        # Use pose_from_file with PDB type (or let it auto-detect)
        return pose_from_file(filepath, False, FileType.PDB_file)
    elif ext == '.cif':
        # Use pose_from_file with CIF type
        return pose_from_file(filepath, False, FileType.CIF_file)
    else:
        # Fallback: try auto-detection
        print(f"[warn] Unknown file extension {ext}, attempting auto-detection")
        return pose_from_file(filepath)

def parse_structure(filepath):
    """
    Universal BioPython structure parser that handles both PDB and CIF formats.
    
    Args:
        filepath: Path to structure file (.pdb or .cif)
    
    Returns:
        Bio.PDB.Structure object
    """
    ext = get_file_extension(filepath)
    
    if ext == '.pdb':
        parser = PDBParser(QUIET=True)
    elif ext == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        # Default to PDB parser
        print(f"[warn] Unknown file extension {ext}, defaulting to PDB parser")
        parser = PDBParser(QUIET=True)
    
    return parser.get_structure("complex", filepath)

def clean_structure(structure_file):
    """
    Clean structure file by keeping only relevant lines.
    For PDB files, keeps ATOM, HETATM, MODEL, TER, END lines.
    For CIF files, no cleaning is performed (CIF files are typically already clean).
    """
    ext = get_file_extension(structure_file)
    
    if ext == '.pdb':
        with open(structure_file, 'r') as f_in:
            relevant = [ln for ln in f_in if ln.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END'))]
        with open(structure_file, 'w') as f_out:
            f_out.writelines(relevant)
    elif ext == '.cif':
        # CIF files don't need cleaning
        pass

def get_binder_id(structure_path):
    """Get binder ID from structure file path (works for both .pdb and .cif)."""
    return os.path.splitext(os.path.basename(structure_path))[0]

# ============================================================
# Hotspots (distance-based)
# ============================================================

def _hotspots_dict(structure_file, chain, partner_chain, cutoff=4.0) -> Dict[int,str]:
    """
    Compute hotspot residues using BioPython (supports both PDB and CIF).
    """
    structure = parse_structure(structure_file)
    chain_obj = structure[0][chain]
    partner_obj = structure[0][partner_chain]
    chain_atoms = Selection.unfold_entities(chain_obj, 'A')
    partner_atoms = Selection.unfold_entities(partner_obj, 'A')

    chain_coords = np.array([atom.coord for atom in chain_atoms])
    partner_coords = np.array([atom.coord for atom in partner_atoms])
    if len(chain_coords) == 0 or len(partner_coords) == 0:
        return {}

    chain_tree = cKDTree(chain_coords)
    partner_tree = cKDTree(partner_coords)
    pairs = chain_tree.query_ball_tree(partner_tree, cutoff)

    hotspots = {}
    for idx, partner_indices in enumerate(pairs):
        if partner_indices:
            residue = chain_atoms[idx].get_parent()
            resnum = residue.id[1]
            if resnum not in hotspots and residue.get_resname() in Polypeptide.standard_aa_names:
                hotspots[resnum] = f"{residue.get_resname()}{resnum}"
    return hotspots

def compute_interface_hotspots(structure_file, binder_chain="A", target_chain="B", cutoff=4.0) -> Tuple[str,str]:
    """Compute interface hotspots (supports both PDB and CIF)."""
    b = _hotspots_dict(structure_file, binder_chain, target_chain, cutoff)
    t = _hotspots_dict(structure_file, target_chain, binder_chain, cutoff)
    binder_str = ":".join([b[k] for k in sorted(b.keys())])
    target_str = ":".join([t[k] for k in sorted(t.keys())])
    return binder_str, target_str

# ============================================================
# SAP (optional) and PyRosetta metrics
# ============================================================

def compute_delta_sap(structure_path, binder_chains='A'):
    """Compute SAP scores (supports both PDB and CIF)."""
    if not HAS_SAP:
        return None, None, None
    pose = load_pose(structure_path)
    pdb_info = pose.pdb_info()
    all_chains = sorted({pdb_info.chain(i) for i in range(1, pose.total_residue()+1)})
    binder_list = list(binder_chains)
    target_list = [c for c in all_chains if c not in binder_list]

    sap_metric = SapScoreMetric()

    def sap_for(chains):
        sel = OrResidueSelector()
        for c in chains:
            sel.add_residue_selector(ChainSelector(c))
        sap_metric.set_score_selector(sel)
        sap_metric.set_sap_calculate_selector(sel)
        sap_metric.set_sasa_selector(sel)
        return sap_metric.calculate(pose)

    sap_binder  = sap_for(binder_list)
    sap_target  = sap_for(target_list)
    sap_complex = sap_for(binder_list + target_list)
    delta_sap   = sap_complex - (sap_binder + sap_target)
    return sap_target, sap_binder, delta_sap

def score_interface(structure_file, binder_chain="A", target_chain="B", binder_hotspots=None, target_hotspots=None):
    """Score interface using PyRosetta (supports both PDB and CIF)."""
    pose = load_pose(structure_file)

    iam = InterfaceAnalyzerMover()
    iam.set_interface(f"{binder_chain}_{target_chain}")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    if binder_hotspots is None or target_hotspots is None:
        binder_hotspots, target_hotspots = compute_interface_hotspots(
            structure_file, binder_chain, target_chain
        )

    # Count binder interface residue types
    standard = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
    interface_AA = {aa: 0 for aa in standard}
    binder_hotspots_list = binder_hotspots.split(":") if binder_hotspots else []
    target_hotspots_list = target_hotspots.split(":") if target_hotspots else []
    for h in binder_hotspots_list:
        rt = h[:3]
        if rt in interface_AA:
            interface_AA[rt] += 1

    interface_nres = len(binder_hotspots_list)
    interface_nres_target = len(target_hotspots_list)

    hydrophobic = {'ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'PRO', 'VAL', 'TRP', 'TYR'}
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic)
    interface_hydrophobicity = (hydrophobic_count / interface_nres * 100) if interface_nres else 0

    data = iam.get_all_data()
    interface_sc = data.sc_value
    interface_interface_hbonds = data.interface_hbonds
    interface_dG = iam.get_interface_dG()
    interface_dSASA = iam.get_interface_delta_sasa()
    interface_packstat = iam.get_interface_packstat()
    interface_dG_SASA_ratio = data.dG_dSASA_ratio * 100

    buns_filter = XmlObjects.static_get_filter(
        '<BuriedUnsatHbonds2 name="Buns2" jump_number="1" cutoff="20" generous_hbonds="true" scorefxn="scorefxn" />'
    )
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    # binder energy and surface SASA fraction
    sel_A = ChainSelector(binder_chain)
    tem = rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(sel_A)
    binder_score = tem.calculate(pose)

    bsasa = rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(sel_A)
    binder_sasa = bsasa.calculate(pose)
    interface_binder_fraction = (interface_dSASA / binder_sasa * 100) if binder_sasa > 0 else 0

    # surface hydrophobicity (binder)
    layer_sel = rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res = layer_sel.apply(pose)
    exp_apol_count = 0
    total_count = 0
    pdb_info = pose.pdb_info()
    for i in range(1, len(surface_res) + 1):
        if surface_res[i] and pdb_info.chain(i) == binder_chain:
            res = pose.residue(i)
            if res.is_apolar() or res.name() in ['PHE','TRP','TYR']:
                exp_apol_count += 1
            total_count += 1
    surface_hydrophobicity = (exp_apol_count / total_count) if total_count else 0

    # SAP (optional)
    try:
        sap_target, sap_binder, delta_sap = compute_delta_sap(structure_file, binder_chain)
    except Exception:
        sap_target, sap_binder, delta_sap = None, None, None

    interface_scores = {
        'rosetta_binder_score': binder_score,
        'rosetta_surface_hydrophobicity': surface_hydrophobicity,
        'rosetta_interface_sc': interface_sc,
        'rosetta_interface_packstat': interface_packstat,
        'rosetta_interface_dG': interface_dG,
        'rosetta_interface_dSASA': interface_dSASA,
        'rosetta_interface_dG_dSASA_ratio': interface_dG_SASA_ratio,
        'rosetta_interface_percent_of_binder_in_intf': interface_binder_fraction,
        'rosetta_interface_binder_hydrophobicity_percent': interface_hydrophobicity,
        'rosetta_interface_nres_binder': interface_nres,
        'rosetta_interface_nres_target': interface_nres_target,
        'rosetta_interface_interface_hbonds': interface_interface_hbonds,
        'rosetta_interface_hbond_percentage': interface_hbond_percentage,
        'rosetta_interface_unsat_hbonds': interface_delta_unsat_hbonds,
        'rosetta_interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage,
        'sap_binder': sap_binder,
        'sap_target': sap_target,
        'sap_delta': delta_sap,
    }
    # round floats for readability
    interface_scores = {k: (round(v, 2) if isinstance(v, float) else v) for k, v in interface_scores.items()}
    return interface_scores, interface_AA, binder_hotspots

# ============================================================
# Relax
# ============================================================

def pr_relax(structure_file, relaxed_pdb_path):
    """Relax structure using PyRosetta (supports both PDB and CIF input)."""
    if os.path.exists(relaxed_pdb_path):
        return
    pose = load_pose(structure_file)
    start_pose = pose.clone()

    mm = MoveMap()
    mm.set_chi(True)
    mm.set_bb(True)
    mm.set_jump(False)

    fr = FastRelax()
    fr.set_scorefxn(pr.get_fa_scorefxn())
    fr.set_movemap(mm)
    fr.max_iter(200)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.constrain_relax_to_start_coords(True)
    fr.apply(pose)

    # Align to starting pose (chain 1 by index)
    aln = AlignChainMover()
    aln.source_chain(0)
    aln.target_chain(0)
    aln.pose(start_pose)
    aln.apply(pose)

    # copy B factors
    for resid in range(1, pose.total_residue() + 1):
        if pose.residue(resid).is_protein():
            b = start_pose.pdb_info().bfactor(resid, 1)
            for atom_id in range(1, pose.residue(resid).natoms() + 1):
                pose.pdb_info().bfactor(resid, atom_id, b)

    # Always output as PDB for relaxed structures
    pose.dump_pdb(relaxed_pdb_path)
    clean_structure(relaxed_pdb_path)

# ============================================================
# Per-file worker
# ============================================================

def worker_init(dalphaball_path: str):
    os.environ["DALPHABALL"] = dalphaball_path
    pr.init(
        f"-ignore_unrecognized_res 1 "
        f"-detect_disulf 0 "
        f"-holes:dalphaball {dalphaball_path} "
        f"-mute all",
        silent=True
    )

def process_one(structure_file, relaxed_dir, binder_chain, target_chain):
    """Process one structure file (PDB or CIF)."""
    binder_id = get_binder_id(structure_file)
    ext = get_file_extension(structure_file)
    
    # Relaxed file is always saved as PDB
    relaxed_pdb_path = os.path.join(relaxed_dir, f"{binder_id}.pdb")
    
    try:
        pr_relax(structure_file, relaxed_pdb_path)
        binder_hotspots, target_hotspots = compute_interface_hotspots(
            relaxed_pdb_path, binder_chain, target_chain
        )
        metrics, interface_AA, _ = score_interface(
            relaxed_pdb_path, binder_chain, target_chain, 
            binder_hotspots, target_hotspots
        )
        data = {
            "metrics": metrics,
            "binder_hotspots": binder_hotspots,
            "target_hotspots": target_hotspots,
        }
        print(f"[ok] {binder_id}{ext}")
        return binder_id, data
    except Exception as e:
        print(f"[fail] {binder_id}{ext}: {e}")
        return binder_id, {
            "metrics": {},
            "binder_hotspots": "",
            "target_hotspots": "",
        }

def merge_metrics_into_run(run_csv: str, results: Dict[str, Dict], prefix: str):
    """
    If run_csv is provided, read it, append/overwrite columns with '{prefix}_*' per binder_id, and write back.
    """
    if not run_csv:
        return  # nothing to merge

    if not os.path.isfile(run_csv):
        raise SystemExit(f"run_csv not found: {run_csv}")
    df_run = pd.read_csv(run_csv)
    if "binder" not in df_run.columns:
        raise SystemExit("run_csv must contain a 'binder' column.")

    df_run = df_run.set_index("binder")

    # build new columns DataFrame
    rows = []
    for binder_id, data in results.items():
        row = {"binder": binder_id}
        for k, v in data["metrics"].items():
            row[f"{prefix}_{k}"] = v
        row[f"{prefix}_distbased_intf_res_binder"] = data.get("binder_hotspots", "")
        row[f"{prefix}_distbased_intf_res_target"] = data.get("target_hotspots", "")
        rows.append(row)
    if not rows:
        print("[merge] no results to merge.")
        return

    df_new = pd.DataFrame(rows).set_index("binder")

    # only merge for binder_ids that exist in run.csv; warn about extras
    missing = [bid for bid in df_new.index if bid not in df_run.index]
    if missing:
        print(f"[merge] WARNING: {len(missing)} binder(s) not found in run.csv; skipping: "
              f"{missing[:5]}{'...' if len(missing)>5 else ''}")
        df_new = df_new[df_new.index.isin(df_run.index)]

    # add/overwrite columns
    for col in df_new.columns:
        df_run[col] = df_new[col]

    # write back
    df_run.reset_index().to_csv(run_csv, index=False)
    print(f"[merge] updated: {run_csv} (added/overwrote {len(df_new.columns)} columns)")


# ============================================================
# CLI
# ============================================================

import argparse, glob, multiprocessing as mp, os


def parse_folder_kv(value):
    """
    Parse --folder flags of the form name:path
    """
    if ":" not in value:
        raise argparse.ArgumentTypeError("Folder must be given as name:path")
    name, path = value.split(":", 1)
    name = name.strip()
    path = os.path.expandvars(os.path.expanduser(path.strip()))
    if not name:
        raise argparse.ArgumentTypeError("Folder name is empty")
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Folder path does not exist: {path}")
    return name, path

def _parse_inputs(argv=None):
    ap = argparse.ArgumentParser(
        description="Compute PyRosetta interface metrics for one or many folders and optionally merge into run.csv.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Repeatable inputs (now using --folder name:path style)
    ap.add_argument(
        "--folder",
        action="append",
        required=True,
        type=parse_folder_kv,
        help='Repeatable. Format "name:path". Example: --folder af3:/path/to/AF3/pdbs',
    )

    # Common
    ap.add_argument("--run-csv", default=None,
                    help="Optional: path to an existing run.csv to update in-place")
    ap.add_argument("--out-csv", default="rosetta_metrics.csv",
                    help="Output CSV filename to write aggregated results (always written)")
    ap.add_argument("--relaxed-dir", default=None,
                    help="Base folder to cache relaxed PDBs. If multiple inputs, a subfolder per prefix is created. "
                         "Default: <pdbs>/relaxed_pdbs for each input.")
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--target-chain", default="B")
    ap.add_argument("--nprocs", type=int, default=mp.cpu_count())
    ap.add_argument("--dalphaball-path", default="./functions/DAlphaBall.gcc")
    ap.add_argument("--update-runcsv",action="store_true",help="If set, update the run.csv in-place. Default: only write per-model CSVs.")

    args = ap.parse_args(argv)

    # Convert to dict mapping
    mapping = dict(args.folder)  # [(name, path), ...] → {name: path}

    # Require at least one input
    if not mapping:
        raise SystemExit("Provide one or more --folder NAME:PATH.")

    return args, mapping

def _relaxed_dir_for(args, structure_dir, prefix, num_inputs):
    if args.relaxed_dir:
        return os.path.join(args.relaxed_dir, prefix) if num_inputs > 1 else args.relaxed_dir
    return os.path.join(structure_dir, "relaxed_pdbs")

def _score_one_folder(structure_dir, prefix, run_csv, relaxed_dir, binder_chain, target_chain, nprocs, dalphaball_path, update_runcsv=False):
    """Score all structure files (PDB and CIF) in a folder."""
    # Find both PDB and CIF files
    pdb_files = sorted(glob.glob(os.path.join(structure_dir, "*.pdb")))
    cif_files = sorted(glob.glob(os.path.join(structure_dir, "*.cif")))
    structure_files = pdb_files + cif_files
    
    if not structure_files:
        print(f"[warn] No structure files (PDB/CIF) found in {structure_dir} (prefix={prefix}); skipping.")
        return {}

    print(f"[info] scoring {len(structure_files)} structure files from {structure_dir} (prefix={prefix})")
    print(f"[info]   - {len(pdb_files)} PDB files")
    print(f"[info]   - {len(cif_files)} CIF files")
    print(f"[info] relaxed pdb cache: {relaxed_dir}")
    os.makedirs(relaxed_dir, exist_ok=True)

    tasks = [(sf, relaxed_dir, binder_chain, target_chain) for sf in structure_files]
    results: Dict[str, Dict] = {}

    if nprocs == 1:
        worker_init(dalphaball_path)
        for t in tasks:
            bid, data = process_one(*t)
            results[bid] = data
    else:
        with mp.Pool(nprocs, initializer=worker_init, initargs=(dalphaball_path,)) as pool:
            for bid, data in pool.starmap(process_one, tasks):
                results[bid] = data

        # Optionally merge into run.csv (only if --update-runcsv set)
    if update_runcsv and run_csv:
        merge_metrics_into_run(run_csv, results, prefix)
    elif run_csv and not update_runcsv:
        print(f"[info] Skipping merge into run.csv for prefix={prefix} (set --update-runcsv to enable)")


    # Return results for global CSV
    return results

def write_rosetta_metrics(all_results_by_prefix: Dict[str, Dict[str, Dict]], out_csv: str):
    """
    Build a wide table with one row per binder_id and columns named {prefix}_{metric}.
    Always includes distbased hotspots as {prefix}_distbased_intf_res_binder/target.
    """
    # Collect union of binder_ids
    all_binder_ids = set()
    for prefix, res in all_results_by_prefix.items():
        all_binder_ids.update(res.keys())
    if not all_binder_ids:
        print(f"[write] No results to write to {out_csv}.")
        return

    rows = []
    for bid in sorted(all_binder_ids):
        row = {"binder": bid}
        for prefix, res in all_results_by_prefix.items():
            if bid not in res:
                continue
            data = res[bid]
            for k, v in data["metrics"].items():
                row[f"{prefix}_{k}"] = v
            row[f"{prefix}_distbased_intf_res_binder"] = data.get("binder_hotspots", "")
            row[f"{prefix}_distbased_intf_res_target"] = data.get("target_hotspots", "")
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"[write] wrote {out_csv} with {len(df)} rows and {len(df.columns)} columns.")

def main():
    args, mapping = _parse_inputs()
    num_inputs = len(mapping)

    all_results_by_prefix: Dict[str, Dict[str, Dict]] = {}
    for prefix, structure_dir in mapping.items():
        rdir = _relaxed_dir_for(args, structure_dir, prefix, num_inputs)
        res = _score_one_folder(
            structure_dir=structure_dir,
            prefix=prefix,
            run_csv=args.run_csv,
            relaxed_dir=rdir,
            binder_chain=args.binder_chain,
            target_chain=args.target_chain,
            nprocs=args.nprocs,
            dalphaball_path=args.dalphaball_path,
            update_runcsv=args.update_runcsv,
        )
        all_results_by_prefix[prefix] = res or {}

    # Always write aggregated CSV (default: rosetta_metrics.csv)
    write_rosetta_metrics(all_results_by_prefix, args.out_csv)
    print("[done] all inputs processed.")
    
if __name__ == "__main__":
    main()