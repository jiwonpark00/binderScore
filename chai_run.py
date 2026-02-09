#!/usr/bin/env python3
"""
chai_run.py - Run Chai-1 inference with optional MSA support

Usage:
    python chai_run.py \
        --fasta <fasta_file> \
        --output_dir <output_dir> \
        --gpu <gpu_id> \
        --num_diffn_samples 5 \
        --num_trunk_recycles 3 \
        --num_diffn_timesteps 200 \
        [--msa_directory <msa_directory>]

Outputs:
    - pred.model_idx_{N}.cif for each sample
    - scores.model_idx_{N}.npz for each sample
    - run.log with inference details
"""

import argparse
import random
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, stream=sys.stderr)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Run Chai-1 inference with optional MSA support"
    )
    parser.add_argument("--fasta", required=True, type=Path,
                        help="Path to input FASTA file")
    parser.add_argument("--output_dir", required=True, type=Path,
                        help="Directory for output files")
    parser.add_argument("--gpu", required=True,
                        help="GPU device ID (e.g., 0, 1, 2)")
    parser.add_argument("--num_diffn_samples", type=int, default=5,
                        help="Number of diffusion samples to generate (default: 5)")
    parser.add_argument("--num_trunk_recycles", type=int, default=3,
                        help="Number of trunk recycles (default: 3)")
    parser.add_argument("--num_diffn_timesteps", type=int, default=200,
                        help="Number of diffusion timesteps (default: 200)")
    parser.add_argument("--msa_directory", type=Path, default=None,
                        help="Optional: directory containing .aligned.pqt MSA files")

    args = parser.parse_args()

    device = f"cuda:{args.gpu}"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate random seed for each run to ensure different results
    random_seed = random.randint(0, 2**32 - 1)

    # Import chai_lab after setting up paths
    from chai_lab.chai1 import run_inference

    # Build inference arguments
    inference_kwargs = {
        "fasta_file": args.fasta,
        "output_dir": args.output_dir,
        "num_trunk_recycles": args.num_trunk_recycles,
        "num_diffn_timesteps": args.num_diffn_timesteps,
        "num_diffn_samples": args.num_diffn_samples,
        "seed": random_seed,
        "device": device,
        "low_memory": False,
    }

    # Add MSA directory if provided
    if args.msa_directory and args.msa_directory.is_dir():
        inference_kwargs["msa_directory"] = args.msa_directory
        logger.info(f"Using MSA directory: {args.msa_directory}")
    else:
        logger.info("Running without pre-computed MSAs (single-sequence mode or MSA server)")

    # Run inference
    logger.info(f"Starting Chai-1 inference for {args.fasta.name}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Number of samples: {args.num_diffn_samples}")
    logger.info(f"Trunk recycles: {args.num_trunk_recycles}")
    logger.info(f"Diffusion timesteps: {args.num_diffn_timesteps}")
    logger.info(f"Random seed: {random_seed}")

    candidates = run_inference(**inference_kwargs)

    logger.info(f"Inference complete. Generated {len(candidates.ranking_data)} samples")
    logger.info(f"Output files written to: {args.output_dir}")


if __name__ == "__main__":
    main()
