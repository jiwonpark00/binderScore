#!/usr/bin/env python3
"""
PARALLEL batch converter for single .a3m files to .aligned.pqt format for Chai-1.
Auto-detects CPU cores and adjusts workers accordingly.

Usage:
    python convert_a3m_parallel.py /path/to/a3m_files /path/to/output [--workers auto]

Example:
    python convert_a3m_parallel.py ./msas ./pqt_output --workers auto
"""

import logging
import os
from pathlib import Path
from typing import Optional, Tuple
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

import pandas as pd

# Import Chai-1 components
from chai_lab.data.parsing.msas.aligned_pqt import (
    a3m_to_aligned_dataframe,
    expected_basename,
    AlignedParquetModel,
)
from chai_lab.data.parsing.msas.data_source import MSADataSource

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ConversionTask:
    """A single conversion task."""
    a3m_path: Path
    output_dir: Path
    source_database: MSADataSource
    insert_pairing_key: bool


def get_optimal_workers(requested_workers: Optional[int] = None) -> int:
    """
    Determine optimal number of worker processes.
    
    Args:
        requested_workers: User-requested workers, or None for auto
    
    Returns:
        Number of workers to use
    """
    try:
        # Get total CPU count
        cpu_count = os.cpu_count() or 1
        
        if requested_workers is not None:
            # User specified a number
            if requested_workers <= 0:
                # Auto mode
                # Use 75% of available cores (leave some for system)
                workers = max(1, int(cpu_count * 0.75))
                logger.info(f"Auto-detected {cpu_count} CPUs, using {workers} workers (75%)")
            else:
                # Use user-specified value
                workers = min(requested_workers, cpu_count)
                if requested_workers > cpu_count:
                    logger.warning(
                        f"Requested {requested_workers} workers but only {cpu_count} CPUs available. "
                        f"Using {workers} workers."
                    )
                else:
                    logger.info(f"Using {workers} workers as requested")
        else:
            # Default: auto mode
            workers = max(1, int(cpu_count * 0.75))
            logger.info(f"Auto-detected {cpu_count} CPUs, using {workers} workers (75%)")
        
        return workers
        
    except Exception as e:
        logger.warning(f"Could not detect CPU count: {e}. Defaulting to 16 workers.")
        return 16


def convert_single_file(task: ConversionTask) -> Tuple[str, bool, str]:
    """
    Worker function to convert a single file.
    
    Returns:
        (filename, success, message)
    """
    try:
        # Convert a3m to DataFrame
        df = a3m_to_aligned_dataframe(
            task.a3m_path,
            source_database=task.source_database,
            insert_pairing_key=task.insert_pairing_key,
        )
        
        # Get query sequence
        query_seq = df.iloc[0]["sequence"]
        
        # Create output directory if needed
        task.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Save as parquet
        output_path = task.output_dir / expected_basename(query_seq)
        df.to_parquet(output_path)
        
        return task.a3m_path.name, True, str(output_path.name)
        
    except Exception as e:
        return task.a3m_path.name, False, str(e)


def parallel_batch_convert(
    input_dir: Path,
    output_dir: Path,
    source_database: MSADataSource = MSADataSource.UNIREF90,
    insert_pairing_key: bool = False,
    max_workers: Optional[int] = None,
    pattern: str = "*.a3m",
) -> dict:
    """
    Convert all .a3m files in parallel.
    
    Args:
        input_dir: Directory containing .a3m files
        output_dir: Directory to save .pqt files
        source_database: Database source for all files
        insert_pairing_key: Whether to extract pairing keys
        max_workers: Number of parallel processes (None or 0 for auto)
        pattern: Glob pattern for finding files
    
    Returns:
        Dictionary with statistics
    """
    # Find all a3m files
    a3m_files = sorted(input_dir.glob(pattern))
    
    if not a3m_files:
        logger.warning(f"No files found matching {pattern} in {input_dir}")
        return {"total": 0, "success": 0, "failed": 0}
    
    # Determine optimal number of workers
    workers = get_optimal_workers(max_workers)
    
    logger.info("="*60)
    logger.info("PARALLEL CONVERSION STARTING")
    logger.info("="*60)
    logger.info(f"Input directory:  {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Files to process: {len(a3m_files)}")
    logger.info(f"Source database:  {source_database.value}")
    logger.info(f"Pairing keys:     {'enabled' if insert_pairing_key else 'disabled (faster)'}")
    logger.info(f"Parallel workers: {workers}")
    logger.info("="*60)
    
    # Create tasks
    tasks = [
        ConversionTask(
            a3m_path=a3m_path,
            output_dir=output_dir,
            source_database=source_database,
            insert_pairing_key=insert_pairing_key,
        )
        for a3m_path in a3m_files
    ]
    
    # Process in parallel
    stats = {"total": len(tasks), "success": 0, "failed": 0}
    failed_files = []
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        # Submit all tasks
        futures = {
            executor.submit(convert_single_file, task): task 
            for task in tasks
        }
        
        # Process results as they complete
        for i, future in enumerate(as_completed(futures), 1):
            filename, success, message = future.result()
            
            if success:
                stats["success"] += 1
                if i % 10 == 0 or i == len(tasks):  # Log every 10 files
                    logger.info(f"[{i}/{len(tasks)}] ✓ {filename} -> {message}")
            else:
                stats["failed"] += 1
                failed_files.append((filename, message))
                logger.error(f"[{i}/{len(tasks)}] ✗ {filename}: {message}")
    
    # Summary
    logger.info("\n" + "="*60)
    logger.info("CONVERSION SUMMARY")
    logger.info("="*60)
    logger.info(f"Total files:    {stats['total']}")
    logger.info(f"Successful:     {stats['success']}")
    logger.info(f"Failed:         {stats['failed']}")
    logger.info(f"Success rate:   {stats['success']/stats['total']*100:.1f}%")
    
    if failed_files:
        logger.info("\nFAILED FILES:")
        for filename, error in failed_files[:10]:  # Show first 10
            logger.info(f"  - {filename}: {error}")
        if len(failed_files) > 10:
            logger.info(f"  ... and {len(failed_files) - 10} more")
    
    logger.info("="*60)
    
    return stats


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Parallel converter for .a3m to .aligned.pqt files"
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing .a3m files"
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory to save .pqt files"
    )
    parser.add_argument(
        "--workers",
        type=str,
        default="auto",
        help="Number of parallel workers ('auto' or number, default: auto)"
    )
    parser.add_argument(
        "--source",
        type=str,
        default="uniref90",
        choices=["uniref90", "uniprot", "mgnify", "bfd_uniclust"],
        help="Source database (default: uniref90)"
    )
    parser.add_argument(
        "--with-pairing-keys",
        action="store_true",
        help="Extract pairing keys (slower, only needed for multimers)"
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.a3m",
        help="Glob pattern for finding files (default: *.a3m)"
    )
    
    args = parser.parse_args()
    
    # Validate input
    if not args.input_dir.exists():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    if not args.input_dir.is_dir():
        logger.error(f"Input path is not a directory: {args.input_dir}")
        sys.exit(1)
    
    # Parse workers argument
    if args.workers.lower() == "auto":
        max_workers = 0  # Auto mode
    else:
        try:
            max_workers = int(args.workers)
            if max_workers < 0:
                logger.error("Workers must be positive or 'auto'")
                sys.exit(1)
        except ValueError:
            logger.error(f"Invalid workers value: {args.workers}. Use 'auto' or a number.")
            sys.exit(1)
    
    # Get source database
    source_database = MSADataSource(args.source)
    
    # Run conversion
    stats = parallel_batch_convert(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        source_database=source_database,
        insert_pairing_key=args.with_pairing_keys,
        max_workers=max_workers if max_workers > 0 else None,
        pattern=args.pattern,
    )
    
    # Exit with error if any failed
    if stats["failed"] > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()