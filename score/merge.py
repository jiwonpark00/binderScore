#!/usr/bin/env python3
"""
merge.py - Merge metrics from multiple structure prediction methods

This script merges metrics CSVs from different prediction methods (AF2-multimer,
AF2-ptm complex, Boltz, Chai, Binder) into a single binderScore.csv file with
proper column organization.

Usage:
    python merge.py <project_dir>
    python merge.py <project_dir> --output custom_scores.csv
    python merge.py <project_dir> --metrics-dir /custom/metrics

Author: Generated for binder screening pipeline
"""

import pandas as pd
import os
import sys
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple


# Configuration constants
COMPLEX_METHODS = ['af2_multimer', 'af2_ptm_complex', 'boltz', 'chai']
PREFIX_MAP = {
    'af2_multimer': 'mul',
    'af2_ptm_complex': 'mon',
    'boltz': 'bol',
    'chai': 'cha',
    'binder': 'bin'
}
PREFIX_ORDER = ['mul', 'mon', 'bol', 'cha']
PIVOT_COLS = ['status', 'attempt', 'num_passing', 'min_distance']
METRICS_FILES = {
    'af2_multimer': 'af2_multimer_metrics.csv',
    'af2_ptm_complex': 'af2_ptm_complex_metrics.csv',
    'boltz': 'boltz_metrics.csv',
    'chai': 'chai_metrics.csv',
    'binder': 'binder_metrics.csv'
}


def validate_inputs(project_dir: Path, metrics_dir: Path) -> Dict[str, Path]:
    """
    Validate input directories and files.
    
    Args:
        project_dir: Path to project directory
        metrics_dir: Path to metrics directory
        
    Returns:
        Dictionary with validated paths
        
    Raises:
        FileNotFoundError: If required files/directories are missing
        ValueError: If no metrics files are found
    """
    if not project_dir.exists():
        raise FileNotFoundError(f"Project directory not found: {project_dir}")
    
    summary_path = project_dir / "summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.csv not found: {summary_path}")
    
    if not metrics_dir.exists():
        raise FileNotFoundError(f"Metrics directory not found: {metrics_dir}")
    
    # Check if at least one metrics file exists
    metrics_found = False
    for filename in METRICS_FILES.values():
        if (metrics_dir / filename).exists():
            metrics_found = True
            break
    
    if not metrics_found:
        raise ValueError(f"No metrics CSV files found in: {metrics_dir}")
    
    return {
        'project_dir': project_dir,
        'metrics_dir': metrics_dir,
        'summary_path': summary_path
    }


def pivot_summary(summary_path: Path, complex_methods: List[str], 
                  pivot_cols: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    """
    Pivot summary.csv to create wide-format columns for each method.
    
    Args:
        summary_path: Path to summary.csv
        complex_methods: List of complex method names
        pivot_cols: List of columns to pivot
        
    Returns:
        Tuple of (pivoted DataFrame, ordered column names)
    """
    print("[INFO] Pivoting summary.csv...")
    summary_df = pd.read_csv(summary_path)
    
    # Filter to only complex methods (binder method not in summary)
    summary_df = summary_df[summary_df['method'].isin(complex_methods)]
    
    # Pivot each column separately
    pivoted_dfs = []
    for col in pivot_cols:
        pivot = summary_df.pivot(index='binder', columns='method', values=col)
        pivot.columns = [f"{method}_{col}" for method in pivot.columns]
        pivoted_dfs.append(pivot)
    
    # Concatenate all pivoted columns
    pivoted_summary = pd.concat(pivoted_dfs, axis=1).reset_index()
    
    # Reorder pivoted columns: status, attempt, num_passing, min_distance
    # (each grouped by all methods)
    ordered_pivot_cols = ['binder']
    for col_type in pivot_cols:
        for method in complex_methods:
            col_name = f"{method}_{col_type}"
            if col_name in pivoted_summary.columns:
                ordered_pivot_cols.append(col_name)
    
    pivoted_summary = pivoted_summary[[c for c in ordered_pivot_cols 
                                       if c in pivoted_summary.columns]]
    
    print(f"[INFO] Pivoted summary: {len(pivoted_summary)} binders")
    
    return pivoted_summary, ordered_pivot_cols


def load_metrics_files(metrics_dir: Path, metrics_files: Dict[str, str], 
                       prefix_map: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """
    Load all available metrics CSV files.
    
    Args:
        metrics_dir: Path to metrics directory
        metrics_files: Dictionary mapping method names to filenames
        prefix_map: Dictionary mapping method names to column prefixes
        
    Returns:
        Dictionary of DataFrames keyed by method name
    """
    print("[INFO] Loading metrics CSVs...")
    metrics_dfs = {}
    
    for method, filename in metrics_files.items():
        filepath = metrics_dir / filename
        if filepath.exists():
            df = pd.read_csv(filepath)
            prefix = prefix_map[method]
            
            # Rename sample_name column to 'binder'
            sample_col = f"{prefix}_sample_name"
            if sample_col in df.columns:
                df = df.rename(columns={sample_col: 'binder'})
            
            metrics_dfs[method] = df
            print(f"[INFO] Loaded {filename}: {len(df)} rows, {len(df.columns)} columns")
        else:
            print(f"[WARN] {filename} not found, skipping")
    
    return metrics_dfs


def merge_dataframes(pivoted_summary: pd.DataFrame, 
                    metrics_dfs: Dict[str, pd.DataFrame],
                    complex_methods: List[str]) -> pd.DataFrame:
    """
    Merge all dataframes into a single DataFrame.
    
    Args:
        pivoted_summary: Pivoted summary DataFrame
        metrics_dfs: Dictionary of metrics DataFrames
        complex_methods: List of complex method names
        
    Returns:
        Merged DataFrame
    """
    print("[INFO] Merging dataframes...")
    
    # Start with pivoted summary
    merged_df = pivoted_summary.copy()
    
    # Merge binder metrics first
    if 'binder' in metrics_dfs:
        merged_df = merged_df.merge(metrics_dfs['binder'], on='binder', how='left')
    
    # Merge complex method metrics
    for method in complex_methods:
        if method in metrics_dfs:
            merged_df = merged_df.merge(metrics_dfs[method], on='binder', how='left')
    
    print(f"[INFO] Merged dataframe: {len(merged_df)} rows, {len(merged_df.columns)} columns")
    
    return merged_df


def reorder_columns(merged_df: pd.DataFrame, ordered_pivot_cols: List[str],
                   prefix_order: List[str]) -> pd.DataFrame:
    """
    Reorder columns in the merged DataFrame for logical grouping.
    
    Column order:
    1. binder (ID)
    2. Pivoted summary columns
    3. Binder metrics (bin_*)
    4. Complex method metrics (grouped by metric suffix, ordered by prefix)
    
    Args:
        merged_df: Merged DataFrame
        ordered_pivot_cols: Ordered list of pivoted column names
        prefix_order: Order of prefixes for complex metrics
        
    Returns:
        DataFrame with reordered columns
    """
    print("[INFO] Reordering columns...")
    
    all_cols = list(merged_df.columns)
    
    # Group 1: binder
    id_cols = ['binder']
    
    # Group 2-5: Pivoted columns (already ordered in pivoted_summary)
    pivot_col_names = [c for c in ordered_pivot_cols 
                      if c != 'binder' and c in all_cols]
    
    # Group 6: Binder metrics (bin_* columns)
    binder_cols = [c for c in all_cols if c.startswith('bin_')]
    
    # Group 7: Complex method metrics (grouped by metric suffix)
    complex_metric_cols = [c for c in all_cols 
                          if any(c.startswith(p + '_') for p in prefix_order) 
                          and c not in pivot_col_names]
    
    # Extract metric suffixes and group
    def get_suffix(col: str) -> str:
        """Extract metric suffix from column name."""
        for prefix in prefix_order:
            if col.startswith(prefix + '_'):
                return col[len(prefix) + 1:]
        return None
    
    # Build suffix -> columns mapping
    suffix_to_cols = defaultdict(list)
    for col in complex_metric_cols:
        suffix = get_suffix(col)
        if suffix:
            suffix_to_cols[suffix].append(col)
    
    # Sort suffixes alphabetically, then order columns within each suffix by PREFIX_ORDER
    grouped_complex_cols = []
    for suffix in sorted(suffix_to_cols.keys()):
        cols_for_suffix = suffix_to_cols[suffix]
        # Sort by prefix order
        cols_for_suffix.sort(
            key=lambda c: next((i for i, p in enumerate(prefix_order) 
                              if c.startswith(p + '_')), 999)
        )
        grouped_complex_cols.extend(cols_for_suffix)
    
    # Final column order
    final_cols = id_cols + pivot_col_names + binder_cols + grouped_complex_cols
    
    # Ensure no duplicates and all columns are included
    final_cols_set = set(final_cols)
    missing_cols = [c for c in all_cols if c not in final_cols_set]
    if missing_cols:
        print(f"[WARN] Adding {len(missing_cols)} unclassified columns at the end")
        final_cols.extend(missing_cols)
    
    # Reorder
    merged_df = merged_df[[c for c in final_cols if c in merged_df.columns]]
    
    return merged_df


def main():
    """Main execution function."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Merge metrics from multiple structure prediction methods',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python merge.py /path/to/project
  python merge.py /path/to/project --output custom_scores.csv
  python merge.py /path/to/project --metrics-dir /custom/metrics
        """
    )
    
    parser.add_argument(
        'project_dir',
        type=Path,
        help='Project directory path'
    )
    
    parser.add_argument(
        '--metrics-dir',
        type=Path,
        default=None,
        help='Metrics directory (default: <project_dir>/metrics)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='binderScore.csv',
        help='Output filename (default: binderScore.csv)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    # Set metrics directory
    metrics_dir = args.metrics_dir if args.metrics_dir else args.project_dir / 'metrics'
    
    try:
        # Step 0: Validate inputs
        paths = validate_inputs(args.project_dir, metrics_dir)
        
        # Step 1: Pivot summary.csv
        pivoted_summary, ordered_pivot_cols = pivot_summary(
            paths['summary_path'], 
            COMPLEX_METHODS, 
            PIVOT_COLS
        )
        
        # Step 2: Load metrics CSVs
        metrics_dfs = load_metrics_files(metrics_dir, METRICS_FILES, PREFIX_MAP)
        
        # Step 3: Merge all dataframes
        merged_df = merge_dataframes(pivoted_summary, metrics_dfs, COMPLEX_METHODS)
        
        # Step 4: Reorder columns
        merged_df = reorder_columns(merged_df, ordered_pivot_cols, PREFIX_ORDER)
        
        # Step 5: Save output
        output_path = metrics_dir / args.output
        merged_df.to_csv(output_path, index=False)
        
        print(f"[INFO] Saved {args.output}: {len(merged_df)} rows, {len(merged_df.columns)} columns")
        print(f"[INFO] Output: {output_path}")
        
        return 0
        
    except FileNotFoundError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())