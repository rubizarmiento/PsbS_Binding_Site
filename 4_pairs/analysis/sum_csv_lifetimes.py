#!/usr/bin/env python3
"""
Sum lifetime values across multiple CSV files.

This script reads multiple CSV files with a common prefix from a directory,
extracts the 'resid' and 'sum_ns' columns, and sums the 'sum_ns' values 
for each residue across all files. The output contains only 'resid' and 'sum_ns'.

Usage:
------
python sum_csv_lifetimes.py -d <directory> -o <output_csv> -prefix <file_prefix>

Arguments:
----------
-d, --dir : str
    Directory containing the CSV files to sum.
-o, --output : str
    Output CSV file path.
-prefix, --prefix : str
    Prefix of CSV files to include (e.g., 'psbs_' will match 'psbs_*.csv').

Input CSV Format:
-----------------
Each CSV file must contain at least:
- resid : int - Residue ID
- sum_ns : float - Lifetime value in nanoseconds

Output CSV Format:
------------------
- resid : int - Residue ID
- sum_ns : float - Sum of lifetime values across all input files

Examples:
---------
# Sum all psbs_*.csv files
python sum_csv_lifetimes.py -d ./lifetimes -o ./summary.csv -prefix psbs_

# Sum all chain_4_*.csv files
python sum_csv_lifetimes.py -d ./lifetimes -o ./chain4_sum.csv -prefix chain_4_

Author: Rubi Zarmiento-Garcia
"""

import os
import sys
import argparse
import glob
import pandas as pd


def find_csv_files(directories, prefix):
    """
    Find all CSV files matching the prefix pattern across multiple directories.
    
    Parameters
    ----------
    directories : str or list
        Directory or list of directories to search for CSV files.
    prefix : str
        Prefix pattern to match (e.g., 'psbs_').
        
    Returns
    -------
    list
        List of matching CSV file paths.
    """
    # Convert single directory to list
    if isinstance(directories, str):
        directories = [directories]
    
    all_files = []
    
    for directory in directories:
        if not os.path.isdir(directory):
            print(f"WARNING: Directory not found: {directory}")
            continue
            
        pattern = os.path.join(directory, f"{prefix}*residue_summary_df.csv")
        files = glob.glob(pattern)
        
        if not files:
            # Try without the residue_summary_df suffix
            pattern = os.path.join(directory, f"{prefix}*.csv")
            files = glob.glob(pattern)
        
        all_files.extend(files)
    
    return sorted(all_files)


def load_and_validate_csv(csv_file):
    """
    Load CSV file and validate required columns.
    
    Parameters
    ----------
    csv_file : str
        Path to CSV file.
        
    Returns
    -------
    pd.DataFrame or None
        DataFrame with 'resid' and 'sum_ns' columns, or None if invalid.
    """
    try:
        df = pd.read_csv(csv_file)
        
        # Check for required columns
        if 'resid' not in df.columns:
            print(f"WARNING: Skipping {csv_file} - missing 'resid' column")
            print(f"  Available columns: {list(df.columns)}")
            return None
        
        if 'sum_ns' not in df.columns:
            print(f"WARNING: Skipping {csv_file} - missing 'sum_ns' column")
            print(f"  Available columns: {list(df.columns)}")
            return None
        
        # Return only the required columns
        return df[['resid', 'sum_ns']]
        
    except Exception as e:
        print(f"ERROR: Failed to read {csv_file}: {e}")
        return None


def sum_lifetimes(csv_files):
    """
    Sum lifetime values across multiple CSV files.
    
    Parameters
    ----------
    csv_files : list
        List of CSV file paths to sum.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with 'resid' and 'sum_ns' columns containing summed values.
    """
    all_data = []
    
    for csv_file in csv_files:
        basename = os.path.basename(csv_file)
        print(f"Reading: {basename}")
        
        df = load_and_validate_csv(csv_file)
        if df is not None:
            all_data.append(df)
            print(f"  Loaded {len(df)} residues")
    
    if not all_data:
        raise ValueError("No valid CSV files found with required columns")
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Group by resid and sum the sum_ns values
    summed_df = combined_df.groupby('resid', as_index=False)['sum_ns'].sum()
    
    # Sort by resid for cleaner output
    summed_df = summed_df.sort_values('resid').reset_index(drop=True)
    
    return summed_df


def main():
    """
    Main function: parse arguments and sum lifetime data.
    """
    parser = argparse.ArgumentParser(
        description="Sum lifetime values across multiple CSV files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -d ./lifetimes -o ./summary.csv -prefix psbs_
  %(prog)s -d ./lifetimes -o ./chain4_sum.csv -prefix chain_4_
        """
    )
    
    parser.add_argument('-d', '--dir', required=True, nargs='+',
                        help='Directory or directories containing CSV files')
    parser.add_argument('-o', '--output', required=True,
                        help='Output CSV file path')
    parser.add_argument('-prefix', '--prefix', required=True,
                        help='Prefix of CSV files to include (e.g., psbs_)')
    
    args = parser.parse_args()
    
    # Validate input directories
    directories = args.dir if isinstance(args.dir, list) else [args.dir]
    valid_dirs = []
    for directory in directories:
        if os.path.isdir(directory):
            valid_dirs.append(directory)
        else:
            print(f"WARNING: Directory not found: {directory}")
    
    if not valid_dirs:
        raise FileNotFoundError(f"No valid directories found in: {args.dir}")
    
    print(f"Input directories: {valid_dirs}")
    print(f"File prefix: {args.prefix}")
    print(f"Output file: {args.output}")
    print()
    
    # Find matching CSV files
    csv_files = find_csv_files(valid_dirs, args.prefix)
    
    if not csv_files:
        raise FileNotFoundError(
            f"No CSV files found matching pattern: {args.prefix}*.csv in {args.dir}"
        )
    
    print(f"Found {len(csv_files)} CSV files matching prefix '{args.prefix}'")
    print()
    
    # Sum lifetime data
    summed_df = sum_lifetimes(csv_files)
    
    print()
    print(f"Total unique residues: {len(summed_df)}")
    print(f"Sum of all lifetimes: {summed_df['sum_ns'].sum():.2f} ns")
    print()
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # Save to CSV
    summed_df.to_csv(args.output, index=False)
    print(f"Saved summed lifetimes to: {args.output}")
    
    # Show first few rows
    print("\nFirst 10 rows:")
    print(summed_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
