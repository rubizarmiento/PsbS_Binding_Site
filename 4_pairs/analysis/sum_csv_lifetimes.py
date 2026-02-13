#!/usr/bin/env python3
"""
Sum lifetime values across multiple CSV files.

This script reads multiple CSV files with a common prefix from a directory,
extracts the 'resid' and 'sum_ns' columns, and sums the 'sum_ns' values 
for each residue across all files. If available, 'resname_i' and 'chainID_i'
are renamed to 'resname' and 'chain' in the output.

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
- resname : str (optional) - Residue name (from resname_i)
- chain : str (optional) - Chain ID (from chainID_i)

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


def find_csv_files(directories, prefix, ignore_missing=False):
    """
    Find all CSV files matching the prefix pattern across multiple directories.
    
    Parameters
    ----------
    directories : str or list
        Directory or list of directories to search for CSV files.
    prefix : str
        Prefix pattern to match (e.g., 'psbs_').
    ignore_missing : bool, optional
        If True, skip missing directories with warning. If False, raise error.
        
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
        DataFrame with 'resid', 'sum_ns', and optionally 'resname_i', 'chainID_i' columns, or None if invalid.
    """
    try:
        # Read CSV with keep_default_na=False to treat "NA" as a string literal, not NaN
        df = pd.read_csv(csv_file, keep_default_na=False, na_values=[''])
        
        # Check for required columns
        if 'resid' not in df.columns:
            print(f"WARNING: Skipping {csv_file} - missing 'resid' column")
            print(f"  Available columns: {list(df.columns)}")
            return None
        
        if 'sum_ns' not in df.columns:
            print(f"WARNING: Skipping {csv_file} - missing 'sum_ns' column")
            print(f"  Available columns: {list(df.columns)}")
            return None
        
        # Start with required columns
        cols_to_keep = ['resid', 'sum_ns']
        
        # Add optional columns if they exist
        if 'resname_i' in df.columns:
            cols_to_keep.append('resname_i')
        if 'chainID_i' in df.columns:
            cols_to_keep.append('chainID_i')
        
        # Return the selected columns
        return df[cols_to_keep]
        
    except Exception as e:
        print(f"ERROR: Failed to read {csv_file}: {e}")
        return None


def sum_lifetimes(csv_files, ignore_missing=False):
    """
    Sum lifetime values across multiple CSV files.
    
    Parameters
    ----------
    csv_files : list
        List of CSV file paths to sum.
    ignore_missing : bool, optional
        If True, return empty dataframe if no valid files found.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with 'resid', 'sum_ns', and optionally 'resname', 'chain' columns.
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
        if ignore_missing:
            print("WARNING: No valid CSV files found with required columns")
            return pd.DataFrame(columns=['resid', 'sum_ns'])
        else:
            raise ValueError("No valid CSV files found with required columns")
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Check if we have resname_i and chainID_i columns
    has_resname = 'resname_i' in combined_df.columns
    has_chain = 'chainID_i' in combined_df.columns
    
    # Group by resid (and resname_i, chainID_i if present) and sum the sum_ns values
    group_cols = ['resid']
    if has_resname:
        group_cols.append('resname_i')
    if has_chain:
        group_cols.append('chainID_i')
    
    # Use dropna=False to keep NaN values in groupby (important when resname_i is NaN for chain-level grouping)
    summed_df = combined_df.groupby(group_cols, as_index=False, dropna=False)['sum_ns'].sum()
    
    # Rename columns: resname_i -> resname, chainID_i -> chain
    rename_map = {}
    if has_resname:
        rename_map['resname_i'] = 'resname'
    if has_chain:
        rename_map['chainID_i'] = 'chain'
    
    if rename_map:
        summed_df = summed_df.rename(columns=rename_map)
    
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
    parser.add_argument('--ignore-missing', action='store_true',
                        help='Skip missing directories/files with warning instead of raising error')
    
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
        if args.ignore_missing:
            print("WARNING: No valid directories found, but continuing due to --ignore-missing flag")
        else:
            raise FileNotFoundError(f"No valid directories found in: {args.dir}")
    
    print(f"Input directories: {valid_dirs}")
    print(f"File prefix: {args.prefix}")
    print(f"Output file: {args.output}")
    print()
    
    # Find matching CSV files
    csv_files = find_csv_files(valid_dirs, args.prefix, ignore_missing=args.ignore_missing)
    
    if not csv_files:
        if args.ignore_missing:
            print(f"WARNING: No CSV files found matching pattern: {args.prefix}*.csv in {args.dir}")
            print("Creating empty output file due to --ignore-missing flag")
            # Create empty dataframe with expected columns
            empty_df = pd.DataFrame(columns=['resid', 'sum_ns'])
            empty_df.to_csv(args.output, index=False)
            print(f"Saved empty file to: {args.output}")
            return
        else:
            raise FileNotFoundError(
                f"No CSV files found matching pattern: {args.prefix}*.csv in {args.dir}"
            )
    
    print(f"Found {len(csv_files)} CSV files matching prefix '{args.prefix}'")
    print()
    
    # Sum lifetime data
    summed_df = sum_lifetimes(csv_files, ignore_missing=args.ignore_missing)
    
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
