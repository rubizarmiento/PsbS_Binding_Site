"""
Database Creation Script for Protein Binding Events

Processes CSV files containing protein residue pair binding events,
aggregates lifetime statistics, and adds structural information
(helix assignments) and labels (chain identifiers).

The script reads multiple CSV files matching a glob pattern, concatenates them,
adds aggregated statistics (mean, median, etc.) for lifetime values grouped by
residue pairs, and maps residue IDs to helix names and chain labels using YAML
configuration files.

FILES REQUIREMENTS
-----
- Input CSV files must contain columns: resid_i, resid_j, chainID_i, chainID_j, lifetime_ns
- YAML helix definition structure: {chain: {helix_name: {start: int, end: int, class: str}}}
- YAML chain label structure: {chain_id: label_name}
- Output includes aggregated statistics: mean, std, median, min, max, count, sum

Author
------
Rubi Zarmiento-Garcia


Usage
-----
Example:

    python write_databases.py \\
        --csv_files "/path/to/data/*respairs_events*.csv" \\
        --output "/path/to/output/database.csv" \\
        --helix_def_yaml_group1 "/path/to/helix_labels_1.yaml" \\
        --helix_def_yaml_group2 "/path/to/helix_labels_2.yaml" \\
        --labels_chain_yaml_group1 "/path/to/chain_labels_1.yaml" \\
        --labels_chain_yaml_group2 "/path/to/chain_labels_2.yaml" \\
        --add_labels_colname "sim_type" \\
        --add_labels_values "pairs"

Required Arguments
------------------
--csv_files : str
    Glob pattern for input CSV files
--output : str
    Output CSV file path
--helix_def_yaml_group1 : str
    YAML file for helix definitions (group 1, typically binding proteins)
--helix_def_yaml_group2 : str
    YAML file for helix definitions (group 2, typically target protein)

Optional Arguments
------------------
--labels_chain_yaml_group1 : str
    YAML file for chain labels (group 1)
--labels_chain_yaml_group2 : str
    YAML file for chain labels (group 2)
--add_labels_colname : list of str
    Column names for additional custom labels
--add_labels_values : list of str
    Values for additional custom labels


"""

import os
import pandas as pd
from pathlib import Path 
import yaml
import glob
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process CSV files and write database.")
    parser.add_argument(
        "--csv_files",
        type=str,
        required=True,
        help="Glob pattern for input CSV files.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV file path."
    )
    parser.add_argument(
        "--add_labels_colname",
        type=str,
        nargs='+',
        required=False,
        help="Column names for additional labels."
    )
    parser.add_argument(
        "--add_labels_values",
        type=str,
        nargs='+',
        required=False,
        help="Values for additional labels."
    )
    parser.add_argument(
        "--helix_def_yaml_group1",
        type=str,
        required=True,
        help="YAML file for helix definitions group 1."
    )
    parser.add_argument(
        "--helix_def_yaml_group2",
        type=str,
        required=True,
        help="YAML file for helix definitions group 2."
    )
    parser.add_argument(
        "--labels_chain_yaml_group1",
        type=str,
        required=False,
        help="YAML file for chain labels group 1."
    )
    parser.add_argument(
        "--labels_chain_yaml_group2",
        type=str,
        required=False,
        help="YAML file for chain labels group 2."
    )
    return parser.parse_args()

def files_exist(files):
    """
    Check if files exist and exit if any are missing.

    Parameters
    ----------
    files : list of str
        List of file paths to check

    Raises
    ------
    SystemExit
        If any file does not exist
    """
    for file in files:
        if not os.path.isfile(file):
            print(f"ERROR: File {file} does not exist")
            exit(1)

def csv_to_dfs(files):
    """
    Convert CSV files to list of pandas DataFrames.

    Parameters
    ----------
    files : list of str
        List of CSV file paths to read

    Returns
    -------
    list of pd.DataFrame
        List of DataFrames, one for each CSV file
    """
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file, keep_default_na=False, na_values=[''])
        except:
            df = pd.read_csv(file, header=None, keep_default_na=False, na_values=[''])
        dfs.append(df)
    
    return dfs

def get_files_in_dir(glob_pattern):
    """
    Get files matching a glob pattern.

    Parameters
    ----------
    glob_pattern : str
        Glob pattern to match files (e.g., '/path/*.csv')

    Returns
    -------
    list of str
        Sorted list of file paths matching the pattern
    """
    return sorted(glob.glob(glob_pattern))

def concat_dfs(dfs):
    """
    Concatenate list of DataFrames into a single DataFrame.

    Parameters
    ----------
    dfs : list of pd.DataFrame
        List of DataFrames to concatenate

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame with reset index
    """
    return pd.concat(dfs, ignore_index=True)

def stats_to_dfs(df, group_by, values, stats=["mean", "std", "median", "min", "max", "count","sum"]):
    """
    Add aggregation statistics to DataFrame using vectorized transform.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    group_by : str or list of str
        Column name(s) to group by
    values : str
        Column name to calculate statistics on
    stats : list of str, optional
        List of statistics to calculate (default: ["mean", "std", "median", "min", "max", "count", "sum"])

    Returns
    -------
    pd.DataFrame
        DataFrame with added statistic columns named as '{stat}_{values}'
    
    Examples
    --------
    >>> df = stats_to_dfs(df, group_by=["resid_i", "resid_j"], values="lifetime_ns")
    """
    group = df.groupby(group_by)[values]
    
    # Vectorized: all stats at once
    for stat in stats:
        df[f'{stat}_{values}'] = group.transform(stat)
    
    return df


def yaml_to_dict(yaml_file):
    """
    Load YAML file into a Python dictionary.

    Parameters
    ----------
    yaml_file : str
        Path to YAML file

    Returns
    -------
    dict
        Dictionary containing YAML contents
    """
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def add_helix_to_dfs(df, helix_dict, resid_col='resid_i', helix_col='helix_i', chain_col='chainID_i'):
    """
    Add helix information to DataFrame using vectorized operations.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    helix_dict : dict
        Dictionary mapping chain names to helix definitions.
        Structure: {chain: {helix_name: {start: int, end: int, class: str}}}
    resid_col : str, optional
        Name of residue ID column (default: 'resid_i')
    helix_col : str, optional
        Name of output helix column (default: 'helix_i')
    chain_col : str, optional
        Name of chain ID column (default: 'chainID_i')

    Returns
    -------
    pd.DataFrame
        DataFrame with added helix information columns:
        - {helix_col} : helix name or 'loop' if not in a helix
        - {helix_col}_class : helix class or 'NA' if not applicable

    Notes
    -----
    Residues not mapped to any helix are assigned 'loop' as helix name
    and 'NA' as class. Chain names in YAML (e.g., 'chain_4') are mapped
    to single character chain IDs in the data (e.g., '4').
    """
    # Check if chain_col exists in the dataframe
    use_chain = chain_col is not None and chain_col in df.columns
    
    # Create mapping dataframe from helix dict
    helix_list = []
    
    if use_chain:
        unique_chains = df[chain_col].unique()
    
    # Parse YAML structure: chain -> helix_name -> helix_info
    for chain_name, helices in helix_dict.items():
        # Ensure chain_name is a string (YAML may load numeric keys as int)
        chain_name = str(chain_name)
        # Map chain names: if data has "4" but YAML has "chain_4", handle the mapping
        chain_name_mapped = chain_name.replace('chain_', '') if chain_name.startswith('chain_') else chain_name
        
        # Only process chains that are present in the data (if using chain filtering)
        if not use_chain or chain_name_mapped in [str(c) for c in unique_chains]:
            for helix_name, helix_info in helices.items():
                start_resid = helix_info['start']
                end_resid = helix_info['end']
                class_value = helix_info.get('class', 'None')  # Get class or 'None' if not found
                
                # Add all residues in this helix range
                for resid in range(start_resid, end_resid + 1):
                    if use_chain:
                        helix_list.append({'resid': resid, 'chain': chain_name_mapped, 'helix': helix_name, 'class': class_value})
                    else:
                        helix_list.append({'resid': resid, 'helix': helix_name, 'class': class_value})
    
    helix_df = pd.DataFrame(helix_list)
    
    if helix_df.empty:
        # If no helices found, create columns but with NaN values
        df[helix_col] = None
        df[f'{helix_col}_class'] = 'NA'
        return df
    
    # Merge based on whether we have chain information
    if use_chain:
        # Ensure chain column has the same type as the dataframe chain column
        helix_df['chain'] = helix_df['chain'].astype(df[chain_col].dtype)
        # VECTORIZED: Merge on both resid and chain to ensure correct matching
        df = df.merge(helix_df, left_on=[resid_col, chain_col], right_on=['resid', 'chain'], how='left')
        df = df.drop(['resid', 'chain'], axis=1)
    else:
        # Merge on resid only if no chain column available
        df = df.merge(helix_df, left_on=resid_col, right_on='resid', how='left')
        df = df.drop('resid', axis=1)
    
    df = df.rename(columns={'helix': helix_col, 'class': f'{helix_col}_class'})
    
    # Fill NaN values for residues not in helices
    # Set helix to 'loop' and class to 'NA' for unmapped residues
    df[helix_col] = df[helix_col].fillna('loop')
    df[f'{helix_col}_class'] = df[f'{helix_col}_class'].fillna('NA')
    
    return df

def add_col_if_match(df, dict_map, key_col, new_col, default_value='NA'):
    """
    Add a new column to DataFrame based on a mapping dictionary.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    dict_map : dict
        Mapping dictionary where keys match values in key_col
    key_col : str
        Name of column to use for mapping lookup
    new_col : str
        Name of new column to create
    default_value : str, optional
        Default value for unmapped entries (default: 'NA')

    Returns
    -------
    pd.DataFrame
        DataFrame with added column
    """
    df[new_col] = df[key_col].map(dict_map).fillna(default_value)
    return df

def write_df(df, output):
    """
    Write DataFrame to CSV file, creating parent directories if needed.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write
    output : str
        Output CSV file path

    Notes
    -----
    Creates parent directories automatically if they don't exist.
    Prints confirmation message upon successful write.
    """
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    print(f"Output file: {output}")

def add_labels_dfs(df, colname, values):
    """
    Add label columns to DataFrame with constant values.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    colname : list of str or None
        Column names for new label columns
    values : list of str or None
        Values to assign to each column (same value for all rows)

    Returns
    -------
    pd.DataFrame
        DataFrame with added label columns

    Raises
    ------
    SystemExit
        If lengths of colname and values don't match

    Notes
    -----
    Returns unmodified DataFrame if colname or values is None.
    Each column receives a constant value for all rows.
    """
    # Handle case where colname or values might be None
    if colname is None or values is None:
        return df
    
    # colname and values should be lists of the same length
    if len(colname) != len(values):
        print("ERROR: Lengths of colname and values do not match")
        exit(1)
    
    # Add each column with its corresponding value (same value for all rows)
    for col, val in zip(colname, values):
        df[col] = val
    return df

def main():
    # Parse arguments
    args = parse_arguments()
    
    # Check if files exist
    input_files = [args.helix_def_yaml_group1, args.helix_def_yaml_group2]
    files_exist(input_files)

    # Get array of files 
    csv_files_arr = get_files_in_dir(args.csv_files)
    
    if not csv_files_arr:
        print(f"ERROR: No files found matching {args.csv_files}")
        exit(1)

    # Convert to dataframes
    dfs = csv_to_dfs(csv_files_arr)

    # Add filenames as columns
    for i, df in enumerate(dfs):
        filename = os.path.basename(csv_files_arr[i])
        df['source_file'] = filename

    # Drop empty dataframes
    dfs = [df for df in dfs if not df.empty]

    # Concatenate dataframes
    concat_df = concat_dfs(dfs)
    
    if concat_df is None:
        print("ERROR: No data to process")
        exit(1)

    # Convert chainID columns to string
    concat_df['chainID_i'] = concat_df['chainID_i'].astype(str)
    concat_df['chainID_j'] = concat_df['chainID_j'].astype(str)

    # Add lifetime stats
    concat_df = stats_to_dfs(concat_df, 
                            group_by=["resid_i", "resid_j"],
                            values="lifetime_ns",
                            stats=["mean", "std", "median", "min", "max", "count", "sum"])

    # Get helix dictionaries
    dict_helix_group1 = yaml_to_dict(args.helix_def_yaml_group1)
    dict_helix_group2 = yaml_to_dict(args.helix_def_yaml_group2)

    # Add helix information
    concat_df = add_helix_to_dfs(concat_df, dict_helix_group1, resid_col='resid_i', helix_col='helix_i', chain_col='chainID_i')
    concat_df = add_helix_to_dfs(concat_df, dict_helix_group2, resid_col='resid_j', helix_col='helix_j', chain_col='chainID_j')

    # Get chain labels
    if args.labels_chain_yaml_group1:
        dict_chain_labels_group1 = yaml_to_dict(args.labels_chain_yaml_group1)
        concat_df = add_col_if_match(concat_df, dict_chain_labels_group1, key_col='chainID_i', new_col='chain_label_i', default_value='NA')
    
    if args.labels_chain_yaml_group2:
        dict_chain_labels_group2 = yaml_to_dict(args.labels_chain_yaml_group2)
        concat_df = add_col_if_match(concat_df, dict_chain_labels_group2, key_col='chainID_j', new_col='chain_label_j', default_value='NA')
    
    # Add additional labels (only if both colname and values are provided)
    if args.add_labels_colname and args.add_labels_values:
        concat_df = add_labels_dfs(concat_df, args.add_labels_colname, args.add_labels_values)

    # Write result
    write_df(concat_df, args.output)

if __name__ == '__main__':
    main()