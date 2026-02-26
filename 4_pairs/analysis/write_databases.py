import os
import pandas as pd
from pathlib import Path 
import yaml
import glob
import argparse

def parse_arguments():
    """Parse command-line arguments."""
    # Default values    
    HELIX_DEFINITIONS_YAML_GROUP1="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix_labels.yaml"
    HELIX_DEFINITIONS_YAML_GROUP2="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_helix_labels_merged.yaml"
    LABELS_CHAIN_YAML_GROUP1="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_labels.yaml"
    LABELS_CHAIN_YAML_GROUP2="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_labels.yaml"
    OUTPUT = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/10_database/database.csv"
    CSV_FILES = '/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/analysis_pairs/chain_4/7_lifetimes_grouped/*respairs_events*.csv'
    ADD_LABELS_VALUES = ["pairs"]
    ADD_LABELS_COLNAME = ["sim_type"]

    parser = argparse.ArgumentParser(description="Process CSV files and write database.")
    parser.add_argument(
        "--csv_files",
        type=str,
        default=CSV_FILES,
        help="Glob pattern for input CSV files.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=OUTPUT,
        help="Output CSV file path."
    )
    parser.add_argument(
        "--add_labels_colname",
        type=str,
        nargs='+',
        default=ADD_LABELS_COLNAME,
        help="Column names for additional labels."
    )
    parser.add_argument(
        "--add_labels_values",
        type=str,
        nargs='+',
        default=ADD_LABELS_VALUES,
        help="Values for additional labels."
    )
    parser.add_argument(
        "--helix_def_yaml_group1",
        type=str,
        default=HELIX_DEFINITIONS_YAML_GROUP1,
        help="YAML file for helix definitions group 1."
    )
    parser.add_argument(
        "--helix_def_yaml_group2",
        type=str,
        default=HELIX_DEFINITIONS_YAML_GROUP2,
        help="YAML file for helix definitions group 2."
    )
    parser.add_argument(
        "--labels_chain_yaml_group1",
        type=str,
        default=LABELS_CHAIN_YAML_GROUP1,
        help="YAML file for chain labels group 1.",
        required=False
    )
    parser.add_argument(
        "--labels_chain_yaml_group2",
        type=str,
        default=LABELS_CHAIN_YAML_GROUP2,
        help="YAML file for chain labels group 2.",
        required=False  
    )
    return parser.parse_args()

def files_exist(files):
    """Check if files exist."""
    for file in files:
        if not os.path.isfile(file):
            print(f"ERROR: File {file} does not exist")
            exit(1)

def csv_to_dfs(files):
    """Convert CSV files to list of dataframes."""
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file, keep_default_na=False, na_values=[''])
        except:
            df = pd.read_csv(file, header=None, keep_default_na=False, na_values=[''])
        dfs.append(df)
    
    return dfs

def get_files_in_dir(glob_pattern):
    """Get files matching glob pattern."""
    return sorted(glob.glob(glob_pattern))

def concat_dfs(dfs):
    """Concatenate list of dataframes."""
    return pd.concat(dfs, ignore_index=True)

def stats_to_dfs(df, group_by, values, stats=["mean", "std", "median", "min", "max", "count","sum"]):
    """Add aggregation statistics using vectorized transform."""
    group = df.groupby(group_by)[values]
    
    # Vectorized: all stats at once
    for stat in stats:
        df[f'{stat}_{values}'] = group.transform(stat)
    
    return df


def yaml_to_dict(yaml_file):
    """Load YAML file to dictionary."""
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def add_helix_to_dfs(df, helix_dict, resid_col='resid_i', helix_col='helix_i', chain_col='chainID_i'):
    """Add helix information using vectorized operations."""
    # Check if chain_col exists in the dataframe
    use_chain = chain_col is not None and chain_col in df.columns
    
    # Create mapping dataframe from helix dict
    helix_list = []
    
    if use_chain:
        unique_chains = df[chain_col].unique()
    
    # Parse YAML structure: chain -> helix_name -> helix_info
    for chain_name, helices in helix_dict.items():
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
    """Add a new column to the dataframe based on a mapping dictionary."""
    df[new_col] = df[key_col].map(dict_map).fillna(default_value)
    return df

def apply_modifier_operation(df, modifier_config):
    """
    Apply conditional operations to dataframe columns based on modifier config.
    
    Parameters:
    -----------
    df : pd.DataFrame
        The dataframe to modify
    modifier_config : dict
        Configuration dictionary with keys:
        - 'col': source column name
        - 'new_col': target column name
        - 'operation': numeric value to apply (e.g., -212 means subtract 212)
        - 'condition': condition string (e.g., '> 212')
        - 'value_if_true': (optional) value to assign when condition is true
        - 'value_if_false': (optional) value to assign when condition is false
        
        If value_if_true/value_if_false are provided, performs conditional assignment.
        Otherwise, performs arithmetic operation.
    
    Returns:
    --------
    pd.DataFrame
        Modified dataframe with new column
    """
    col = modifier_config['col']
    new_col = modifier_config['new_col']
    operation = modifier_config.get('operation', 0)
    condition = modifier_config['condition']
    
    # Check if this is a conditional value assignment
    has_conditional_values = 'value_if_true' in modifier_config and 'value_if_false' in modifier_config
    
    # Parse condition (e.g., '> 212' -> operator='>' and value=212)
    condition_parts = condition.strip().split()
    if len(condition_parts) != 2:
        raise ValueError(f"Invalid condition format: {condition}. Expected format: '> 212'")
    
    operator = condition_parts[0]
    threshold = float(condition_parts[1])
    
    # Apply operation based on condition
    if operator == '>':
        mask = df[col] > threshold
    elif operator == '>=':
        mask = df[col] >= threshold
    elif operator == '<':
        mask = df[col] < threshold
    elif operator == '<=':
        mask = df[col] <= threshold
    elif operator == '==':
        mask = df[col] == threshold
    elif operator == '!=':
        mask = df[col] != threshold
    else:
        raise ValueError(f"Unsupported operator: {operator}")
    
    if has_conditional_values:
        # Conditional value assignment mode
        value_if_true = modifier_config['value_if_true']
        value_if_false = modifier_config['value_if_false']
        
        # Initialize new column with proper dtype to avoid warnings
        # Use object dtype to handle any type (string, int, float)
        df[new_col] = None
        df[new_col] = df[new_col].astype(object)
        df.loc[mask, new_col] = value_if_true
        df.loc[~mask, new_col] = value_if_false
    else:
        # Arithmetic operation mode
        # Create new column, starting with original values
        df[new_col] = df[col]
        # Apply the operation (add the operation value, which can be negative)
        df.loc[mask, new_col] = df.loc[mask, col] + operation
    
    return df

def write_df(df, output):
    """Write dataframe to CSV file."""
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    print(f"Saved to {output}")

def add_labels_dfs(df, colname, values):
    """Add label columns to dataframe."""
    # colname and values should be lists of the same length
    if len(colname) != len(values):
        print("ERROR: Lengths of colname and values do not match")
        exit(1)
    
    # Add each column with its corresponding value (same value for all rows)
    for col, val in zip(colname, values):
        df[col] = val
    return df

def main():
    # Parse command-line arguments
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
    if args.labels_chain_yaml_group1:
        concat_df = add_col_if_match(concat_df, dict_chain_labels_group1, key_col='chainID_i', new_col='chain_label_i', default_value='NA')
    if args.labels_chain_yaml_group2:
        dict_chain_labels_group2 = yaml_to_dict(args.labels_chain_yaml_group2)
        concat_df = add_col_if_match(concat_df, dict_chain_labels_group2, key_col='chainID_j', new_col='chain_label_j', default_value='NA')
    

    # Add additional labels
    concat_df = add_labels_dfs(concat_df, args.add_labels_colname, args.add_labels_values)

    # Write result
    write_df(concat_df, args.output)

if __name__ == '__main__':
    main()