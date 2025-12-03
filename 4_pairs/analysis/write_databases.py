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
            df = pd.read_csv(file)
        except:
            df = pd.read_csv(file, header=None)
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


def get_restype_dict():
    """Get residue type dictionary."""
    return {
    # Hydrophobic (aliphatic)
    'ALA': 'hydrophobic',
    'VAL': 'hydrophobic',
    'LEU': 'hydrophobic',
    'ILE': 'hydrophobic',
    'MET': 'hydrophobic',

    # Aromatic
    'PHE': 'aromatic',
    'TYR': 'aromatic',
    'TRP': 'aromatic',

    # Polar uncharged
    'SER': 'polar',
    'THR': 'polar',
    'ASN': 'polar',
    'GLN': 'polar',

    # Negative (acidic)
    'ASP': 'acidic',
    'GLU': 'acidic',

    # Positive (basic)
    'LYS': 'basic',
    'ARG': 'basic',
    'HIS': 'basic',

    # Special cases
    'GLY': 'special',
    'PRO': 'special',
    'CYS': 'special',

    # Chlorophylls
    'CLA': 'chlorophyll',
    'CLB': 'chlorophyll',
    'CHL': 'chlorophyll',

    # Carotenoids
    'LUT': 'carotenoid',
    'VIO': 'carotenoid',
    'XAT': 'carotenoid',
    'NEO': 'carotenoid',
    'BCR': 'carotenoid',
    'NEX': 'carotenoid'
}

def get_type_dict():
    """Get type protein or cofactor dictionary."""
    return {
    # Proteins
    'ALA': 'protein',
    'VAL': 'protein',
    'LEU': 'protein',
    'ILE': 'protein',
    'MET': 'protein',
    'PHE': 'protein',
    'TYR': 'protein',
    'TRP': 'protein',
    'SER': 'protein',
    'THR': 'protein',
    'ASN': 'protein',
    'GLN': 'protein',
    'ASP': 'protein',
    'GLU': 'protein',
    'LYS': 'protein',
    'ARG': 'protein',
    'HIS': 'protein',
    'GLY': 'protein',
    'PRO': 'protein',
    'CYS': 'protein',

    # Chlorophylls
    'CLA': 'chlorophyll',
    'CLB': 'chlorophyll',
    'CHL': 'chlorophyll',

    # Carotenoids
    'LUT': 'LUT',
    'VIO': 'VIO',
    'XAT': 'XAT',
    'NEO': 'NEO',
    'BCR': 'BCR',
    'NEX': 'NEO'

}

def add_restypes_to_dfs(df, resname_col1='resname_i', resname_col2='resname_j'):
    """Add restype columns using vectorized .map()."""
    restype_dict = get_restype_dict()
    
    # VECTORIZED: Apply to entire column at once
    df['restype_i'] = df[resname_col1].map(restype_dict).fillna('unknown')
    df['restype_j'] = df[resname_col2].map(restype_dict).fillna('unknown')

    return df

def add_types_to_dfs(df, resname_col1='resname_i', resname_col2='resname_j'):
    """Add type columns using vectorized .map()."""
    type_dict = get_type_dict()
    
    # VECTORIZED: Apply to entire column at once
    df['type_i'] = df[resname_col1].map(type_dict).fillna('unknown')
    df['type_j'] = df[resname_col2].map(type_dict).fillna('unknown')

    return df

def yaml_to_dict(yaml_file):
    """Load YAML file to dictionary."""
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def add_helix_to_dfs(df, helix_dict, resid_col='resid_i', helix_col='helix_i', type_col='type_i', chain_col='chainID_i'):
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
    
    # Fill NaN values based on residue type
    # For protein residues not in helices, set helix to 'loop'
    # For cofactor residues not in helices, set helix to 'NA'
    if type_col in df.columns:
        mask_protein = (df[helix_col].isna()) & (df[type_col] == 'protein')
        mask_cofactor = (df[helix_col].isna()) & (df[type_col] == 'cofactor')
        
        df.loc[mask_protein, helix_col] = 'loop'
        df.loc[mask_cofactor, helix_col] = 'NA'
    
    # For class column, set to 'NA' for residues not in helices
    df[f'{helix_col}_class'] = df[f'{helix_col}_class'].fillna('NA')
    
    return df

def add_col_if_match(df, dict_map, key_col, new_col, default_value='NA'):
    """Add a new column to the dataframe based on a mapping dictionary."""
    df[new_col] = df[key_col].map(dict_map).fillna(default_value)
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

    # Add resid+resname column
    concat_df['resid_resname_i'] = concat_df['resid_i'].astype(str) + concat_df['resname_i'].astype(str)
    concat_df['resid_resname_j'] = concat_df['resid_j'].astype(str) + concat_df['resname_j'].astype(str)

    #Add resid_resname_i+resid_resname_j column
    concat_df['resid_resname_pair'] = concat_df['resid_resname_i'] + "-" + concat_df['resid_resname_j']

    # Add lifetime stats
    concat_df = stats_to_dfs(concat_df, 
                            group_by=["resid_i", "resid_j"],
                            values="lifetime_ns",
                            stats=["mean", "std", "median", "min", "max", "count", "sum"])

    # Add restype columns (commented out - no resname columns in CSV)
    concat_df = add_restypes_to_dfs(concat_df)

    # Add types (protein/cofactor) to dataframe
    concat_df = add_types_to_dfs(concat_df)

    # Get helix dictionaries
    dict_helix_group1 = yaml_to_dict(args.helix_def_yaml_group1)
    dict_helix_group2 = yaml_to_dict(args.helix_def_yaml_group2)

    # Add helix information
    concat_df = add_helix_to_dfs(concat_df, dict_helix_group1, resid_col='resid_i', helix_col='helix_i', type_col='type_i', chain_col='chainID_i')
    concat_df = add_helix_to_dfs(concat_df, dict_helix_group2, resid_col='resid_j', helix_col='helix_j', type_col='type_j', chain_col='chainID_j')

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