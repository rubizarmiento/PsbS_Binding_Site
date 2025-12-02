import os
import pandas as pd
from pathlib import Path 
import yaml

HELIX_DEFINITIONS_YAML_GROUP1="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_helix.yaml"
HELIX_DEFINITIONS_YAML_GROUP2="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix.yaml"
BASENAMES_CSV="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/4_trj_cluster/binding_basenames.csv"
OUTPUT = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/10_database/database.csv"
CSV_FILES = '/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/7_lifetimes_grouped/psbs*events*.csv'
ADD_LABELS_VALUES = ["pairs"]
ADD_LABELS_COLNAME = ["sim_type"]

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
    return sorted(Path('.').glob(glob_pattern))

def concat_dfs(dfs):
    """Concatenate list of dataframes."""
    return pd.concat(dfs, ignore_index=True)

def stats_to_dfs(df, group_by, values, stats=["mean", "std", "median", "min", "max", "count"]):
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
    'CYS': 'special'
}

def add_restypes_to_dfs(df, resname_col1='resname1', resname_col2='resname2'):
    """Add restype columns using vectorized .map()."""
    restype_dict = get_restype_dict()
    
    # VECTORIZED: Apply to entire column at once
    df['restype1'] = df[resname_col1].map(restype_dict).fillna('unknown')
    df['restype2'] = df[resname_col2].map(restype_dict).fillna('unknown')
    
    return df

def yaml_to_dict(yaml_file):
    """Load YAML file to dictionary."""
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def add_helix_to_dfs(df, helix_dict, resid_col='resid1', helix_col='helix1'):
    """Add helix information using vectorized operations."""
    # Create mapping dataframe from helix dict
    helix_list = []
    for helix_name, residues in helix_dict.items():
        for resid in residues:
            helix_list.append({'resid': resid, 'helix': helix_name})
    
    helix_df = pd.DataFrame(helix_list)
    
    # VECTORIZED: Merge instead of looping
    df = df.merge(helix_df, left_on=resid_col, right_on='resid', how='left')
    df = df.rename(columns={'helix': helix_col})
    
    return df

def write_df(df, output):
    """Write dataframe to CSV file."""
    df.to_csv(output, index=False)
    print(f"Saved to {output}")

def add_labels_dfs(df, colname, values):
    # Check if lengths match
    if len(colname) != len(values):
        print("ERROR: Lengths of colname and values do not match")
        exit(1)
    
    for i, val in enumerate(values):
        df.at[i, colname] = val
    return df

def main():
    # Check if files exist
    input_files = [HELIX_DEFINITIONS_YAML_GROUP1, HELIX_DEFINITIONS_YAML_GROUP2, BASENAMES_CSV]
    files_exist(input_files)

    # Get array of files 
    csv_files_arr = get_files_in_dir(CSV_FILES)
    
    if not csv_files_arr:
        print(f"ERROR: No files found matching {CSV_FILES}")
        exit(1)

    # Convert to dataframes
    dfs = csv_to_dfs(csv_files_arr)

    # Drop empty dataframes
    dfs = [df for df in dfs if not df.empty]

    # Concatenate dataframes
    concat_df = concat_dfs(dfs)
    
    if concat_df is None:
        print("ERROR: No data to process")
        exit(1)

    # Add lifetime stats
    concat_df = stats_to_dfs(concat_df, 
                            group_by=["resid1", "resid2", "chain1", "chain2"],
                            values="lifetimes",
                            stats=["mean", "std", "median", "min", "max", "count"])

    # Add restype columns
    concat_df = add_restypes_to_dfs(concat_df)

    # Get helix dictionaries
    dict_helix_group1 = yaml_to_dict(HELIX_DEFINITIONS_YAML_GROUP1)
    dict_helix_group2 = yaml_to_dict(HELIX_DEFINITIONS_YAML_GROUP2)

    # Add helix information
    concat_df = add_helix_to_dfs(concat_df, dict_helix_group1, resid_col='resid1', helix_col='helix1')
    concat_df = add_helix_to_dfs(concat_df, dict_helix_group2, resid_col='resid2', helix_col='helix2')

    # Add additional labels
    concat_df = add_labels_dfs(concat_df, ADD_LABELS_COLNAME, ADD_LABELS_VALUES)

    # Write result
    write_df(concat_df, OUTPUT)

if __name__ == '__main__':
    main()