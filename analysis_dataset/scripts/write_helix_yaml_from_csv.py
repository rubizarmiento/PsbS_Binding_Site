"""
Writes a Helix Definitions yaml file from a csv file and a pdb file.
A equivalent_chains.yaml can be used to map chain IDs in the PDB file to the chain IDs used in the CSV file.

yaml format:
  chain_ID:
    helix_name:
      start: residue_number  # int
      end: residue_number    # int
      start_res_name: "RES"  # str, three-letter residue name
      end_res_name: "RES"    # str
      length: helix_length   # int
      class: helix_class     # int (typically 1 for alpha helices)

csv format:
chain_ID, helix_name, class, start, end

equivalent_chains.yaml format:
  chain_ID:
    - "chain_id1"
    - "chain_id2"
  ```
e.g.
  A:
    - "A"
    - "B"

Arguments:
-f --pdb: Path to the PDB file.
-c --csv: Path to the CSV file with helix definitions.
-o --output: Path to the output YAML file (e.g., helix_definitions.yaml
-yaml --equivalent_chains: Optional path to a YAML file containing equivalent chain mappings.

Example usage:
python3 write_helix_yaml_from_csv.py \
  --pdb path/to/structure.pdb \
  --csv path/to/helix_definitions.csv \
  --output path/to/helix_definitions.yaml
  

"""
import os
import argparse
import yaml
import pandas as pd
import MDAnalysis as mda

def parse_arguments():
    parser = argparse.ArgumentParser(description="Write helix definitions YAML from CSV and PDB files.")
    parser.add_argument("-f", "--pdb", required=True, help="Path to the PDB file.")
    parser.add_argument("-c", "--csv", required=True, help="Path to the CSV file with helix definitions.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output YAML file (e.g., helix_definitions.yaml).")
    parser.add_argument("-y", "--equivalent_chains", help="Optional path to a YAML file containing equivalent chain mappings.")
    return parser.parse_args()

def read_helix_definitions(csv_file):
    """
    Read helix definitions from a CSV file.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing helix definitions.

    Returns
    -------
    list of dict
        List of helix definitions, where each definition is a dictionary with keys:
        'chain_ID', 'helix_name', 'class', 'start', 'end'.
    """
    df = pd.read_csv(csv_file)
    helix_definitions = df.to_dict(orient='records')
    print(f"Read {len(helix_definitions)} helix definitions from {csv_file}")
    return helix_definitions

def get_universe(pdb_file):
    """
    Load a PDB file into an MDAnalysis Universe.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    -------
    MDAnalysis.Universe
        The loaded Universe object.
    """
    try:
        u = mda.Universe(pdb_file)
        print(f"Successfully loaded PDB file: {pdb_file}")
        return u
    except Exception as e:
        print(f"Error loading PDB file {pdb_file}: {e}")
        raise

def read_equivalent_chains(equivalent_chains_file):
    """
    Read equivalent chains from a YAML file.

    Parameters
    ----------
    equivalent_chains_file : str
        Path to the YAML file containing equivalent chain mappings.

    Returns
    -------
    dict
        A dictionary mapping chain IDs in the PDB file to lists of equivalent chain IDs.
    """
    if equivalent_chains_file is None:
        return {}

    if not os.path.exists(equivalent_chains_file):
        print(f"Equivalent chains file not found: {equivalent_chains_file}")
        return {}
    
    with open(equivalent_chains_file, 'r') as f:
        equivalent_chains = yaml.safe_load(f)
    
    print(f"Read equivalent chains from {equivalent_chains_file}: {equivalent_chains}")
    return equivalent_chains

def add_resnames_and_chain_lenghts_to_df(u,df):
    """
    Add start_res_name, end_res_name, and helix_length columns to the DataFrame based on the PDB structure.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The loaded Universe object containing the structure.
    df : pandas.DataFrame
        DataFrame with columns 'chain_ID', 'helix_name', 'class', 'start', 'end'.

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with additional columns for residue names and helix lengths.
    """
    start_res_names = []
    end_res_names = []
    helix_lengths = []

    for index, row in df.iterrows():
        chain_id = row['chain_ID']
        start_resid = row['start']
        end_resid = row['end']

        # Select residues in the specified chain and range
        selection_str = f"chainID {chain_id} and resid {start_resid}:{end_resid}"
        selected_residues = u.select_atoms(selection_str).residues

        if len(selected_residues) == 0:
            print(f"Warning: No residues found for selection: {selection_str}")
            start_res_names.append("UNK")
            end_res_names.append("UNK")
            helix_lengths.append(0)
            continue

        start_res_names.append(selected_residues[0].resname)
        end_res_names.append(selected_residues[-1].resname)
        helix_lengths.append(len(selected_residues))

    df['start_res_name'] = start_res_names
    df['end_res_name'] = end_res_names
    df['length'] = helix_lengths

    return df

def duplicate_df_for_equivalent_chains(df, equivalent_chains):
    """
    Duplicate helix definitions for equivalent chains.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with helix definitions, including a 'chain_ID' column.
    equivalent_chains : dict
        A dictionary mapping chain IDs to lists of equivalent chain IDs.

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with duplicated rows for equivalent chains.
    """
    duplicated_rows = []

    for index, row in df.iterrows():
        chain_id = row['chain_ID']
        if chain_id in equivalent_chains:
            for eq_chain_id in equivalent_chains[chain_id]:
                new_row = row.copy()
                new_row['chain_ID'] = eq_chain_id
                duplicated_rows.append(new_row)

    if duplicated_rows:
        dup_df = pd.DataFrame(duplicated_rows)
        df = pd.concat([df, dup_df], ignore_index=True)

    return df

def write_yaml(df, output_file):
    """
    Write the helix definitions to a YAML file organized by chain and helix name.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with helix definitions, including 'chain_ID' and 'helix_name' columns.
    output_file : str
        Path to the output YAML file.
    """
    yaml_data = {}

    for index, row in df.iterrows():
        chain_id = row['chain_ID']
        helix_name = row['helix_name']

        if chain_id not in yaml_data:
            yaml_data[chain_id] = {}

        yaml_data[chain_id][helix_name] = {
            'start': int(row['start']),
            'end': int(row['end']),
            'start_res_name': row['start_res_name'],
            'end_res_name': row['end_res_name'],
            'length': int(row['length']),
            'class': int(row['class'])
        }

    with open(output_file, 'w') as f:
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)

    print(f"Helix definitions written to {output_file}")

def main():
    args = parse_arguments()
    
    # Read helix definitions from CSV
    helix_definitions = read_helix_definitions(args.csv)
    
    # Load PDB structure
    universe = get_universe(args.pdb)
    
    # Read equivalent chains if available and duplicate rows accordingly
    equivalent_chains_file = args.equivalent_chains if args.equivalent_chains else None
    equivalent_chains = read_equivalent_chains(equivalent_chains_file)

    # Duplicate rows for equivalent chains BEFORE looking up residue names
    df = pd.DataFrame(helix_definitions)
    df = duplicate_df_for_equivalent_chains(df, equivalent_chains)
    
    # Ensure start/end are integers after concat in duplicate step
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    
    # Add residue names and helix lengths to the DataFrame (now covers all chains)
    df = add_resnames_and_chain_lenghts_to_df(universe, df)
    
    # Write the YAML file with helix definitions
    write_yaml(df, args.output)

if __name__ == "__main__":
    main()


    
