"""
Workflow:

1) Load dataframes with pandas in the path:
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/lifetimes_chain_{chain}/lifetime_analysis_chain_{chain}_{id}_resname_residue_summary_df.csv
2) Sort by the columns "median_ns"
3) If one of the top 3 resnames is CLA or CLB they are classified as chl. If not as non-chl.
4) Inside the path: 
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs
   we looked for the trajectories with the name chain_{chain}_{id}.xtc and joined with MDAnalysis them as:
     /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/chain_{chain}_chl.xtc and 
     /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/chain_{chain}_non_chl.xtc
   according to the category they belong.
5) The name of the chl trajectories and non_chl trajectories is written in 
     /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/chain_{chain}_chl.txt and 
     /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/chain_{chain}_non_chl.txt
"""
import os
# Fix OMP_NUM_THREADS to avoid numexpr error (disable parallelism as requested)
import os
# Fix OMP_NUM_THREADS to avoid numexpr error (disable parallelism as requested)
os.environ['OMP_NUM_THREADS'] = '1'

import glob
import pandas as pd
import MDAnalysis as mda
import shutil
import warnings
warnings.filterwarnings("ignore")

def concatenate_trajectories(traj_files, output_path, topology):
    if not traj_files:
        print(f"No trajectories to concatenate for {output_path}")
        return
    if not os.path.exists(topology):
        print(f"Topology not found: {topology}")
        return
    universes = []
    for traj in traj_files:
        u = mda.Universe(topology, traj)
        universes.append(u)
    
    if universes:
        with mda.Writer(output_path, universes[0].atoms.n_atoms) as w:
            for u in universes:
                for ts in u.trajectory:
                    w.write(u.atoms)
        print(f"Concatenated trajectories into {output_path}")
        
        # Copy topology
        topology_dst = output_path.replace('.xtc', '.pdb')
        shutil.copy(topology, topology_dst)
        print(f"Copied topology to {topology_dst}")

# Set the chain
chains = ["4","r","s"]
for chain in chains:
    print(f"Processing chain {chain}")
    # Base paths
    base_path = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs"
    joined_dir = f"{base_path}/joined"
    os.makedirs(joined_dir, exist_ok=True)
    lifetimes_dir = f"{base_path}/lifetimes_chain_{chain}"
    topology = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/initial_fit_merged.pdb"
    # Find all CSV files
    csv_pattern = f"{lifetimes_dir}/lifetime_analysis_chain_{chain}_*_resname_residue_summary_df.csv"
    csv_files = glob.glob(csv_pattern)
    print(f"Found {len(csv_files)} CSV files")

    chl_ids = []
    non_chl_ids = []

    for csv_file in csv_files:
        # Extract id from filename
        filename = os.path.basename(csv_file)
        # Assuming format: lifetime_analysis_chain_{chain}_{id}_resname_residue_summary_df.csv
        parts = filename.split('_')
        id_str = parts[4]  # e.g., '1' from chain_4_1_resname...
        id_num = int(id_str)
        
        # Load dataframe
        df = pd.read_csv(csv_file)
        
        # Sort by median_ns descending
        df_sorted = df.sort_values(by='median_ns', ascending=False)
        
        # Get top 3 resnames
        top_resnames = df_sorted['resid'].head(5).tolist()
        
        # Classify
        if 'CLA' in top_resnames or 'CLB' in top_resnames:
            chl_ids.append(id_num)
        else:
            non_chl_ids.append(id_num)

    print(f"Chl ids: {chl_ids}")
    print(f"Non-chl ids: {non_chl_ids}")

    # Chl trajectories
    chl_traj_files = []
    for id in chl_ids:
        pattern = f"{base_path}/chain_{chain}_{id}*.xtc"
        matches = glob.glob(pattern)
        if matches:
            chl_traj_files.extend(matches)
    print(f"Chl traj files: {chl_traj_files}")
    concatenate_trajectories(chl_traj_files, f"{joined_dir}/chain_{chain}_chl.xtc", topology)

    # Non-chl trajectories
    non_chl_traj_files = []
    for id in non_chl_ids:
        pattern = f"{base_path}/chain_{chain}_{id}*.xtc"
        matches = glob.glob(pattern)
        if matches:
            non_chl_traj_files.extend(matches)
    print(f"Non-chl traj files: {non_chl_traj_files}")
    concatenate_trajectories(non_chl_traj_files, f"{joined_dir}/chain_{chain}_non_chl.xtc", topology)

    # Write lists to txt files
    with open(f"{joined_dir}/chain_{chain}_chl.txt", 'w') as f:
        for id in chl_ids:
            f.write(f"chain_{chain}_{id}.xtc\n")

    with open(f"{joined_dir}/chain_{chain}_non_chl.txt", 'w') as f:
        for id in non_chl_ids:
            f.write(f"chain_{chain}_{id}.xtc\n")

    print("Classification and concatenation completed.")

import glob
import pandas as pd
import MDAnalysis as mda
import shutil
import warnings
warnings.filterwarnings("ignore")

# Set the chain
chains = ["4","r","s"]
for chain in chains:
    print(f"Processing chain {chain}")
    # Base paths
    base_path = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs"
    joined_dir = f"{base_path}/joined"
    os.makedirs(joined_dir, exist_ok=True)
    lifetimes_dir = f"{base_path}/lifetimes_chain_{chain}"
    topology = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/initial_fit_merged.pdb"
    # Find all CSV files
    csv_pattern = f"{lifetimes_dir}/lifetime_analysis_chain_{chain}_*_resname_residue_summary_df.csv"
    csv_files = glob.glob(csv_pattern)
    print(f"Found {len(csv_files)} CSV files")

    chl_ids = []
    non_chl_ids = []

    for csv_file in csv_files:
        # Extract id from filename
        filename = os.path.basename(csv_file)
        # Assuming format: lifetime_analysis_chain_{chain}_{id}_resname_residue_summary_df.csv
        parts = filename.split('_')
        id_str = parts[4]  # e.g., '1' from chain_4_1_resname...
        id_num = int(id_str)
        
        # Load dataframe
        df = pd.read_csv(csv_file)
        
        # Sort by median_ns descending
        df_sorted = df.sort_values(by='median_ns', ascending=False)
        
        # Get top 3 resnames
        top_resnames = df_sorted['resid'].head(3).tolist()
        
        # Classify
        if 'CLA' in top_resnames or 'CLB' in top_resnames:
            chl_ids.append(id_num)
        else:
            non_chl_ids.append(id_num)

    print(f"Chl ids: {chl_ids}")
    print(f"Non-chl ids: {non_chl_ids}")

    # Now, concatenate trajectories
def concatenate_trajectories(traj_files, output_path, topology):
    if not traj_files:
        print(f"No trajectories to concatenate for {output_path}")
        return
    if not os.path.exists(topology):
        print(f"Topology not found: {topology}")
        return
    universes = []
    for traj in traj_files:
        u = mda.Universe(topology, traj)
        universes.append(u)
    
    if universes:
        with mda.Writer(output_path, universes[0].atoms.n_atoms) as w:
            for u in universes:
                for ts in u.trajectory:
                    w.write(u.atoms)
        print(f"Concatenated trajectories into {output_path}")
        
        # Copy topology
        topology_dst = output_path.replace('.xtc', '.pdb')
        shutil.copy(topology, topology_dst)
        print(f"Copied topology to {topology_dst}")    # Chl trajectories
    chl_traj_files = []
    for id in chl_ids:
        pattern = f"{base_path}/chain_{chain}_{id}*.xtc"
        matches = glob.glob(pattern)
        if matches:
            chl_traj_files.extend(matches)
    print(f"Chl traj files: {chl_traj_files}")
    concatenate_trajectories(chl_traj_files, f"{joined_dir}/chain_{chain}_chl.xtc", topology)

    # Non-chl trajectories
    non_chl_traj_files = []
    for id in non_chl_ids:
        pattern = f"{base_path}/chain_{chain}_{id}*.xtc"
        matches = glob.glob(pattern)
        if matches:
            non_chl_traj_files.extend(matches)
    print(f"Non-chl traj files: {non_chl_traj_files}")
    concatenate_trajectories(non_chl_traj_files, f"{joined_dir}/chain_{chain}_non_chl.xtc", topology)

    # Write lists to txt files
    with open(f"{joined_dir}/chain_{chain}_chl.txt", 'w') as f:
        for id in chl_ids:
            f.write(f"chain_{chain}_{id}.xtc\n")

    with open(f"{joined_dir}/chain_{chain}_non_chl.txt", 'w') as f:
        for id in non_chl_ids:
            f.write(f"chain_{chain}_{id}.xtc\n")

    print("Classification and concatenation completed.")