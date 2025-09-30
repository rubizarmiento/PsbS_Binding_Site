"""
Workflow:
1) Reads /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/clustering/members.csv with the format (ignoring the header):

cluster_id	member_indices
0	2	1	7	3
1	6	5	4
2	0

2) Per cluster id searches the trajectories with the path: 
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/chain{chain}_{member_id}*.xtc

3) Concatenates the trajectories with MDAnalsysis and saves them in the directory:
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/clustering/cluster_{cluster_id}.xtc

4) Per cluster id searches the cluster center with the path: 
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/chain_chain_{member_id}/clust_c075/centers.pdb

5) Concatenated the first frame of each PDB file in the directory: 
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/clustering/cluster_{cluster_id}_center.pdb

"""

import os
from sys import argv
# Fix OMP_NUM_THREADS to avoid numexpr error (disable parallelism as requested)
os.environ['OMP_NUM_THREADS'] = '1'

import glob
import MDAnalysis as mda
import pandas as pd
#Ignore warnings
import warnings
warnings.filterwarnings("ignore")

# Set the chain
chain = argv[1]

# Base paths
base_path = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis"
members_csv = f"{base_path}/binding_poses/chain_{chain}/clustering/members.csv"
clustering_dir = f"{base_path}/binding_poses/chain_{chain}/clustering"
top_lifetimes_dir = f"{base_path}/top_lifetimes_pdbs"
topology = f"{base_path}/chain_{chain}/initial_fit_merged.pdb"
binding_poses_dir = f"{base_path}/binding_poses/chain_{chain}"

# Ensure clustering directory exists
os.makedirs(clustering_dir, exist_ok=True)

# Read members.csv
with open(members_csv, 'r') as f:
    lines = f.readlines()

# Parse header and data
header = lines[0].strip().split('\t')
data = []
for line in lines[1:]:
    parts = line.strip().split('\t')
    cluster_id = int(parts[0])
    member_indices = list(map(int, parts[1:]))
    data.append({'cluster_id': cluster_id, 'member_indices': member_indices})

df = pd.DataFrame(data)

# Process each cluster
for idx, row in df.iterrows():
    cluster_id = row['cluster_id']
    member_indices = row['member_indices']
    
    # Step 2 & 3: Concatenate trajectories
    traj_files = []
    for member_id in member_indices:
        pattern = f"{top_lifetimes_dir}/chain_{chain}_{member_id}*.xtc"
        traj_files.extend(glob.glob(pattern))
    
    if traj_files:
        # Load trajectories (assuming same topology)
        universes = []
        for traj in traj_files:
            u = mda.Universe(topology, traj)
            universes.append(u)
        
        # Concatenate by writing sequentially
        output_xtc = f"{clustering_dir}/cluster_{cluster_id}.xtc"
        with mda.Writer(output_xtc, universes[0].atoms.n_atoms) as w:
            for u in universes:
                for ts in u.trajectory:
                    w.write(u.atoms)
        # Copy the topology file
        os.system(f"cp {topology} {clustering_dir}/cluster_{cluster_id}_topology.pdb")
        print(f"Saved concatenated trajectory for cluster {cluster_id} to {output_xtc}")
    
    # Step 4 & 5: Concatenate PDB centers (as multi-frame PDB)
    pdb_files = []
    for member_id in member_indices:
        pdb_path = f"{binding_poses_dir}/chain_{chain}_{member_id}/clust_c075/centers.pdb"
        if os.path.exists(pdb_path):
            pdb_files.append(pdb_path)
    
    if pdb_files:
        # Load PDBs and concatenate as multi-frame
        universes_pdb = [mda.Universe(pdb) for pdb in pdb_files]
        output_pdb = f"{clustering_dir}/cluster_{cluster_id}_center.pdb"
        with mda.Writer(output_pdb, universes_pdb[0].atoms.n_atoms) as w:
            for u in universes_pdb:
                w.write(u.atoms)
        print(f"Saved concatenated PDB centers for cluster {cluster_id} to {output_pdb}")

print("Clustering trajectories and PDBs completed.")