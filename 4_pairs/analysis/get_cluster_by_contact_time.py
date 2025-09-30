"""
Workflow:
Read the data frame:
df= ${an1}/csv_files/lifetime/lifetime_protein_chain_${chain}_events_df.csv
with the format:
resid_i,resid_j,start_frame,end_frame,frames,lifetime_ns
4,A,7090,7747,658,1316.00

Add a column with the original index, sort it by lfetime_ns, then read the trajectories with the format:
/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_${chain}/chain_${chain}_${index}/centers_aligned.pdb. 
The first frame of each pdb (biggest cluster) is saved in the the file /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_${chain}/all_centers_by_time.xtc
"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
import MDAnalysis as mda

chains = ["4", "r", "s"]

for chain in chains:
    # Path to the CSV file
    csv_path = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/lifetime/lifetime_protein_chain_{chain}_events_df.csv"
    
    if not os.path.exists(csv_path):
        print(f"CSV file not found: {csv_path}")
        continue
    
    # Read the dataframe
    df = pd.read_csv(csv_path)
    
    # Add original index column (starting from 1)
    df['original_index'] = df.index + 1
    
    # Sort by lifetime_ns descending
    df = df.sort_values('lifetime_ns', ascending=False)
    
    # Output XTC path
    output_xtc = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/all_centers_by_time.xtc"
    
    # List to hold universes
    universes = []
    
    for idx, row in df.iterrows():
        original_index = row['original_index']
        pdb_path = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/chain_{chain}_{original_index}/clust_c075/centers_aligned.pdb"
        
        if os.path.exists(pdb_path):
            u = mda.Universe(pdb_path)
            universes.append(u)
        else:
            print(f"PDB file not found: {pdb_path}")
    
    if universes:
        # Write to XTC
        with mda.Writer(output_xtc, universes[0].atoms.n_atoms) as w:
            for i, u in enumerate(universes):
                # Set time for the frame
                u.trajectory.ts.time = i * 1.0  # Assuming 1 ps per frame or adjust as needed
                w.write(u.atoms)
        
        print(f"Created XTC file: {output_xtc} with {len(universes)} frames")
    else:
        print(f"No valid PDB files found for chain {chain}")

    # Write a csv file with the pdb file paths and the lifetimes
    output_csv = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/all_centers_by_time.csv"
    df['pdb_path'] = df['original_index'].apply(lambda idx: f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/chain_{chain}_{idx}/clust_c075/centers_aligned.pdb")
    df.to_csv(output_csv, index=False)
    print(f"Wrote CSV file: {output_csv}")

print("Script completed.")