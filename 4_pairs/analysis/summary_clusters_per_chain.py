import os
import MDAnalysis as mda
import pandas as pd
chains=["4", "r", "s"]
for chain in chains:
    dir_path=f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}"
    # Find all the pdbs with the name centers.pdb
    pdb_files = []
    for root, dirs, files in os.walk(dir_path):
        for f in files:
            if f.endswith("centers.pdb"):
                # Store relative path from dir_path for consistency with rest of code
                rel_path = os.path.relpath(os.path.join(root, f), dir_path)
                pdb_files.append(rel_path)
    if not pdb_files:
        print(f"No pdb files found in {dir_path} for chain {chain}. Skipping...")
        continue
    universes = [mda.Universe(os.path.join(dir_path, pdb)) for pdb in pdb_files]
    # The filename format is /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/chain_{chain}_{id}/clust_c{cutoff}
    if not universes:
        print(f"No valid universes found for chain {chain}. Skipping...")
        continue
    pdb_files_parent = [os.path.dirname(os.path.join(dir_path, pdb)) for pdb in pdb_files]
    headers = [os.path.basename(pdb).replace("clust_c", "clust_") for pdb in pdb_files_parent]
    values = [[universe.trajectory.n_frames] for universe in universes]
    # Flatten values to a simple list (since each inner list has only 1 element)
    n_frames_values = [v[0] for v in values]
    df = pd.Series(data=n_frames_values, index=headers, name='n_frames')
    df = df.sort_index()
    # Save the DataFrame to a CSV file
    output_dir = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses/chain_{chain}/"
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, f"summary_chain_{chain}.csv"), index=True, sep='\t')
    print(f"Summary for chain {chain} saved to {output_dir}/summary_chain_{chain}.csv")

