"""
Reads a csv file and copies the representative pdbs and tpr files for the type of binding event

With the columns:
tag original     new          old_chains count tag_number
n_s sim_4_A4_N_S sim_4_A4_n_s N_S        6     1


Then simply copies the first {original}.*pdb and {original}.*tpr files in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
With the name 
{tag_number}_{tag}

"""

import MDAnalysis as mda
import os
import sys
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
#Ignore warnings
import warnings
import glob
import shutil
from pathlib import Path
warnings.filterwarnings("ignore")

pdb_dir=sys.argv[1]
csv_path=sys.argv[2]
odir = sys.argv[3]

print(f"Processing directory: {pdb_dir} and CSV file: {csv_path}")

# Create output directory if it doesn't exist
if not os.path.exists(odir):
    os.makedirs(odir)

# Delete contents of the output directory
for filename in os.listdir(odir):
    file_path = os.path.join(odir, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
    except Exception as e:
        print(f'Failed to delete {file_path}. Reason: {e}') 



# Read the CSV
df = pd.read_csv(csv_path,header=0,sep=' ')

# Drop duplicated tag rows, keeping the first occurrence
df = df.drop_duplicates(subset=['tag'], keep='first')

for index, row in df.iterrows():
    chains_str = row['tag']  # First column: chains (e.g., 'A1_A2')
    basename = row['original']    # Second column: basename (e.g., 'sim_1_A2_Z')
    tag_number = row['tag_number']
    old_chains_str = row['old_chains']
    output = row['tag']
    
    # Input PDB path
    pdb_pattern = os.path.join(pdb_dir, f'{basename}.pdb')
    pdb_files = glob.glob(pdb_pattern)
    if not pdb_files:
        print(f"Warning: No files found matching pattern {pdb_pattern}, skipping.")
        continue
    print(f"Found PDB files: {pdb_files[0]}")

    pdb_path = pdb_files[0]  # Use the first matching file

    # Output grouped PDB path
    grouped_pdb_path = os.path.join(odir, f'{tag_number}_{output}.pdb')

    
    if not os.path.exists(pdb_path):
        print(f"Warning: {pdb_path} not found, skipping.")
        continue
    
    # Load the universe
    u = mda.Universe(pdb_path)

    if u.atoms.n_atoms == 0:
        print(f"Warning: {pdb_path} has no atoms, skipping.")
        continue

    sel_new = chains_str.split('_')
    sel_old = old_chains_str.split('_')

    # Dictionary to map new chains to old chains
    chain_map = dict(zip(sel_old, sel_new))
    #print(f"Mapping chains: {chain_map}")
    # Change chainIDs
    for atom in u.atoms:
        orig_chainID = atom.chainID
        if orig_chainID in chain_map:
            new_chainID = chain_map[orig_chainID]
            atom.chainID = new_chainID
    
    # Write the selected atoms to a new PDB
    u.select_atoms("all").write(grouped_pdb_path)
    
    print(f"Created grouped PDB: {grouped_pdb_path}")

    # Copy the corresponding TPR file
    tpr_pattern = os.path.join(pdb_dir, f'{basename}.tpr')
    tpr_files = glob.glob(tpr_pattern)
    if not tpr_files:
        print(f"Warning: No TPR files found matching pattern {tpr_pattern}, skipping TPR copy.")
        continue
    tpr_path = tpr_files[0]  # Use the first matching file
    if not os.path.exists(tpr_path):
        print(f"Warning: {tpr_path} not found, skipping TPR copy.")
        continue
    grouped_tpr_path = os.path.join(odir, f'{tag_number}_{output}.tpr')
    
    tpr_cofactors = tpr_path.replace('.tpr','_cofactors.tpr')
    if os.path.exists(tpr_cofactors):
        shutil.copyfile(tpr_cofactors, grouped_tpr_path.replace('.tpr','_cofactors.tpr'))
        print(f"Copied grouped TPR cofactors: {grouped_tpr_path.replace('.tpr','_cofactors.tpr')}")

    shutil.copyfile(tpr_path, grouped_tpr_path)
    print(f"Copied grouped TPR: {grouped_tpr_path}")

print("Grouping completed.")