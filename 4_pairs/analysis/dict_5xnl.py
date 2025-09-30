"""
This script processes PDB files in a specified directory, assigns chain IDs based on a reference PDB file,
and saves modified PDB files with updated chain IDs.

"""

import MDAnalysis as mda
import os
import re
import warnings
warnings.filterwarnings("ignore")
import sys
dir=sys.argv[1]
print(f"Processing directory: {dir}")




ref_pdb="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb"

u = mda.Universe(ref_pdb)

chains = u.atoms.chainIDs

# Get unique values without sorting
unique_chains = list(set(chains))

n_residues = []
for chain in unique_chains:
    sel = u.select_atoms(f"chainID {chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB CHL)")
    n_r = sel.n_residues
    n_residues.append(n_r)

keys = unique_chains
values = n_residues
chain_dict = {keys[i]: values[i] for i in range(len(keys))}

# PDB files in the directory
pdb_files = [f for f in os.listdir(dir) if f.endswith('.pdb') and not f.startswith('fix') and not f.endswith('aligned.pdb')]
#pdb_files = [f for f in os.listdir(dir) if f.endswith('.pdb')]
#print(pdb_files)

#add path to the files
pdb_files = [os.path.join(dir, f) for f in pdb_files]
arr_of_arr_chains = []
arr_of_arr_residues = []
for file in pdb_files:
    u = mda.Universe(file)
    sel =  u.select_atoms(f"not resname *GG* *SQ* *PG* W* HOH *HEME* *HG* *MG* PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB CHL")
    chains = sel.atoms.chainIDs
    resnames = sel.atoms.resnames
    n_residues = sel.atoms.n_residues
    unique_chains = list(set(chains))
    unique_resnames = list(set(resnames))
    print(file)
    print(unique_chains)
    n_residues = [u.select_atoms(f"(not resname *GG* *SQ* *PG* W* HOH *HG* *MG* PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB CHL) and chainID {chain}").n_residues for chain in unique_chains]
    print(n_residues)
    # Get the basename after "A1_", "A2_", "A3_" or "A4_    
    basename = os.path.basename(file)
    # Get the basename after "A1_", "A2_", "A3_" or "A4_    
    name_without_ext = os.path.splitext(basename)[0]
    # Find the pattern A1, A2, A3, or A4 followed by chain IDs
    pattern = r'A[1-4]_(.+)'
    match = re.search(pattern, name_without_ext)
    if match:
        chains_part = match.group(1)
        # Split by underscore to get individual chain IDs
        chain_ids = chains_part.split('_')

    n_res_arr = []
    for chain in chain_ids:
        n_res = dict(chain_dict).get(chain, 0)
        print(f"Chain {chain}: {n_res} residues (from reference)")
        n_res_arr.append(n_res)

    keys = chain_ids
    values = n_residues
    dict_actual = {keys[i]: values[i] for i in range(len(keys))}

    # Create a dictionary with n_residues and the chainID assigned from the reference.
    # By default, if n_residues = 212 or 213, assign all to PsbS chain '9', else search for the closest value in chain_ids
    assigned_chain_dict = {}
    counter = 0
    for i in range(len(n_residues)):
        n = n_residues[i]
        c = unique_chains[i]
        if n in [212, 213] and counter == 0:
            assigned_chain_dict[c] = 'A'
            counter += 1
        elif n in [212, 213] and counter == 1:
            assigned_chain_dict[c] = 'B'
            counter += 1
        else:
            closest_chain = min(dict_actual.keys(), key=lambda k: abs(dict_actual[k] - n))
            assigned_chain_dict[c] = closest_chain
    print(f"Assigned chains: {assigned_chain_dict}")
    # Change the chainID to the assigned chainID
    # Work on a copy to avoid conflicts when chains overlap
    u_copy = u.copy()
    for old_chain, new_chain in assigned_chain_dict.items():
        print(f"Changing chain {old_chain} to {new_chain}")
        # Select atoms with the current chain ID from the original universe
        atoms_to_change = u.select_atoms(f"chainID {old_chain}")
        if len(atoms_to_change) > 0:
            # Find corresponding atoms in the copy and change their chainIDs
            copy_atoms = u_copy.atoms[atoms_to_change.indices]
            copy_atoms.chainIDs = [new_chain] * len(copy_atoms)
    # Save the modified PDB file
    new_filename = os.path.join(dir, f"fix_{basename}")
    u_copy.atoms.write(new_filename)
    print(f"Modified PDB saved as {new_filename}")
