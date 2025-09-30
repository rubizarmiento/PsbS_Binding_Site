"""
Print the dictionary for the chainI

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

# Get unique values preserving order
unique_chains = list(dict.fromkeys(chains))

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
    # Preserve order of chains as they appear in the file (important for A/B at the end)
    unique_chains = list(dict.fromkeys(chains))  # dict.fromkeys preserves order
    unique_resnames = list(set(resnames))
    print(file)
    print(unique_chains)
    n_residues = [u.select_atoms(f"(not resname *GG* *SQ* *PG* W* HOH *HG* *MG* PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB CHL) and chainID {chain}").n_residues for chain in unique_chains]
    print(n_residues)
    # Get the basename after "A1_", "A2_", "A3_" or "A4_    
    basename = os.path.basename(file)
    # Get the basename after "A1_", "A2_", "A3_" or "A4_    
    name_without_ext = os.path.splitext(basename)[0]
    # Find the pattern: number followed by underscore and chain IDs
    pattern = r'(\d+)_(.+)'
    match = re.search(pattern, name_without_ext)
    chain_ids = []  # Initialize as empty list
    if match:
        chains_part = match.group(2)
        # Split by underscore to get individual chain IDs
        chain_ids = chains_part.split('_')
    else:
        print(f"Warning: Could not parse chain IDs from filename {name_without_ext}")
        # Try to use unique chains from the PDB file as fallback
        chain_ids = unique_chains

    n_res_arr = []
    for chain in chain_ids:
        n_res = dict(chain_dict).get(chain, 0)
        print(f"Chain {chain}: {n_res} residues (from reference)")
        n_res_arr.append(n_res)

    # Create a dictionary mapping actual chains to their residue counts
    dict_actual = {unique_chains[i]: n_residues[i] for i in range(len(unique_chains))}

    # Identify PsbS candidates (chains with 212 or 213 residues)
    psbs_candidates = [i for i, n in enumerate(n_residues) if n in [212, 213]]

    # Create a dictionary with n_residues and the chainID assigned from the reference.
    # Prioritize last two PsbS candidates as A and B
    assigned_chain_dict = {}
    if len(psbs_candidates) >= 2:
        # Assign the last two PsbS candidates: second last as A, last as B
        a_idx = psbs_candidates[-2]
        b_idx = psbs_candidates[-1]
        assigned_chain_dict[unique_chains[a_idx]] = 'A'
        assigned_chain_dict[unique_chains[b_idx]] = 'B'
        if len(psbs_candidates) > 2:
            warnings.warn(f"More than 2 PsbS candidates found ({len(psbs_candidates)}), assigning last two as A and B")
    elif len(psbs_candidates) == 1:
        # Only one PsbS candidate, assign as A
        assigned_chain_dict[unique_chains[psbs_candidates[0]]] = 'A'
        warnings.warn("Only one PsbS candidate found, assigned as A")
    else:
        warnings.warn("No PsbS candidates found (212/213 residues)")

    # For remaining chains, assign closest match from reference
    for i in range(len(n_residues)):
        c = unique_chains[i]
        if c not in assigned_chain_dict:
            n = n_residues[i]
            if chain_ids:  # Only if we have chain_ids from filename
                closest_chain = min(chain_ids, key=lambda k: abs(dict(chain_dict).get(k, 0) - n))
                assigned_chain_dict[c] = closest_chain
            else:
                # Fallback: keep original chain ID
                assigned_chain_dict[c] = c
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
