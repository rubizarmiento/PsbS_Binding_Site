import os
import sys
from lib_contacts import *

#Input
pdb_dir = sys.argv[1]
lifetimes_dir = sys.argv[2]
odir = sys.argv[3]

# Create dir
os.makedirs(odir, exist_ok=True)
cases = [d for d in os.listdir(pdb_dir) if os.path.isdir(os.path.join(pdb_dir, d))] # Get all directories in the cases_dir

# Check if cases were found
if len(cases) == 0:
    raise ValueError(f"No case directories found in {pdb_dir}. Please check the path.")
    exit()
else:
    print(f"Found {len(cases)} case directories in {pdb_dir}.")

for i in range(len(cases)):
    case=cases[i]
    #Split chain name
    chains=case.split("_")[1:]
    file = f"{pdb_dir}/{case}/FINAL/final.pdb" 
    #Check if file exists
    if not os.path.isfile(file):
        print(f"File {file} does not exist. Skipping...")
        continue
    u = mda.Universe(file)
    
    # Initialize
    chain0 = chains[0]
    contact_protein = pd.read_csv(f"{lifetimes_dir}/chain_{chain0}_{case}_residue_summary_df.csv") 

    selelections = u.select_atoms(f'chainID {chain0}')
    resids_chain = contact_protein["resid"]
    bfactors_chain = contact_protein["median_ns"]
    
    print(len(resids_chain))
    print(resids_chain.values)
    exit()
    # If resid is A1, A2, A3, A4, replace it for 9
    if 'A' in resids_chain.values:
        print("Replacing A1, A2, A3, A4 with 9")
        resids_chain = resids_chain.replace({'A1': '9', 'A2': '9', 'A3': '9', 'A4': '9'})


    # Assign B-factors to the universe for the specific chain
    selections = assign_bfactor_to_universe(selelections, resids_chain, bfactors_chain)


    for chain in chains[1:]:
        #-----Getting DATA---------
        contact_protein = pd.read_csv(f"{lifetimes_dir}/chain_{chain}_{case}_residue_summary_df.csv") 

        sel1 = u.select_atoms(f'chainID {chain}')

        resids_chain = contact_protein["resid"]
        bfactors_chain = contact_protein["median_ns"]
        # Assign B-factors to the universe for the specific chain
        sel1 = assign_bfactor_to_universe(sel1, resids_chain, bfactors_chain)
        selections += sel1

    # Save the universe to a PDB file
    output_filename = f"{odir}/lifetimes_{case}.pdb" #TODO
    save_universe_to_pdb(selections, output_filename)
    print(f"B-factors assigned and saved for chain {chains} in {output_filename}")


    #--- SAVE COFACTORS ---

    # Initialize
    chain0 = chains[0]
    f2 = f"{pdb_dir}/{case}_cofactors_cg.pdb"

    # Check if file exists
    if not os.path.isfile(f2):
        print(f"File {f2} does not exist. Skipping...")
        continue
    else:
        u2 = mda.Universe(f2)
        contact_protein = pd.read_csv(f"{lifetimes_dir}/chain_{chain0}_{case}_residue_summary_df.csv")
        resids_chain = contact_protein["resid"]
        bfactors_chain = contact_protein["median_ns"]
        selections2 = u2.select_atoms(f'chainID {chain0}')
        # Assign B-factors to the universe for the specific chain
        selections2 = assign_bfactor_to_universe(selections2, resids_chain, bfactors_chain)


        for chain in chains[1:]:
            #-----Getting DATA---------
            contact_protein = pd.read_csv(f"{lifetimes_dir}/chain_{chain}_{case}_residue_summary_df.csv")
            sel1 = u2.select_atoms(f'chainID {chain}')

            resids_chain = contact_protein["resid"]
            bfactors_chain = contact_protein["median_ns"]
            # Assign B-factors to the universe for the specific chain
            sel1 = assign_bfactor_to_universe(sel1, resids_chain, bfactors_chain)
            selections2 += sel1

        # Save the universe to a PDB file
        output_filename = f"{odir}/lifetimes_cofactors_{case}.pdb"
        save_universe_to_pdb(selections2, output_filename)
        print(f"B-factors assigned and saved for chain {chains} in {output_filename}")


