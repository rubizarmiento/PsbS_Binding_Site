from lib_contacts import *

#Reference files
dir="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/biggest_clusters_c075"
cases = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]


for i in range(len(cases)):
    case=cases[i]
    #Split chain name
    chains=case.split("_")[1:]
    file = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/biggest_clusters_c075/{case}/FINAL/final.pdb"
    #Check if file exists
    if not os.path.isfile(file):
        print(f"File {file} does not exist. Skipping...")
        continue
    u = mda.Universe(file)
    
    # Initialize
    chain0 = chains[0]
    contact_protein = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes/chain_{chain0}_{case}_residue_summary_df.csv")
    selelections = u.select_atoms(f'chainID {chain0}')
    resids_chain = contact_protein["resid"]
    bfactors_chain = contact_protein["median_ns"]
    # Assign B-factors to the universe for the specific chain
    selections = assign_bfactor_to_universe(selelections, resids_chain, bfactors_chain)


    for chain in chains[1:]:
        #-----Getting DATA---------
        contact_protein = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes/chain_{chain}_{case}_residue_summary_df.csv")
        sel1 = u.select_atoms(f'chainID {chain}')

        resids_chain = contact_protein["resid"]
        bfactors_chain = contact_protein["median_ns"]
        # Assign B-factors to the universe for the specific chain
        sel1 = assign_bfactor_to_universe(sel1, resids_chain, bfactors_chain)
        selections += sel1

    # Save the universe to a PDB file
    output_filename = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/pdbs_lifetime_psii/lifetimes_{case}.pdb"
    save_universe_to_pdb(selections, output_filename)
    print(f"B-factors assigned and saved for chain {chains} in {output_filename}")


    #--- SAVE COFACTORS ---

    # Initialize
    chain0 = chains[0]
    f2 = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/biggest_clusters_c075/{case}_cofactors.pdb"
    u2 = mda.Universe(f2)
    contact_protein = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes/chain_{chain0}_{case}_residue_summary_df.csv")
    resids_chain = contact_protein["resid"]
    bfactors_chain = contact_protein["median_ns"]
    selections2 = u2.select_atoms(f'chainID {chain0}')
    # Assign B-factors to the universe for the specific chain
    selections2 = assign_bfactor_to_universe(selections2, resids_chain, bfactors_chain)


    for chain in chains[1:]:
        #-----Getting DATA---------
        contact_protein = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes/chain_{chain}_{case}_residue_summary_df.csv")
        sel1 = u2.select_atoms(f'chainID {chain}')

        resids_chain = contact_protein["resid"]
        bfactors_chain = contact_protein["median_ns"]
        # Assign B-factors to the universe for the specific chain
        sel1 = assign_bfactor_to_universe(sel1, resids_chain, bfactors_chain)
        selections2 += sel1

    # Save the universe to a PDB file
    output_filename = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/pdbs_lifetime_psii/lifetimes_cofactors_{case}.pdb"
    save_universe_to_pdb(selections2, output_filename)
    print(f"B-factors assigned and saved for chain {chains} in {output_filename}")


