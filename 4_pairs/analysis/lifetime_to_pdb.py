from lib_contacts import *

chains=["4","r","c","s"]
sames=["8","R","C","S"]
#chains=["4"]

#Reference files
ref_gro = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb"
# Save the AA chains with B-factors
u_ref = mda.Universe(ref_gro)

for i in range(len(chains)):
    chain=chains[i]
    same=sames[i]
    #-----Getting DATA---------
    contact_protein = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/lifetime/lifetime_chain{chain}_residue_summary_df.csv")
    contact_cofactors = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/lifetime/lifetime_cofactors_chain_{chain}_residue_summary_df.csv")

    sel1 = u_ref.select_atoms(f'chainID {chain}')

    resids_chain = contact_protein["resid"]
    resids_chain_cofactors = contact_cofactors["resid"]

    bfactors_chain = contact_protein["median_ns"]
    bfactors_chain_cofactors = contact_cofactors["median_ns"]

    print("Max chain: ", bfactors_chain.max())
    print("Max chain cofactors: ", bfactors_chain_cofactors.max())

    # Assign B-factors to the universe for the specific chain
    sel1 = assign_bfactor_to_universe(sel1, resids_chain, bfactors_chain)
    sel1 = assign_bfactor_to_universe(sel1, resids_chain_cofactors, bfactors_chain_cofactors)
    # Save the universe to a PDB file
    output_filename = f'/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/lifetime_5000ns.pdb'
    save_universe_to_pdb(sel1, output_filename)
    print(f"B-factors assigned and saved for chain {chain} in {output_filename}")

    #---Write for the other identical chain
    sel2 = u_ref.select_atoms(f'chainID {same}')
    # Assign B-factors to the universe for the specific chain
    sel2 = assign_bfactor_to_universe(sel2, resids_chain, bfactors_chain)
    sel2 = assign_bfactor_to_universe(sel2, resids_chain_cofactors, bfactors_chain_cofactors)
    # Save the universe to a PDB file
    output_filename = f'/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/lifetime_5000ns_{same}.pdb'
    save_universe_to_pdb(sel2, output_filename)
    print(f"B-factors assigned and saved for chain {same} in {output_filename}")

    #----SAVE DATA TO PSBS------
    psbs_gro = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_4_0_dimer.pdb"
    u_psbs = mda.Universe(psbs_gro)
    contact_chain_psbs = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/lifetime/lifetime_psbs_chain_{chain}_residue_summary_df.csv")    
    bfactors_chain =  contact_chain_psbs["median_ns"]
    print("Max PsbS: ", bfactors_chain.max())
    resids_psbs =  contact_chain_psbs["resid"]
    
    # Assign B-factors to the universe for the specific chain
    u_psbs = assign_bfactor_to_universe(u_psbs, resids_psbs, bfactors_chain)

    # Save the universe to a PDB file
    output_filename = f'/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/psbs_chain{chain}_lifetime_5000ns.pdb'

    save_universe_to_pdb(u_psbs, output_filename)
    print(f"B-factors assigned and saved for PsbS chain {chain} in {output_filename}")




    
