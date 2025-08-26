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
    contact_protein = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/contact_matrix_protein_chain{chain}_5000ns.csv", index_col=0)
    contact_cofactors = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/contact_matrix_cofactors_chain{chain}_5000ns.csv", index_col=0)
    contact_cofactors_resids = pd.read_csv(f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/contact_matrix_cofactors_resids_chain{chain}_5000ns.csv", index_col=0)
    resids_psbs = contact_protein.index.values

    contact_chain_psbs = contact_protein.sum(axis=1).values #Sum the rows, to get PsbS contacts
    contact_chain = contact_protein.sum(axis=0).values #Sum the columns, to chain contacts
    resids_chain = contact_protein.columns.values
    contact_chain_cofactors_psbs = contact_cofactors_resids.sum(axis=1).values #Sum the rows, to get PsbS contacts
    resids_chain_cofactors = contact_cofactors_resids.columns.values
    contact_chain_cofactors = contact_cofactors_resids.sum(axis=0).values #Sum the columns, to chain contacts


     #----Saving info in CHAIN 4, S, R, C------
    sel1 = u_ref.select_atoms(f'chainID {chain}')
    print(f"Selections had {sel1.n_residues} residues for chain {chain}")
    # Use variables directly instead of globals
    bfactors_chain = contact_chain
    print(f"Max value of contact_chain: ", bfactors_chain.max())
    bfactors_chain_cofactors = contact_chain_cofactors

    # Assign B-factors to the universe for the specific chain
    sel1 = assign_bfactor_to_universe(sel1, resids_chain, bfactors_chain)
    sel1 = assign_bfactor_to_universe(sel1, resids_chain_cofactors, bfactors_chain_cofactors)
    # Save the universe to a PDB file
    output_filename = f'/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/contacts_5000ns.pdb'
    save_universe_to_pdb(sel1, output_filename)
    print(f"B-factors assigned and saved for chain {chain} in {output_filename}")

    #---Write for the other identical chain
    sel2 = u_ref.select_atoms(f'chainID {same}')
    print(f"Selections had {sel2.n_residues} residues for chain {same}")
    # Use variables directly instead of globals
    bfactors_chain = contact_chain
    print(f"Max value of contact_chain: ", bfactors_chain.max())
    bfactors_chain_cofactors = contact_chain_cofactors

    # Assign B-factors to the universe for the specific chain
    sel2 = assign_bfactor_to_universe(sel2, resids_chain, bfactors_chain)
    sel2 = assign_bfactor_to_universe(sel2, resids_chain_cofactors, bfactors_chain_cofactors)
    # Save the universe to a PDB file
    output_filename = f'/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/contacts_5000ns_{same}.pdb'
    save_universe_to_pdb(sel2, output_filename)
    print(f"B-factors assigned and saved for chain {same} in {output_filename}")

    #----SAVE DATA TO PSBS------
    psbs_gro = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_4_0_dimer.pdb"
    u_psbs = mda.Universe(psbs_gro)
    # Use variables directly instead of globals
    #Sum bfactors chain and cofactors
    bfactors_chain_psbs = contact_chain_psbs
    bfactors_chain_cofactors_psbs = contact_chain_cofactors_psbs
    #Sum the two arrays
    bfactors_chain = bfactors_chain_psbs + bfactors_chain_cofactors_psbs
    #bfactors_chain = bfactors_chain_psbs 

    print(f"Max value of contact_chain{chain}_psbs: ", bfactors_chain.max())

    # Assign B-factors to the universe for the specific chain
    u_psbs = assign_bfactor_to_universe(u_psbs, resids_psbs, bfactors_chain)

    # Save the universe to a PDB file
    output_filename = f'/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/psbs_chain{chain}_contacts_5000ns.pdb'

    save_universe_to_pdb(u_psbs, output_filename)
    print(f"B-factors assigned and saved for PsbS chain {chain} in {output_filename}")




    
