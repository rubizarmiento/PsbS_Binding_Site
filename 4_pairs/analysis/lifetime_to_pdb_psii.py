import os
import sys
from lib_contacts import *

#Input
pdb_dir = sys.argv[1]
lifetimes_dir = sys.argv[2]
odir = sys.argv[3]
sel_protein = sys.argv[4]
sel_cofactors = sys.argv[5]
sel_psbs = sys.argv[6]


# Create dir
os.makedirs(odir, exist_ok=True)
cases = [d for d in os.listdir(pdb_dir) if os.path.isdir(os.path.join(pdb_dir, d))] # Get all directories in the cases_dir

# Check if cases were found
if len(cases) == 0:
    raise ValueError(f"No case directories found in {pdb_dir}. Please check the path.")
    exit()
else:
    print(f"Found {len(cases)} case directories in {pdb_dir}.")
#cases = ["4_7_8"]




dict = { 
    "type": ["PROTEIN", "COFACTORS", "PSBS"],
    "input_pdb": [f"{pdb_dir}/{{case}}/FINAL/final_aligned.pdb", f"{pdb_dir}/{{case}}_cofactors_cg.pdb", "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_4_0_dimer_aligned.pdb"],
    "input_csv": [f"{lifetimes_dir}/chain_{{chain}}_{{case}}_residue_summary_df.csv", f"{lifetimes_dir}/chain_{{chain}}_{{case}}_residue_summary_df.csv", f"{lifetimes_dir}/psbs_{{case}}_residue_summary_df.csv"],
    "selections": [f'chainID {{chain}} and {sel_protein}', f'chainID {{chain}} and {sel_cofactors}', sel_psbs],
    "output": [f"{odir}/lifetimes_{{case}}_protein.pdb", f"{odir}/lifetimes_{{case}}_cofactors.pdb", f"{odir}/lifetimes_{{case}}_psbs.pdb"]
}



# Proteins
for i in range(len(cases)):
    case=cases[i]
    for j in range(len(dict["input_pdb"])):
        type = dict["type"][j]
        pdb = dict["input_pdb"][j].format(case=case)
        output = dict["output"][j].format(case=case)

        print(f"\n\n---TYPE: {type} {case}---")

        #Split chain name
        if type == "PSBS":
            chains = ["A B"]
        else:
            chains=case.split("_")[1:]

        if os.path.isfile(pdb):
            u = mda.Universe(pdb)
        else:
            print(f"{pdb} not found. Skipping...")
            continue 

        sel = [u.select_atoms("all")]

        for chain in chains:
            print(f"---Chain {chain}---")
            selection = dict["selections"][j].format(chain=chain)
            csv_file = dict["input_csv"][j].format(case=case,chain=chain)
            #-----Getting DATA---------
            if not os.path.isfile(csv_file):
                print(f"{csv_file} does not exist. Skipping...")
            else:
                print(f"Analyzing {csv_file} file")
                contacts = pd.read_csv(csv_file) 

                resids_chain = contacts["resid"]
                bfactors_chain = contacts["median_ns"]


                sel_ = u.select_atoms(selection) 
                if len(sel_.atoms) == 0:
                    print(f" Selection {selection} is empty. Skipping...")
                else:
                    # Separate resids
                    resids = set(sel_.atoms.resids)

                    filtered_resids = []
                filtered_bfactors = []

                for i in range(len(resids_chain)):
                    if resids_chain[i] in resids:
                        filtered_resids.append(resids_chain[i])
                        filtered_bfactors.append(bfactors_chain[i])
                print(f"Found {len(filtered_resids)} proteins residues in chain {chain} csv")

                # Assign B-factors to the universe for the specific chain
                sel_ = u.select_atoms(selection) 
                if len(sel_.atoms) == 0:
                    print(f" Selection {selection} is empty. Skipping...") 

                sel_ = assign_bfactor_to_universe(sel_, filtered_resids, filtered_bfactors)
                sel.append(sel_)


        if len(sel) > 1:
            sels= sel[1]
            for i in range(len(sel)):
                if i == 0:
                    continue
                else:
                    sels += sel[i]
            # Save the universe to a PDB file
            save_universe_to_pdb(sels, output)
