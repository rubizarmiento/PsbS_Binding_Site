"""
1) Loads the dataframes for the chains 4, r, s:
/martini/rubiz/Github/PsbS_Binding_Site/4_pair# Pa# Rename dicts
partner_rename = {'AA1': 'TM1', 'AA2': 'TM2', 'AA3': 'TM3', 'AA4': 'D'}
psbs_rename = {'1': 'TM1', '2': 'H1', '3': 'TM2', '4': 'TM3', '5': 'H2', '6': 'TM4'}er PDB
partner_pdb = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb"
partner_helices = get_helices(partner_pdb, '1')
partner_helices_renamed = {partner_rename.get(k, k): v for k, v in partner_helices.items()}

# PsbS PDB
psbs_pdb = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/4ri2.pdb"
psbs_helices = get_helices(psbs_pdb, 'A')
psbs_helices_renamed = {psbs_rename.get(k, k): v for k, v in psbs_helices.items()}s/top_lifetimes_pdbs/lifetimes_chain_{chain}/chain_{chain}_psbs_chl_residue_summary_df.csv
/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/lifetimes_chain_{chain}/chain_{chain}_psbs_non_chl_residue_summary_df.csv

2) Uses Biopython to parse the helix section the PDB:
/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb
 and stores the helixid as keys and the first and last resid as values of helices of the {chain}.

3) Renames the keys of the dictionary depending on the chain:
chain_4: TM1, TM2, TM3, D 
chain_r: TM1, TM2, TM3, D
chain_s: TM1, TM2, TM3, D

4) Generates two figures of a bar plot of the "resid" vs "median_ns" of the {chain} and as background, adds a grey block with the helix labels and the limits given by the resids ranges. Uses the first three colors of cmap Set3, one per chain. Adds labels to the top 5 residues by getting the resname from the PDB file using mdanalysis and the chain selection.

5) Loads the dataframes:
/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/lifetimes_chain_{chain}/psbs_chain_{chain}_chl_residue_summary_df.csv
/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/lifetimes_chain_{chain}/psbs_chain_{chain}_non_chl_residue_summary_df.csv

6) Uses Biopython to parse the helix section the PDB:
/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/4ri2.pdb
 and stores the helixid as keys and the first and last resid as values of helices of the {chain}.

7) Renames the keys of the dictionary:
psbs: TM1, H1, TM2, TM3, H2, TM4 

8) Generates two figures of a bar plot of the "resid" vs "median_ns" of the psbs and as background, adds a grey block with the helix labels and the limits given by the resids ranges. Uses the fourth color of cmap Set3. Adds labels to the top 5 residues by getting the resname from the PDB file using mdanalysis and the chain selection.

"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
import MDAnalysis as mda
from adjustText import adjust_text
#Ignore warnings
import warnings
warnings.filterwarnings("ignore")
# Function to get helices for a chain
def get_helices(pdb_path, chain_id):
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_path)
    helices = {}
    if 'HELIX' not in structure.header:
        print(f"No helix information found in {pdb_path} header.")
        return helices
    for helix in structure.header['HELIX']:
        if helix['chain_id'] == chain_id:
            helices[helix['helix_id']] = (helix['start'], helix['end'])
    return helices

# Function to get resname

# Function to get resname
def get_resname(u, chain_id, resid):
    #chain_id is an array
    joined_chain_id = ' '.join(chain_id)
    sel = u.select_atoms(f"chainID {joined_chain_id} and resid {resid}")
    return sel.residues[0].resname

# Function to plot bar data
def plot_bar_data(df,  title, color, u, chain_id, output_file,helices=False):
    print(f"Plotting {title}: {len(df)} residues, max median_ns: {df['median_ns'].max() if not df.empty else 'N/A'}")
    fig, ax = plt.subplots(figsize=(12, 4))

    df = df.dropna(subset=['median_ns'])
    min_median = 0  # Start from 0 to avoid covering small bars
    max_median = df['median_ns'].max()
    if helices:
        for name, (start, end) in helices.items():
            ax.add_patch(plt.Rectangle((start, min_median), end - start, max_median - min_median, color='grey', alpha=0.1))
            ax.text((start + end) / 2, max_median, name, ha='center', va='bottom', fontsize=14)
    
    bars = ax.bar(df['resid'], df['median_ns'], color=color, alpha=1.0, width=1)
    
    # Top 5 labels
    top5 = df.nlargest(10, 'median_ns')
    # Filter out the only the ones larger than 1000 ns
    top5 = top5[top5['median_ns'] > 1000]
    texts = []
    for _, row in top5.iterrows():
        resname = get_resname(u, chain_id, int(row['resid']))
        if resname:
            texts.append(ax.text(row['resid'], row['median_ns'], f"{int(row['resid'])}{resname}", ha='center', va='bottom', fontsize=14, rotation=90))
        else:
            print(f"Resname not found for resid {int(row['resid'])} in {u} chain {chain_id}")
    #adjust_text(texts,arrowprops=dict(arrowstyle="-", color='k', lw=0.5), only_move={"text": '+y',"explode": "x"},expand=(3, 3), force_text=10);

    adjust_text(texts,arrowprops=dict(arrowstyle="-", color='k', lw=0.5),expand=(1.5, 1.5),  only_move={"text": '+x',"explode": "xy"} , force_text=10);

    
    # Horizontal line at 1 microsecond (1000 ns)
    ax.axhline(y=2000, color='lightgrey', linestyle='--', linewidth=1, label='1 μs')
    
    ax.set_xlabel('Resid')
    ax.set_ylabel('Median Lifetime (ns)')
    ax.set_title(title)
    if 'PsbS' in title:
        ax.set_ylim(0, 3000)
    else:
        ax.set_ylim(0, 3000)
    #ax.set_xlim(1, 1000)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close(fig)

# INPUT USER

dir = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/biggest_clusters_c075"
output_dir = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_per_binding_figures_psii"
ref_pdb_chains = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb"
ref_pdb_psbs = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_dimer_cg.pdb"
u_chains = mda.Universe(ref_pdb_chains)
u_psbs = mda.Universe(ref_pdb_psbs)
# Colors from Set3
partner_colors = ['blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red',]
psbs_color = 'purple'



# Output directory
os.makedirs(output_dir, exist_ok=True)
cases = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]

for i in range(len(cases)):
    case=cases[i]
    chains=case.split("_")[1:]
    for i, chain in enumerate(chains):
        chain_path = f"/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes/chain_{chain}_{case}_residue_summary_df.csv"
        # Used for labelling top residues
        
        df = pd.read_csv(chain_path)
        plot_bar_data(df, f'{chain}', partner_colors[i], u_chains, [chain], f'{output_dir}/chain_{chain}_{case}.png')
        
    # PsbS dataframes
    psbs_path = f"/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes/psbs_{case}_residue_summary_df.csv"
    df_psbs = pd.read_csv(psbs_path)
    plot_bar_data(df_psbs, f'Case {case} PsbS', psbs_color, u_psbs, ['A','B'], f'{output_dir}/psbs_{case}.png')
