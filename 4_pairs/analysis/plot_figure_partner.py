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
def get_resname(pdb_path, chain_id, resid):
    u = mda.Universe(pdb_path)
    for ch in [chain_id, 'A', 'B']:
        sel = u.select_atoms(f"chainID {ch} and resid {resid}")
        if len(sel) > 0:
            return sel.residues[0].resname
    return None

# Function to plot bar data
def plot_bar_data(df, helices, title, color, pdb_path, chain_id, output_file):
    print(f"Plotting {title}: {len(df)} residues, max median_ns: {df['median_ns'].max() if not df.empty else 'N/A'}")
    fig, ax = plt.subplots(figsize=(12, 4))
    if df.empty:
        print(f"No data for {title}")
        plt.close(fig)
        return
    df = df.dropna(subset=['median_ns'])
    min_median = 0  # Start from 0 to avoid covering small bars
    max_median = df['median_ns'].max()
    for name, (start, end) in helices.items():
        ax.add_patch(plt.Rectangle((start, min_median), end - start, max_median - min_median, color='grey', alpha=0.1))
        ax.text((start + end) / 2, max_median, name, ha='center', va='bottom', fontsize=14)
    
    ax.bar(df['resid'], df['median_ns'], color=color, alpha=1.0, width=1)
    
    # Top 3 labels
    top3 = df.nlargest(3, 'median_ns')
    for _, row in top3.iterrows():
        if row['median_ns'] > 1000:
            resname = get_resname(pdb_path, chain_id, int(row['resid']))
            if resname:
                ax.text(row['resid'], row['median_ns'] + 50, f"{int(row['resid'])}{resname}", ha='center', va='bottom', fontsize=14, rotation=90)
            else:
                print(f"Resname not found for resid {int(row['resid'])} in {pdb_path} chain {chain_id}")
    
    # Horizontal line at 1 microsecond (1000 ns)
    ax.axhline(y=1000, color='red', linestyle='--', linewidth=1, label='1 μs')
    
    ax.set_xlabel('Resid')
    ax.set_ylabel('Median Lifetime (ns)')
    ax.set_title(title)
    if 'PsbS' in title:
        ax.set_ylim(0, 7000)
    else:
        ax.set_ylim(0, 3000)
    #ax.set_xlim(1, 1000)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close(fig)

# INPUT USER
chains = ["4", "r", "s"]

# Output directory
output_dir = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_per_binding_figures"
os.makedirs(output_dir, exist_ok=True)

# Colors from Set3
partner_colors = ['blue', 'green', 'red']
psbs_color = 'purple'

# Rename dicts
partner_rename = {'1': 'TM1', '2': 'TM2', '3': 'TM3', '4': 'D'}
psbs_rename = {'1': 'TM1', '2': 'H1', '3': 'TM2', '4': 'TM3', '5': 'H2', '6': 'TM4'}

# Partner PDB
partner_pdb = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb"
partner_helices = get_helices(partner_pdb, 'A')
partner_helices_renamed = {partner_rename.get(k, k): v for k, v in partner_helices.items()}

# PsbS PDB
psbs_pdb = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_4_0_dimer.pdb"
psbs_helices = get_helices(psbs_pdb, 'A')
psbs_helices_renamed = {psbs_rename.get(k, k): v for k, v in psbs_helices.items()}

for i, chain in enumerate(chains):
    base_path = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/lifetimes_chain_{chain}"
    
    # Partner dataframes
    chl_path = f"{base_path}/chain_{chain}_psbs_chl_residue_summary_df.csv"
    non_chl_path = f"{base_path}/chain_{chain}_psbs_non_chl_residue_summary_df.csv"
    partner_pdb_chain = f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_{chain}/initial_fit_merged.pdb"
    if os.path.exists(chl_path):
        df_chl = pd.read_csv(chl_path)
        plot_bar_data(df_chl, partner_helices_renamed, f'Chain {chain} Partner Chl', partner_colors[i], partner_pdb_chain, chain, f'{output_dir}/chain_{chain}_partner_chl_bar.png')
    if os.path.exists(non_chl_path):
        df_non_chl = pd.read_csv(non_chl_path)
        plot_bar_data(df_non_chl, partner_helices_renamed, f'Chain {chain} Partner Non-Chl', partner_colors[i], partner_pdb_chain, chain, f'{output_dir}/chain_{chain}_partner_non_chl_bar.png')
    
    # PsbS dataframes
    psbs_chl_path = f"{base_path}/psbs_chain_{chain}_chl_residue_summary_df.csv"
    psbs_non_chl_path = f"{base_path}/psbs_chain_{chain}_non_chl_residue_summary_df.csv"
    if os.path.exists(psbs_chl_path):
        df_psbs_chl = pd.read_csv(psbs_chl_path)
        plot_bar_data(df_psbs_chl, psbs_helices_renamed, f'Chain {chain} PsbS Chl', psbs_color, psbs_pdb, 'A', f'{output_dir}/chain_{chain}_psbs_chl_bar.png')
    if os.path.exists(psbs_non_chl_path):
        df_psbs_non_chl = pd.read_csv(psbs_non_chl_path)
        plot_bar_data(df_psbs_non_chl, psbs_helices_renamed, f'Chain {chain} PsbS Non-Chl', psbs_color, psbs_pdb, 'A', f'{output_dir}/chain_{chain}_psbs_non_chl_bar.png')

print("Bar plotting completed, output saved as PNG files.")