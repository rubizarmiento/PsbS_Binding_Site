import yaml
import os   
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
import MDAnalysis as mda
from Bio.PDB import MMCIFParser
import matplotlib.pyplot as plt
import numpy as np

#Ignore warnings
import warnings
warnings.filterwarnings("ignore")

chain_labels_yaml = "/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_labels.yaml"
basenames_csv = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv"
cif_protein = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes/1_n_s_protein.cif"
pdb_cofactors = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes/1_n_s_cofactors.pdb"
color_definitions_yaml = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/color_definitions.yaml"
output_figure = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/test.png"

# If CIF file contains protein, get the case
basename_cif = os.path.basename(cif_protein)
case_plot = basename_cif.split("_protein.cif")[0]
chains = case_plot.split("_")[1:]  # Extract chains from case_plot

if basename_cif.endswith("_protein.cif"):
    type_plot = "one_letter"


# Read YAML file
with open(chain_labels_yaml, 'r') as f:
    data = yaml.safe_load(f)

# Read color definitions YAML file
with open(color_definitions_yaml, 'r') as f:
    color_config = yaml.safe_load(f)

# Helper function to extract nested config values
def get_config_value(config, path, default=None):
    """Extract nested value from config using dot-separated path.
    
    Example: get_config_value(config, 'text.sequence.letters.fontsize')
    """
    keys = path.split('.')
    value = config
    try:
        for key in keys:
            value = value[key]
        return value
    except (KeyError, TypeError):
        return default

# Define the config structure we need (variable_name -> yaml_path)
config_mapping = {
    # Colormap
    'cmap_name': 'colormap.name',
    'vmin': 'colormap.normalization.vmin', 
    'vmax': 'colormap.normalization.vmax',
    
    # Layout
    'box_width': 'layout.box_width',
    'box_height': 'layout.box_height', 
    'box_spacing': 'layout.box_spacing',
    'sequence_y': 'layout.positions.sequence_y',
    'resid_labels_y': 'layout.positions.resid_labels_y',
    'helix_boxes_y': 'layout.positions.helix_boxes_y',
    'helix_labels_y': 'layout.positions.helix_labels_y',
    'plot_title_y': 'layout.positions.plot_title_y',
    'cofactor_y': 'layout.positions.cofactor_y',
    
    # Text styling - sequence letters
    'seq_letter_fontsize': 'text.sequence.letters.fontsize',
    'seq_letter_color': 'text.sequence.letters.color',
    'seq_letter_fontweight': 'text.sequence.letters.fontweight',
    
    # Text styling - residue labels  
    'resid_label_fontsize': 'text.sequence.resid_labels.fontsize',
    'resid_label_color': 'text.sequence.resid_labels.color',
    'resid_label_fontweight': 'text.sequence.resid_labels.fontweight',
    
    # Text styling - cofactors
    'cofactor_fontsize': 'text.cofactors.labels.fontsize',
    'cofactor_color': 'text.cofactors.labels.color',
    'cofactor_fontweight': 'text.cofactors.labels.fontweight',
    
    # Text styling - helix labels
    'helix_label_fontsize': 'text.annotations.helix_labels.fontsize',
    'helix_label_color': 'text.annotations.helix_labels.color',
    'helix_label_fontweight': 'text.annotations.helix_labels.fontweight',
    
    # Text styling - region labels
    'region_label_fontsize': 'text.annotations.region_labels.fontsize',
    'region_label_color': 'text.annotations.region_labels.color',
    'region_label_fontweight': 'text.annotations.region_labels.fontweight',
    
    # Text styling - plot title
    'plot_title_fontsize': 'text.annotations.plot_title.fontsize',
    'plot_title_color': 'text.annotations.plot_title.color',
    'plot_title_fontweight': 'text.annotations.plot_title.fontweight',
    
    # Box styling
    'helix_facecolor': 'boxes.helices.facecolor',
    'helix_edgecolor': 'boxes.helices.edgecolor',
    'helix_alpha': 'boxes.helices.alpha'
}

# Extract all config variables at once
globals().update({var_name: get_config_value(color_config, path) 
                  for var_name, path in config_mapping.items()})

# Map chain IDs to their labels if chain not in yaml add Psb before the chain ID
chain_id_to_label = {}
for chain_id in chains:
    if chain_id in data:
        chain_id_to_label[chain_id] = data[chain_id]
    else:
        chain_id_to_label[chain_id] = f"Psb{chain_id}"

# Read basenames.csv file as dictionary
df = pd.read_csv(basenames_csv, header=0, sep=' ')
chains_case_dict = pd.Series(df['unique_basename'].values, index=df['tag']).to_dict()
#TODO: Function, dict

# Load the CIF file using Biopython
# One-letter code mapping
def resnames_arr_to_one_letter_arr(resnames):
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    one_letter_codes = [three_to_one.get(resname, 'X') for resname in resnames]  # X for unknown
    return one_letter_codes
#TODO: Function, dict




parser = MMCIFParser()
structure = parser.get_structure('protein', cif_protein)
one_letter_arr_of_arr = []
bfactors_arr_of_arr = []
chains = []
for model in structure:
    for chain in model:
        resnames_arr = []
        bfactors_arr = []
        for residue in chain:
            resname = residue.get_resname()
            resnames_arr.append(resname)
            # Get B-factors from all atoms in the residue and average them
            bfactors = [atom.get_bfactor() for atom in residue]
            avg_bfactor = sum(bfactors) / len(bfactors) if bfactors else 0.0
            bfactors_arr.append(avg_bfactor)
        one_letter_seq_arr = resnames_arr_to_one_letter_arr(resnames_arr)
        one_letter_seq = ''.join(one_letter_seq_arr)
        one_letter_arr_of_arr.append(one_letter_seq)
        bfactors_arr_of_arr.append(bfactors_arr)
        chains.append(chain.id)

dict_chain_oneletter_bfactor = {}
for chain_id, one_letter_seq, bfactor_arr in zip(chains, one_letter_arr_of_arr, bfactors_arr_of_arr):
    dict_chain_oneletter_bfactor[chain_id] = {
        'sequence': one_letter_seq,
        'bfactors': bfactor_arr
    }



#TODO: Function, dict
u = mda.Universe(pdb_cofactors)

# Organize cofactors by chain
cofactors_by_chain = {}
for residue in u.residues:
    # Get chain from first atom in residue
    chain_id = residue.atoms[0].chainID if hasattr(residue.atoms[0], 'chainID') else 'unknown'
    
    if chain_id not in cofactors_by_chain:
        cofactors_by_chain[chain_id] = {'labels': [], 'bfactors': []}
    
    label = str(residue.resid) + residue.resname
    bfactor = residue.atoms.tempfactors.mean()
    
    cofactors_by_chain[chain_id]['labels'].append(label)
    cofactors_by_chain[chain_id]['bfactors'].append(bfactor)
#TODO: Function, dict






# Function to plot a single sequence in a given axes
def plot_sequence(ax, one_letter_seq, lifetime_values, label, cmap, norm):
    """Plot a single protein sequence with B-factor coloring."""
    
    for i in range(-1, len(one_letter_seq)):
        if i == -1:
            # Draw empty padding box at the beginning
            rect_x = i * box_spacing - 0.8
            rect_y = sequence_y
            rect_width = box_width
            rect_height = box_height
            ax.add_patch(plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                      facecolor='white', edgecolor='none', alpha=0.0))
            continue

        # Get the letter and value for this position
        letter = one_letter_seq[i]
        value = lifetime_values[i]

        # Skip drawing if this is a padding character (dash)
        if letter == '-':
            continue

        # Convert RGBA tuple to hex color for better compatibility
        rgba = cmap(norm(value))
        hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

        # Draw fixed-size rectangle
        rect_x = i * box_spacing - 0.8  # Center the rectangle
        rect_y = sequence_y
        rect_width = box_width
        rect_height = box_height

        # Draw rectangle with color based on value
        ax.add_patch(plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                  facecolor=hex_color, edgecolor='black', alpha=0.8))

        # Place letter centered in the rectangle
        ax.text(i * box_spacing, rect_y + rect_height/2, letter, ha='center', va='center', 
                fontsize=seq_letter_fontsize, color=seq_letter_color, fontweight=seq_letter_fontweight)

    # Find the actual sequence length (excluding padding)
    actual_seq_length = len(one_letter_seq.rstrip('-'))
    
    # Add position labels above the boxes (positioned at resid_labels_y)
    # Label "1" above the first box
    ax.text(0 * box_spacing, resid_labels_y, '1', ha='center', va='bottom', 
            fontsize=resid_label_fontsize, color=resid_label_color, fontweight=resid_label_fontweight)

    # Add labels every 20 residues (1-indexed)
    for pos in range(20, actual_seq_length, 20):
        ax.text(pos * box_spacing, resid_labels_y, str(pos), ha='center', va='bottom', 
                fontsize=resid_label_fontsize, color=resid_label_color, fontweight=resid_label_fontweight)

    # Label sequence length above the last actual residue (not padding)
    ax.text((actual_seq_length - 1) * box_spacing, resid_labels_y, str(actual_seq_length), ha='center', va='bottom', 
            fontsize=resid_label_fontsize, color=resid_label_color, fontweight=resid_label_fontweight)

    # Add region rectangles above the sequence
    # Rectangle for residues 1-20 with label H1
    rect1_left = 0 * box_spacing - 0.8
    rect1_right = 19 * box_spacing + 0.8
    rect1_width = rect1_right - rect1_left
    rect1_center = (rect1_left + rect1_right) / 2
    ax.add_patch(plt.Rectangle((rect1_left, helix_boxes_y), rect1_width, 0.2,
                              facecolor=helix_facecolor, edgecolor=helix_edgecolor, alpha=helix_alpha))
    ax.text(rect1_center, helix_labels_y, 'H1', ha='center', va='center', 
            fontsize=helix_label_fontsize, color=helix_label_color, fontweight=helix_label_fontweight)

    # Rectangle for residues 100-150 with label "100 to 150"
    if len(one_letter_seq) > 150:  # Only draw if sequence is long enough 
        rect2_left = 99 * box_spacing - 0.8
        rect2_right = 149 * box_spacing + 0.8
        rect2_width = rect2_right - rect2_left
        rect2_center = (rect2_left + rect2_right) / 2
        ax.add_patch(plt.Rectangle((rect2_left, helix_boxes_y), rect2_width, 0.2,
                                  facecolor=helix_facecolor, edgecolor=helix_edgecolor, alpha=helix_alpha))
        ax.text(rect2_center, helix_labels_y, '100 to 150', ha='center', va='center', 
                fontsize=region_label_fontsize, color=region_label_color, fontweight=region_label_fontweight)

    # Set axis limits and styling
    ax.set_ylim(-0.5, plot_title_y + 0.5)
    ax.set_xlim(-2.0, (len(one_letter_seq) - 0.5) * box_spacing)
    ax.set_yticks([])
    ax.set_xticks([])
    # Do not show the axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Add plot title/label in top left
    ax.text(-1.5, plot_title_y, label, ha='left', va='top', 
            fontsize=plot_title_fontsize, color=plot_title_color, fontweight=plot_title_fontweight)

# Function to plot cofactors as a heatmap row
def plot_cofactors(ax, cofactor_labels, cofactor_bfactors, cmap, norm, max_seq_length):
    """Plot only non-zero cofactors with B-factor coloring below the sequence, within sequence bounds."""
    
    # Filter cofactors with non-zero B-factors
    active_cofactors = [(label, bfactor) for label, bfactor in zip(cofactor_labels, cofactor_bfactors) 
                        if bfactor > 0]
    
    if not active_cofactors:
        return  # No active cofactors to plot
    
    # Define y position for cofactors (below sequence)
    cofactor_y = 0.3  # Fixed y position for cofactor row
    
    # Limit number of cofactors to display to max_seq_length (to fit within sequence width)
    num_cofactors_to_plot = min(len(active_cofactors), max_seq_length)
    
    # Track x position as we place cofactors
    current_x = 0.0
    
    for i in range(num_cofactors_to_plot):
        cofactor_label, bfactor = active_cofactors[i]
        
        # Scale box width based on label length (number of letters)
        label_length = len(cofactor_label)
        cofactor_box_width = box_width * label_length
        
        # Position each cofactor at current_x (no longer at i * box_spacing)
        rect_x = current_x
        rect_y = cofactor_y
        rect_width = cofactor_box_width
        rect_height = box_height
        
        # Color based on B-factor value
        rgba = cmap(norm(bfactor))
        hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
        
        # Draw rectangle
        ax.add_patch(plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                  facecolor=hex_color, edgecolor='black', alpha=0.8))
        
        # Place cofactor label centered in the rectangle
        ax.text(current_x + cofactor_box_width/2, rect_y + rect_height/2, cofactor_label, ha='center', va='center', 
                fontsize=cofactor_fontsize, color=cofactor_color, fontweight=cofactor_fontweight)
        
        # Move to next position (add spacing between cofactors)
        current_x += cofactor_box_width + box_spacing
    
    # Update y-axis limits to include cofactors
    ax.set_ylim(-0.3, plot_title_y + 0.3)

# Get number of chains to plot
num_chains = len(dict_chain_oneletter_bfactor)
if num_chains == 0:
    raise ValueError("No chains found in the CIF file")

# Find the maximum sequence length for figure width calculation
max_seq_length = max(len(data['sequence']) for data in dict_chain_oneletter_bfactor.values())

# Pad sequences with trailing empty spaces to match max length
for chain_id in dict_chain_oneletter_bfactor:
    seq_len = len(dict_chain_oneletter_bfactor[chain_id]['sequence'])
    if seq_len < max_seq_length:
        # Pad sequence with spaces (represented as gaps)
        padding_needed = max_seq_length - seq_len
        dict_chain_oneletter_bfactor[chain_id]['sequence'] += '-' * padding_needed
        # Pad b-factors with zeros (invisible in plot)
        dict_chain_oneletter_bfactor[chain_id]['bfactors'] += [0] * padding_needed

total_width = (max_seq_length - 1) * box_spacing + box_width
figure_width = max(20, total_width * 0.1)  # Scale factor to make it reasonable
figure_height = 2.5 * num_chains  # Height per subplot

# Create figure with subplots
fig, axes = plt.subplots(num_chains, 1, figsize=(figure_width, figure_height), sharex=False)
if num_chains == 1:
    axes = [axes]  # Ensure axes is always a list

# Setup colormap and normalization (shared across all subplots)
cmap = getattr(plt.cm, cmap_name)
norm = plt.Normalize(vmin=vmin, vmax=vmax)

# Plot each chain in its subplot
for idx, (chain_id, chain_data) in enumerate(dict_chain_oneletter_bfactor.items()):
    one_letter_seq = chain_data['sequence']
    lifetime_values = chain_data['bfactors']
    label = chain_id_to_label[chain_id]
    
    plot_sequence(axes[idx], one_letter_seq, lifetime_values, label, cmap, norm)
    
    # Plot cofactors specific to this chain below the sequence
    if chain_id in cofactors_by_chain:
        chain_cofactors = cofactors_by_chain[chain_id]
        plot_cofactors(axes[idx], chain_cofactors['labels'], chain_cofactors['bfactors'], 
                       cmap, norm, max_seq_length)
    else:
        # No cofactors for this chain, just update ylim
        axes[idx].set_ylim(-0.3, plot_title_y + 0.3)

plt.tight_layout()
plt.savefig(output_figure, dpi=300, bbox_inches='tight')
print(f"Saved sequence plots for {num_chains} chains to {output_figure}")
