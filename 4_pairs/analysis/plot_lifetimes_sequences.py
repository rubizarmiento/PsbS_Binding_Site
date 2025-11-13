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
fasta_aligned = "/martini/rubiz/Github/PsbS_Binding_Site/fasta/5xnl_lhc_set_aligned.fasta"
output_figure = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/test.png"

# Read YAML file
with open(chain_labels_yaml, 'r') as f:
    data = yaml.safe_load(f)

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
cofactor_resnames = [residue.resname for residue in u.residues]
cofactor_bfactors = [residue.atoms.tempfactors.mean() for residue in u.residues]
#TODO: Function, dict

fasta_arr = []
fasta_labels = []
with open(fasta_aligned, 'r') as fasta_file:
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if sequence:
                fasta_arr.append(sequence)
            fasta_labels.append(line[1:])  # Remove '>' and keep the label
            sequence = ''
        else:
            sequence += line
    if sequence:
        fasta_arr.append(sequence)





# If CIF file contains protein, get the case
basename_cif = os.path.basename(cif_protein)
if basename_cif.endswith("_protein.cif"):
    type_plot = "one_letter"



# Test first sequence in dict CIF
first_chain_id = list(dict_chain_oneletter_bfactor.keys())[0]
one_letter_seq = dict_chain_oneletter_bfactor[first_chain_id]['sequence']
lifetime_values = dict_chain_oneletter_bfactor[first_chain_id]['bfactors']  # Placeholder for actual lifetime values

# Calculate figure width dynamically based on sequence length for consistent spacing
box_width = 1.6
spacing = 3.0  # Fixed spacing between box centers
total_width = (len(one_letter_seq) - 1) * spacing + box_width
figure_width = max(20, total_width * 0.1)  # Scale factor to make it reasonable
plt.figure(figsize=(figure_width, 2))

cmap = plt.cm.RdPu
norm = plt.Normalize(vmin=0, vmax=2000)

for i in range(-1, len(one_letter_seq)):
    if i == -1:
        # Draw empty padding box at the beginning
        rect_x = i * spacing - 0.8
        rect_y = 1.0
        rect_width = 1.6
        rect_height = 0.3
        plt.gca().add_patch(plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                         facecolor='white', edgecolor='none', alpha=0.0))
        continue

    # Get the letter and value for this position
    letter = one_letter_seq[i]
    value = lifetime_values[i]

    # Convert RGBA tuple to hex color for better compatibility
    rgba = cmap(norm(value))
    hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

    # Draw fixed-size rectangle
    rect_x = i * spacing - 0.8  # Center the rectangle
    rect_y = 1.0
    rect_width = 1.6
    rect_height = 0.3

    # Draw rectangle with color based on value
    plt.gca().add_patch(plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                     facecolor=hex_color, edgecolor='black', alpha=0.8))

    # Place letter centered in the rectangle
    plt.text(i * spacing, rect_y + rect_height/2, letter, ha='center', va='center', fontsize=10, color='black')

# Add position labels below the boxes
# Label "1" under the first box
plt.text(0 * spacing, 0.7, '1', ha='center', va='center', fontsize=8, color='black')

# Add labels every 20 residues (1-indexed)
for pos in range(20, len(one_letter_seq), 20):
    plt.text(pos * spacing, 0.7, str(pos), ha='center', va='center', fontsize=8, color='black')

# Label sequence length under the last box
seq_length = len(one_letter_seq)
plt.text((len(one_letter_seq) - 1) * spacing, 0.7, str(seq_length), ha='center', va='center', fontsize=8, color='black')

# Add region rectangles above the sequence
# Rectangle for residues 1-20 with label H1
rect1_left = 0 * spacing - 0.8
rect1_right = 19 * spacing + 0.8
rect1_width = rect1_right - rect1_left
rect1_center = (rect1_left + rect1_right) / 2
plt.gca().add_patch(plt.Rectangle((rect1_left, 1.5), rect1_width, 0.2,
                                 facecolor='lightblue', edgecolor='blue', alpha=0.5))
plt.text(rect1_center, 1.6, 'H1', ha='center', va='center', fontsize=10, color='blue', fontweight='bold')

# Rectangle for residues 100-150 with label "100 to 150"
if len(one_letter_seq) > 150:  # Only draw if sequence is long enough
    rect2_left = 99 * spacing - 0.8
    rect2_right = 149 * spacing + 0.8
    rect2_width = rect2_right - rect2_left
    rect2_center = (rect2_left + rect2_right) / 2
    plt.gca().add_patch(plt.Rectangle((rect2_left, 1.5), rect2_width, 0.2,
                                     facecolor='lightgreen', edgecolor='green', alpha=0.5))
    plt.text(rect2_center, 1.6, '100 to 150', ha='center', va='center', fontsize=10, color='green', fontweight='bold')

plt.ylim(0.5, 1.8)  # Adjusted to include space for region rectangles above boxes
plt.xlim(-2.0, (len(one_letter_seq) - 0.5) * spacing)  # Adjust x-axis limits for spacing with padding before first box
plt.yticks([])
plt.xticks([])
# Do not show the axes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)

# Add a colorbar to show the scale
# Note: Colorbar removed due to text-only plot - colors are visible in the letters
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = plt.colorbar(sm, orientation='vertical', shrink=0.8)
# cbar.set_label('Random Value (0-2000)')

plt.tight_layout()
plt.savefig(output_figure, dpi=300, bbox_inches='tight')
print(f"Saved sequence plot to {output_figure}")
