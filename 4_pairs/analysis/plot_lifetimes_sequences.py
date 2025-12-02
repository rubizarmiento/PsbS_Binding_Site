"""
Protein sequence visualization with B-factor coloring, cofactors, and helix annotations.

This module generates publication-quality sequence plots for PSII and PsbS proteins,
showing amino acid sequences with B-factor coloring, cofactor information organized
by chain, and helix region annotations dynamically loaded from YAML configurations.
"""
import yaml
import os   
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
from Bio.PDB import MMCIFParser
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

#Ignore warnings
import warnings
warnings.filterwarnings("ignore")

# ============================================================================
# Configuration Paths
# ============================================================================

CHAIN_LABELS_YAML = "/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_labels.yaml"
BASENAMES_CSV = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv"
CIFS_LIFETIMES_DIR = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes/"
COLOR_DEFINITIONS_YAML = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/color_definitions.yaml"
HELIX_DEFINITIONS_YAML_PSII = "/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix.yaml"
HELIX_DEFINITIONS_YAML_PSBS = "/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_helix.yaml"
OUTPUT_FIGURES_DIR = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/lifetimes/"

# ============================================================================
# Configuration Loading Functions
# ============================================================================

def load_yaml_file(file_path):
    """Load and parse a YAML file.
    
    Parameters
    ----------
    file_path : str
        Path to the YAML file to load.
    
    Returns
    -------
    dict
        Parsed YAML content as a dictionary.
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    yaml.YAMLError
        If the file is not valid YAML.
    """
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)


def get_config_value(config, path, default=None):
    """Extract nested value from config using dot-separated path.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary to extract from.
    path : str
        Dot-separated path to the value (e.g., 'text.sequence.letters.fontsize').
    default : any, optional
        Default value to return if path not found. Default is None.
    
    Returns
    -------
    any
        The value at the specified path, or default if not found.
    
    Examples
    --------
    >>> config = {'text': {'sequence': {'letters': {'fontsize': 12}}}}
    >>> get_config_value(config, 'text.sequence.letters.fontsize')
    12
    >>> get_config_value(config, 'missing.path', 'default_val')
    'default_val'
    """
    keys = path.split('.')
    value = config
    try:
        for key in keys:
            value = value[key]
        return value
    except (KeyError, TypeError):
        return default


def merge_helix_configs(psii_yaml_path, psbs_yaml_path):
    """Load and merge PSII and PsbS helix configuration files.
    
    Parameters
    ----------
    psii_yaml_path : str
        Path to PSII helix definitions YAML file.
    psbs_yaml_path : str
        Path to PsbS helix definitions YAML file.
    
    Returns
    -------
    dict
        Merged helix configuration with both PSII and PsbS chains.
    """
    helix_config_psii = load_yaml_file(psii_yaml_path)
    helix_config_psbs = load_yaml_file(psbs_yaml_path)
    
    helix_config = {}
    helix_config.update(helix_config_psii)
    helix_config.update(helix_config_psbs)
    
    return helix_config


def extract_config_variables(color_config):
    """Extract all configuration variables from color_config YAML.
    
    Parameters
    ----------
    color_config : dict
        Loaded color configuration dictionary.
    
    Returns
    -------
    dict
        Dictionary mapping variable names to their values from config.
    """
    config_mapping = {
        # Colormap
        'cmap_name': 'colormap.name',
        'vmin': 'colormap.normalization.vmin', 
        'vmax': 'colormap.normalization.vmax',
        
        # Layout
        'box_width': 'layout.box_width',
        'box_height': 'layout.box_height',
        'helix_box_height': 'layout.helix_box_height',
        'cofactor_box_width': 'layout.cofactor_box_width',
        'cofactor_box_height': 'layout.cofactor_box_height',
        'cofactor_box_spacing': 'layout.cofactor_box_spacing',
        'cofactor_box_padding': 'layout.cofactor_box_padding',
        'left_margin': 'layout.left_margin',
        'right_margin': 'layout.right_margin',
        'upper_margin': 'layout.upper_margin',
        'sequence_y': 'layout.positions.sequence_y',
        'resid_labels_y': 'layout.positions.resid_labels_y',
        'annotations_y': 'layout.positions.annotations_y',
        'helix_labels_y': 'layout.positions.helix_labels_y',
        'cofactor_y': 'layout.positions.cofactor_y',
        
        # Tight layout settings
        'tight_layout_enabled': 'layout.tight_layout.enabled',
        'tight_layout_pad': 'layout.tight_layout.pad',
        'tight_layout_w_pad': 'layout.tight_layout.w_pad',
        'tight_layout_h_pad': 'layout.tight_layout.h_pad',
        
        # Text styling - sequence letters
        'seq_letter_fontsize': 'text.sequence.letters.fontsize',
        'seq_letter_color': 'text.sequence.letters.color',
        'seq_letter_fontweight': 'text.sequence.letters.fontweight',
        
        # Box styling - sequence boxes
        'seq_box_edgecolor': 'text.sequence.boxes.edgecolor',
        'seq_box_edgewidth': 'text.sequence.boxes.edgewidth',
        'seq_box_alpha': 'text.sequence.boxes.alpha',
        
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
        
        # Text styling - plot title
        'plot_title_fontsize': 'text.annotations.plot_title.fontsize',
        'plot_title_color': 'text.annotations.plot_title.color',
        'plot_title_fontweight': 'text.annotations.plot_title.fontweight',
        
        # Box styling
        'helix_facecolor': 'boxes.helices.facecolor',
        'helix_edgecolor': 'boxes.helices.edgecolor',
        'helix_alpha': 'boxes.helices.alpha',
        'helix_side_padding': 'boxes.helices.side_padding'
    }
    
    return {var_name: get_config_value(color_config, path) 
            for var_name, path in config_mapping.items()}


# Load configuration files
chain_labels = load_yaml_file(CHAIN_LABELS_YAML)
color_config = load_yaml_file(COLOR_DEFINITIONS_YAML)
helix_config = merge_helix_configs(HELIX_DEFINITIONS_YAML_PSII, HELIX_DEFINITIONS_YAML_PSBS)

# Extract and set all configuration variables globally
config_vars = extract_config_variables(color_config)
globals().update(config_vars)

# Load tags from basenames CSV
basenames_df = pd.read_csv(BASENAMES_CSV, header=0, sep=' ')
try:
    all_tags = basenames_df['unique_basename'].values
except KeyError:
    print(f"Warning: Column 'unique_basename' not found in {BASENAMES_CSV}")
    print(f"         Falling back to first column as basenames")
    all_tags = basenames_df.iloc[:, 0].values

# Helper function to extract nested config values
def get_config_value(config, path, default=None):
    """Extract nested value from config using dot-separated path.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary to extract from.
    path : str
        Dot-separated path to the value (e.g., 'text.sequence.letters.fontsize').
    default : any, optional
        Default value to return if path not found. Default is None.
    
    Returns
    -------
    any
        The value at the specified path, or default if not found.
    
    Examples
    --------
    >>> config = {'text': {'sequence': {'letters': {'fontsize': 12}}}}
    >>> get_config_value(config, 'text.sequence.letters.fontsize')
    12
    >>> get_config_value(config, 'missing.path', 'default_val')
    'default_val'
    """
    keys = path.split('.')
    value = config
    try:
        for key in keys:
            value = value[key]
        return value
    except (KeyError, TypeError):
        return default




# ============================================================================
# Sequence Processing Functions
# ============================================================================

def resnames_to_one_letter(resnames):
    """Convert three-letter residue names to single-letter code.
    
    Parameters
    ----------
    resnames : list of str
        List of three-letter residue names (e.g., ['ALA', 'GLY', 'VAL']).
    
    Returns
    -------
    list of str
        List of single-letter amino acid codes (e.g., ['A', 'G', 'V']).
        Unknown residues are represented as 'X'.
    
    Examples
    --------
    >>> resnames_to_one_letter(['ALA', 'GLY', 'VAL'])
    ['A', 'G', 'V']
    """
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return [three_to_one.get(resname, 'X') for resname in resnames]


def extract_chain_sequences_from_cif(cif_path):
    """Extract protein sequences and B-factors from a CIF file.
    
    Parameters
    ----------
    cif_path : str
        Path to the CIF structure file.
    
    Returns
    -------
    dict
        Dictionary with keys as chain IDs and values as dicts containing:
        - 'sequence': str, single-letter amino acid codes
        - 'bfactors': list of float, B-factors per residue
    
    Examples
    --------
    >>> chains_data = extract_chain_sequences_from_cif('protein.cif')
    >>> chains_data['A']['sequence']
    'MKVLSPAD...'
    >>> len(chains_data['A']['bfactors'])
    350
    """
    parser = MMCIFParser()
    structure = parser.get_structure('protein', cif_path)
    
    chains_data = {}
    for model in structure:
        for chain in model:
            resnames_arr = []
            bfactors_arr = []
            
            for residue in chain:
                resnames_arr.append(residue.get_resname())
                bfactors = [atom.get_bfactor() for atom in residue]
                avg_bfactor = sum(bfactors) / len(bfactors) if bfactors else 0.0
                bfactors_arr.append(avg_bfactor)
            
            one_letter_seq = ''.join(resnames_to_one_letter(resnames_arr))
            chains_data[chain.id] = {
                'sequence': one_letter_seq,
                'bfactors': bfactors_arr
            }
    
    return chains_data


def split_psbs_sequence(sequence, bfactors, split_point=212):
    """Split PsbS sequence into two segments.
    
    Parameters
    ----------
    sequence : str
        Single-letter amino acid sequence.
    bfactors : list of float
        B-factor values for each residue.
    split_point : int, optional
        Position to split the sequence (default 212).
    
    Returns
    -------
    dict
        Dictionary with keys 'n_psbs_1' and 'n_psbs_2', each containing:
        - 'sequence': str, single-letter codes for that segment
        - 'bfactors': list of float, B-factors for that segment
    
    Examples
    --------
    >>> seq_data = split_psbs_sequence(psbs_seq, psbs_bfactors, split_point=212)
    >>> len(seq_data['n_psbs_1']['sequence'])
    212
    >>> len(seq_data['n_psbs_2']['sequence'])
    212
    """
    return {
        'n_psbs_1': {
            'sequence': sequence[:split_point],
            'bfactors': bfactors[:split_point]
        },
        'n_psbs_2': {
            'sequence': sequence[split_point:],
            'bfactors': bfactors[split_point:]
        }
    }


def extract_cofactors_by_chain(cif_path):
    """Extract cofactor information organized by chain from CIF file.
    
    Parameters
    ----------
    cif_path : str
        Path to the CIF file containing cofactors.
    
    Returns
    -------
    dict
        Dictionary with keys as chain IDs and values as dicts containing:
        - 'labels': list of str, cofactor labels (resid + resname)
        - 'bfactors': list of float, B-factors for each cofactor
    
    Examples
    --------
    >>> cofactors = extract_cofactors_by_chain('cofactors.cif')
    >>> cofactors['A']['labels']
    ['1CLH', '2HEM', ...]
    >>> len(cofactors['A']['bfactors'])
    2
    """
    parser = MMCIFParser()
    structure = parser.get_structure('cofactors', cif_path)
    
    cofactors_by_chain = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            if chain_id not in cofactors_by_chain:
                cofactors_by_chain[chain_id] = {'labels': [], 'bfactors': []}
            
            for residue in chain:
                # Get B-factor from first atom in the residue
                bfactors = [atom.get_bfactor() for atom in residue]
                avg_bfactor = sum(bfactors) / len(bfactors) if bfactors else 0.0
                
                label = str(residue.get_id()[1]) + residue.get_resname()
                cofactors_by_chain[chain_id]['labels'].append(label)
                cofactors_by_chain[chain_id]['bfactors'].append(avg_bfactor)
    
    return cofactors_by_chain


def build_chain_mapping(chains, chain_labels):
    """Build mapping of chain IDs to display labels.
    
    Parameters
    ----------
    chains : list of str
        Chain identifiers from CIF file, may include segment suffixes (e.g., 'c_seg1').
    chain_labels : dict
        Dictionary of known chain labels from YAML config.
    
    Returns
    -------
    dict
        Mapping of chain ID (including segment ID) to display label.
        Split segments use the same label as the base chain (without segment suffix).
    
    Examples
    --------
    >>> mapping = build_chain_mapping(['n', 's', 'c_seg1'], {'n': 'PSII-n'})
    >>> mapping['n']
    'PSII-n'
    >>> mapping['c_seg1']
    'PsbC'
    """
    chain_id_to_label = {}
    for chain_id in chains:
        # Check if this is a split segment
        if '_seg' in chain_id:
            # Extract base chain ID from segment (e.g., 'c_seg1' → 'c' or 'n_psbs_1_seg1' → 'n_psbs_1')
            base_id = chain_id.rsplit('_seg', 1)[0]
            
            # Check if this is a PsbS segment
            if base_id == 'n_psbs_1':
                base_label = "PsbS (Chain A)"
            elif base_id == 'n_psbs_2':
                base_label = "PsbS (Chain B)"
            else:
                # Get base label from chain_labels if available
                if base_id in chain_labels:
                    base_label = chain_labels[base_id]
                else:
                    # Generate label: capitalize chain ID (e.g., 'c' → 'PsbC')
                    base_label = f"Psb{base_id.upper()}"
            
            chain_id_to_label[chain_id] = base_label
        else:
            # Regular chain, no splitting
            if chain_id in chain_labels:
                chain_id_to_label[chain_id] = chain_labels[chain_id]
            else:
                # Generate label: capitalize chain ID (e.g., 'c' → 'PsbC')
                chain_id_to_label[chain_id] = f"Psb{chain_id.upper()}"
    
    # Add labels for PsbS (in case they were not added via split segments)
    chain_id_to_label['n_psbs_1'] = "PsbS (Chain A)"
    chain_id_to_label['n_psbs_2'] = "PsbS (Chain B)"
    
    return chain_id_to_label


def split_long_sequences(chains_data, split_threshold):
    """
    Split protein sequences longer than threshold into multiple segments.
    
    Parameters:
    -----------
    chains_data : dict
        Dictionary with chain IDs as keys and dicts containing 'sequence' and 'bfactors' as values
    split_threshold : int or None
        Maximum residues per segment. If None, no splitting is performed
        
    Returns:
    --------
    dict
        Modified chains_data with split sequences. Long sequences are split into
        segments with keys like "chain_A_seg1", "chain_A_seg2", etc.
        Each segment contains the split sequence and corresponding bfactors.
    """
    if split_threshold is None or split_threshold <= 0:
        return chains_data
    
    split_chains_data = {}
    
    for chain_id, chain_data in chains_data.items():
        seq = chain_data['sequence']
        bfactors = chain_data['bfactors']
        
        if len(seq) <= split_threshold:
            # Sequence is short enough, keep as is
            split_chains_data[chain_id] = chain_data
        else:
            # Split into segments
            num_segments = (len(seq) + split_threshold - 1) // split_threshold
            for seg_idx in range(num_segments):
                start_idx = seg_idx * split_threshold
                end_idx = min((seg_idx + 1) * split_threshold, len(seq))
                
                seg_seq = seq[start_idx:end_idx]
                seg_bfactors = bfactors[start_idx:end_idx]
                seg_key = f"{chain_id}_seg{seg_idx + 1}"
                
                split_chains_data[seg_key] = {
                    'sequence': seg_seq,
                    'bfactors': seg_bfactors
                }
    
    return split_chains_data




# Function to plot a single sequence in a given axes
def plot_sequence(ax, one_letter_seq, lifetime_values, label, cmap, norm, helix_definitions=None, residue_offset=0, show_label=True, debug=False):
    """Plot a single protein sequence with B-factor coloring.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    one_letter_seq : str
        Single-letter amino acid sequence.
    lifetime_values : list of float
        B-factor values for each residue.
    label : str
        Label to display for this sequence.
    cmap : matplotlib.cm.Colormap
        Colormap for coloring.
    norm : matplotlib.colors.Normalize
        Normalization for color mapping.
    helix_definitions : dict, optional
        Helix definitions with 'start' and 'end' keys.
    residue_offset : int, optional
        Offset to add to residue numbering (for split segments). Default: 0.
    show_label : bool, optional
        Whether to display the label (True for first segment only). Default: True.
    debug : bool, optional
        Whether to show debug visualization. Default: False.
    """
    
    # Use fixed box dimensions (centered on letter position)
    seq_box_width = globals()['box_width']
    seq_box_height = globals()['box_height']
    seq_letter_fontsize = globals()['seq_letter_fontsize']
    seq_letter_color = globals()['seq_letter_color']
    seq_letter_fontweight = globals()['seq_letter_fontweight']
    seq_box_edgecolor = globals()['seq_box_edgecolor']
    seq_box_edgewidth = globals()['seq_box_edgewidth']
    seq_box_alpha = globals()['seq_box_alpha']
    sequence_y = globals()['sequence_y']
    resid_labels_y = globals()['resid_labels_y']
    resid_label_fontsize = globals()['resid_label_fontsize']
    resid_label_color = globals()['resid_label_color']
    resid_label_fontweight = globals()['resid_label_fontweight']
    helix_side_padding = globals()['helix_side_padding']
    helix_label_fontsize = globals()['helix_label_fontsize']
    helix_labels_y = globals()['helix_labels_y']
    helix_facecolor = globals()['helix_facecolor']
    helix_edgecolor = globals()['helix_edgecolor']
    helix_alpha = globals()['helix_alpha']
    helix_label_color = globals()['helix_label_color']
    helix_label_fontweight = globals()['helix_label_fontweight']
    annotations_y = globals()['annotations_y']
    upper_margin = globals()['upper_margin']
    right_margin = globals()['right_margin']
    left_margin = globals()['left_margin']
    plot_title_fontsize = globals()['plot_title_fontsize']
    plot_title_color = globals()['plot_title_color']
    plot_title_fontweight = globals()['plot_title_fontweight']
    
    for i in range(-1, len(one_letter_seq)):
        if i == -1:
            # Draw empty padding box at the beginning
            rect_x = i * seq_box_width - seq_box_width / 2
            rect_y = sequence_y - seq_box_height / 2
            ax.add_patch(plt.Rectangle((rect_x, rect_y), seq_box_width, seq_box_height,
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

        # Draw fixed-size rectangle (aligned with text position at i * seq_box_width)
        rect_x = i * seq_box_width - seq_box_width / 2  # Position box centered on text
        rect_y = sequence_y - seq_box_height / 2
        
        # Draw rectangle with color based on value
        ax.add_patch(plt.Rectangle((rect_x, rect_y), seq_box_width, seq_box_height,
                                  facecolor=hex_color, edgecolor=seq_box_edgecolor, 
                                  linewidth=seq_box_edgewidth, alpha=seq_box_alpha))

        # Place letter centered in the rectangle
        ax.text(i * seq_box_width, sequence_y, letter, ha='center', va='center', 
                fontsize=seq_letter_fontsize, color=seq_letter_color, fontweight=seq_letter_fontweight)

    # Find the actual sequence length (excluding padding)
    actual_seq_length = len(one_letter_seq.rstrip('-'))
    
    # Determine if this is a split segment
    is_split_segment = '_seg' in one_letter_seq  # This will be False, need to pass it as param
    # For now, check if chain_id has _seg (we need access to it)
    # We'll determine this by checking if residue_offset > 0
    is_intermediate_segment = residue_offset > 0
    
    # Add position labels above the boxes (positioned at resid_labels_y)
    # Label first residue number (with offset) - always show
    first_resid = 1 + residue_offset
    ax.text(0 * seq_box_width, resid_labels_y, str(first_resid), ha='center', va='bottom', 
            fontsize=resid_label_fontsize, color=resid_label_color, fontweight=resid_label_fontweight)

    # Add labels every 20 residues (at positions 19, 39, 59, etc. which correspond to residues 20, 40, 60, etc.)
    labeled_positions = {0}  # Track which positions we've labeled
    for pos in range(19, actual_seq_length, 20):
        resid = pos + 1 + residue_offset
        ax.text(pos * seq_box_width, resid_labels_y, str(resid), ha='center', va='bottom', 
                fontsize=resid_label_fontsize, color=resid_label_color, fontweight=resid_label_fontweight)
        labeled_positions.add(pos)

    # Label sequence length above the last actual residue (not padding)
    # Only show on the final segment of a chain (not on intermediate split segments)
    final_pos_in_segment = actual_seq_length - 1
    if final_pos_in_segment not in labeled_positions and final_pos_in_segment > 0 and not is_intermediate_segment:
        final_resid = actual_seq_length + residue_offset
        ax.text(final_pos_in_segment * seq_box_width, resid_labels_y, str(final_resid), ha='center', va='bottom', 
                fontsize=resid_label_fontsize, color=resid_label_color, fontweight=resid_label_fontweight)

    # Add region rectangles above the sequence
    # Draw helices from definitions if available
    if helix_definitions:
        for helix_id, helix_info in helix_definitions.items():
            start_res = helix_info['start']
            end_res = helix_info['end']
            
            # Adjust for segment offset - only draw if helix is within this segment
            start_in_segment = start_res - residue_offset
            end_in_segment = end_res - residue_offset
            
            # Check if helix overlaps with this segment
            if end_in_segment > 0 and start_in_segment <= actual_seq_length:
                # Clamp to segment bounds
                start_in_segment = max(1, start_in_segment)
                end_in_segment = min(actual_seq_length, end_in_segment)
                
                # Get box dimensions from config
                box_width_val = globals()['box_width']
                
                # Calculate helix width based on number of residues
                # Helix width = num_residues * box_width (exact span, no padding)
                num_residues_in_helix = end_in_segment - start_in_segment + 1
                helix_width = num_residues_in_helix * box_width_val
                
                # Calculate left position: first residue box left edge
                first_residue_center = (start_in_segment - 1) * seq_box_width
                helix_left = first_residue_center - box_width_val / 2
                # Center is at the middle of the helix width
                helix_center = helix_left + helix_width / 2
                
                # Use fixed helix box height (from config)
                helix_box_height_val = globals()['helix_box_height']
                
                # Center helix box at helix_labels_y (box center aligned with label anchor point)
                helix_box_bottom = helix_labels_y - helix_box_height_val / 2
                ax.add_patch(plt.Rectangle((helix_left, helix_box_bottom), helix_width, helix_box_height_val,
                                          facecolor=helix_facecolor, edgecolor=helix_edgecolor, alpha=helix_alpha))
                ax.text(helix_center, helix_labels_y, helix_id, ha='center', va='center', 
                        fontsize=helix_label_fontsize, color=helix_label_color, fontweight=helix_label_fontweight)

    # Set axis limits and styling
    ax.set_ylim(0.0, annotations_y + upper_margin)
    # Right limit: account for the last box width (seq_box_width/2 on each side) + right margin
    right_xlim = (len(one_letter_seq) - 1) * seq_box_width + seq_box_width / 2 + right_margin
    ax.set_xlim(-left_margin, right_xlim)
    ax.set_yticks([])
    ax.set_xticks([])
    # Do not show the axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Add plot title/label in top left (only if show_label is True)
    if show_label:
        # Position label at the left edge of the first sequence box
        label_x = -seq_box_width / 2
        ax.text(label_x, helix_labels_y, label, ha='left', va='center', 
                fontsize=plot_title_fontsize, color=plot_title_color, fontweight=plot_title_fontweight)
    
    # DEBUG: Draw subplot bounds (red for sequences) - only if debug is enabled
    if debug:
        plot_subplot_limits_debug(ax, 'sequence', alpha=0.1)

# Function to plot cofactors as a heatmap row
def plot_cofactors(ax, cofactor_labels, cofactor_bfactors, cmap, norm, max_seq_length, debug=False):
    """Plot only non-zero cofactors with B-factor coloring below the sequence, within sequence bounds."""
    
    # Extract config variables from globals
    cofactor_y = globals()['cofactor_y']
    cofactor_fontsize = globals()['cofactor_fontsize']
    cofactor_color = globals()['cofactor_color']
    cofactor_fontweight = globals()['cofactor_fontweight']
    cofactor_box_padding = globals()['cofactor_box_padding']
    cofactor_box_spacing = globals()['cofactor_box_spacing']
    seq_box_width = globals()['box_width']
    annotations_y = globals()['annotations_y']
    upper_margin = globals()['upper_margin']
    left_margin = globals()['left_margin']
    right_margin = globals()['right_margin']
    
    # Filter cofactors with non-zero B-factors
    active_cofactors = [(label, bfactor) for label, bfactor in zip(cofactor_labels, cofactor_bfactors) 
                        if bfactor > 0]
    
    if not active_cofactors:
        return  # No active cofactors to plot
    
    # Get y position for cofactors from configuration
    cofactor_y_pos = cofactor_y  # Use global config variable (from color_definitions.yaml)
    
    # Use fixed cofactor box dimensions from configuration
    cofactor_box_height_val = globals()['cofactor_box_height']
    cofactor_box_width_val = globals()['cofactor_box_width']
    
    # Limit number of cofactors to display to max_seq_length (to fit within sequence width)
    num_cofactors_to_plot = min(len(active_cofactors), max_seq_length)
    
    # Track x position as we place cofactors
    current_x = 0.0
    
    for i in range(num_cofactors_to_plot):
        cofactor_label, bfactor = active_cofactors[i]
        
        # Use fixed box width from configuration
        cofactor_box_width = cofactor_box_width_val
        
        # Position each cofactor at current_x (no longer at i * box_spacing)
        rect_x = current_x
        rect_y = cofactor_y_pos - cofactor_box_height_val / 2
        rect_width = cofactor_box_width
        rect_height = cofactor_box_height_val
        
        # Color based on B-factor value
        rgba = cmap(norm(bfactor))
        hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
        
        # Draw rectangle
        ax.add_patch(plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                  facecolor=hex_color, edgecolor='black', alpha=0.8))
        
        # Place cofactor label centered in the rectangle
        ax.text(current_x + cofactor_box_width/2, cofactor_y_pos, cofactor_label, ha='center', va='center', 
                fontsize=cofactor_fontsize, color=cofactor_color, fontweight=cofactor_fontweight)
        
        # Move to next position (add cofactor spacing between boxes)
        current_x += cofactor_box_width + cofactor_box_spacing
    
    # Update axis limits and styling to match sequence plots
    ax.set_ylim(0.0, annotations_y + upper_margin)
    # Set x-axis limits to match the sequence plot dimensions (boxes are adjacent with no spacing)
    right_xlim = (max_seq_length - 1) * seq_box_width + seq_box_width / 2 + right_margin
    ax.set_xlim(-left_margin, right_xlim)
    ax.set_yticks([])
    ax.set_xticks([])
    # Hide all axis spines like in sequence plots
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # DEBUG: Draw subplot bounds (blue for cofactors) - only if debug is enabled
    if debug:
        plot_subplot_limits_debug(ax, 'cofactor', alpha=0.1)


def plot_subplot_limits_debug(ax, subplot_type, alpha=0.1):
    """
    Plot the subplot limits as a debug visualization.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to debug.
    subplot_type : str
        Type of subplot: 'sequence' (red) or 'cofactor' (blue).
    alpha : float, optional
        Transparency of the debug box (default: 0.1).
    """
    # Get current limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Choose color based on type
    color = 'red' if subplot_type == 'sequence' else 'blue' if subplot_type == 'cofactor' else 'green'
    
    # Draw a rectangle showing the subplot bounds
    rect = plt.Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0],
                         facecolor=color, alpha=alpha, edgecolor=color, linewidth=2, zorder=-1)
    ax.add_patch(rect)


# ============================================================================
# Plotting Support Functions
# ============================================================================

def get_helix_key_for_chain(chain_id):
    """Map chain ID to its helix definition key.
    
    Parameters
    ----------
    chain_id : str
        Chain identifier from CIF file, may include segment suffix.
    
    Returns
    -------
    str
        Key to use for looking up helix definitions.
        PsbS segments are mapped to chain_A and chain_B.
        Split segments use base chain ID (e.g., 'c_seg1' → 'chain_c').
    
    Examples
    --------
    >>> get_helix_key_for_chain('n')
    'chain_n'
    >>> get_helix_key_for_chain('n_psbs_1')
    'chain_A'
    >>> get_helix_key_for_chain('c_seg1')
    'chain_c'
    >>> get_helix_key_for_chain('n_psbs_1_seg1')
    'chain_A'
    """
    # Remove segment suffix if present (e.g., 'n_psbs_1_seg1' → 'n_psbs_1')
    base_chain_id = chain_id.rsplit('_seg', 1)[0] if '_seg' in chain_id else chain_id
    
    # Check for PsbS chains
    if base_chain_id == 'n_psbs_1':
        return 'chain_A'
    elif base_chain_id == 'n_psbs_2':
        return 'chain_B'
    else:
        return f"chain_{base_chain_id}"


def get_residue_offset_for_chain(chain_id, original_chains_data, split_threshold):
    """Calculate the residue offset for a split segment.
    
    Parameters
    ----------
    chain_id : str
        Current chain ID (may include segment number).
    original_chains_data : dict
        Original chains_data before splitting (to get original sequence lengths).
    split_threshold : int or None
        The split threshold used (None if no splitting).
    
    Returns
    -------
    int
        Residue offset for this chain/segment (0 for non-split chains).
    """
    if '_seg' not in chain_id or split_threshold is None:
        return 0
    
    # Extract base chain ID and segment number
    base_id, seg_part = chain_id.rsplit('_seg', 1)
    seg_num = int(seg_part)
    
    # Calculate offset: (seg_num - 1) * split_threshold
    return (seg_num - 1) * split_threshold


def pad_sequences_to_max_length(chains_data):
    """Pad all sequences to the same maximum length.
    
    Parameters
    ----------
    chains_data : dict
        Dictionary of chain sequences and B-factors.
        Modified in-place to add padding.
    
    Returns
    -------
    int
        The maximum sequence length after padding.
    """
    max_seq_length = max(len(data['sequence']) for data in chains_data.values())
    
    for chain_id in chains_data:
        seq_len = len(chains_data[chain_id]['sequence'])
        if seq_len < max_seq_length:
            padding_needed = max_seq_length - seq_len
            chains_data[chain_id]['sequence'] += '-' * padding_needed
            chains_data[chain_id]['bfactors'] += [0] * padding_needed
    
    return max_seq_length


def create_figure_and_axes(num_chains, max_seq_length, box_width, num_cofactor_chains=0, height_ratios=None):
    """Create matplotlib figure and axes for multi-chain plotting.
    
    Parameters
    ----------
    num_chains : int
        Number of chains to plot.
    max_seq_length : int
        Maximum length of sequences.
    box_width : float
        Width of each sequence box (boxes are adjacent with no spacing).
    num_cofactor_chains : int, optional
        Number of chains with cofactors (default: 0, will be estimated as num_chains/2).
    height_ratios : list of float, optional
        Height ratios for each subplot. If None, all subplots have equal height.
    
    Returns
    -------
    tuple
        (fig, axes) - matplotlib Figure and Axes objects.
    """
    total_width = (max_seq_length - 1) * box_width + box_width
    figure_width = max(20, total_width * 0.1)
    
    # Estimate number of axes: each chain + estimate of cofactor chains
    if num_cofactor_chains == 0:
        num_cofactor_chains = max(1, num_chains // 2)  # Rough estimate
    
    num_axes = num_chains + num_cofactor_chains
    figure_height = 2.5 * num_axes
    
    # Create axes for sequences and cofactors with optional height ratios
    if height_ratios is not None:
        gridspec_kw = {'height_ratios': height_ratios}
        fig, axes = plt.subplots(num_axes, 1, figsize=(figure_width, figure_height), 
                                sharex=False, gridspec_kw=gridspec_kw)
    else:
        fig, axes = plt.subplots(num_axes, 1, figsize=(figure_width, figure_height), sharex=False)
    
    if num_axes == 1:
        axes = [axes]
    
    return fig, axes


def calculate_subplot_heights_and_count(chains_data, cofactors_data, box_height, helix_box_height, 
                                        resid_label_fontsize, cofactor_fontsize, 
                                        split_threshold=None):
    """Calculate subplot heights dynamically based on box sizes and font sizes.
    
    Parameters
    ----------
    chains_data : dict
        Dictionary of all chains with their sequences.
    cofactors_data : dict
        Dictionary of cofactors by base chain ID.
    box_height : float
        Height of sequence boxes (in data coordinates, typically 0.5).
    helix_box_height : float
        Height of helix boxes (in data coordinates, typically 0.5).
    resid_label_fontsize : float
        Font size for residue labels (in points).
    cofactor_fontsize : float
        Font size for cofactor labels (in points).
    split_threshold : int or None, optional
        The split threshold used for segments.
    
    Returns
    -------
    tuple of (list of float, int, float)
        Height ratios for GridSpec, exact number of axes needed, and total figure height in inches.
    """
    # Calculate sequence subplot height based on content:
    # - Helix boxes: ~0.3-0.4 inches
    # - Residue labels (10pt font): ~0.2 inches
    # - Sequence boxes (10pt font): ~0.2 inches
    # - Padding and spacing: ~0.3 inches
    # Total: ~1.0-1.2 inches per sequence
    sequence_height = 1.1
    
    # Calculate cofactor subplot height based on content:
    # - Cofactor labels (10pt font): ~0.2 inches
    # - Padding: ~0.2 inches
    # Total: ~0.4 inches per cofactor
    cofactor_height = 0.4
    
    heights = []
    previous_base_chain_id = None
    
    for chain_id in chains_data.keys():
        # Extract base chain ID
        base_chain_id = chain_id.split('_seg')[0] if '_seg' in chain_id else chain_id
        
        # If we've moved to a new base chain, add cofactor height for the previous chain
        if previous_base_chain_id is not None and base_chain_id != previous_base_chain_id:
            if previous_base_chain_id in cofactors_data:
                heights.append(cofactor_height)
        
        # Add sequence height for current chain
        heights.append(sequence_height)
        
        previous_base_chain_id = base_chain_id
    
    # Add cofactor height for the last chain
    if previous_base_chain_id is not None and previous_base_chain_id in cofactors_data:
        heights.append(cofactor_height)
    
    num_axes = len(heights)
    total_height = sum(heights)  # No extra margin, use exact sum of subplot heights
    
    # Convert absolute heights to ratios for GridSpec
    height_ratios = [h / sum(heights) for h in heights]
    
    return height_ratios, num_axes, total_height


def plot_all_chains(axes, chains_data, chain_id_to_label, cofactors_data, 
                    helix_config, cmap, norm, max_seq_length, split_threshold=None, debug=False):
    """Plot all chains on their respective axes with configurable spacing.
    
    Parameters
    ----------
    axes : list of matplotlib.axes.Axes
        List of axes for each chain.
    chains_data : dict
        Sequences and B-factors for all chains.
    chain_id_to_label : dict
        Mapping of chain IDs to display labels.
    cofactors_data : dict
        Cofactor information by chain.
    helix_config : dict
        Helix definitions by chain.
    cmap : matplotlib.cm.Colormap
        Colormap for B-factor visualization.
    norm : matplotlib.colors.Normalize
        Normalization for B-factor values.
    max_seq_length : int
        Maximum sequence length.
    split_threshold : int or None, optional
        The split threshold used for segments.
    debug : bool, optional
        Whether to show debug visualization. Default: False.
    
    Returns
    -------
    int
        Number of axes actually used.
    """
    axis_idx = 0
    previous_base_chain_id = None
    
    for chain_id, chain_data in chains_data.items():
        one_letter_seq = chain_data['sequence']
        lifetime_values = chain_data['bfactors']
        label = chain_id_to_label[chain_id]
        
        # Calculate residue offset
        residue_offset = get_residue_offset_for_chain(chain_id, chains_data, split_threshold)
        show_label = (residue_offset == 0)
        
        helix_key = get_helix_key_for_chain(chain_id)
        chain_helices = helix_config.get(helix_key, {})
        
        # Extract base chain ID
        base_chain_id = chain_id.split('_seg')[0] if '_seg' in chain_id else chain_id
        
        # If we've moved to a new base chain, plot cofactors for the previous chain
        if previous_base_chain_id is not None and base_chain_id != previous_base_chain_id:
            if previous_base_chain_id in cofactors_data:
                chain_cofactors = cofactors_data[previous_base_chain_id]
                plot_cofactors(axes[axis_idx], chain_cofactors['labels'], chain_cofactors['bfactors'],
                              cmap, norm, max_seq_length, debug=debug)
                axis_idx += 1
        
        # Plot sequence
        plot_sequence(axes[axis_idx], one_letter_seq, lifetime_values, label, cmap, norm,
                     helix_definitions=chain_helices, residue_offset=residue_offset, show_label=show_label, debug=debug)
        axis_idx += 1
        
        previous_base_chain_id = base_chain_id
    
    # Plot cofactors for the last chain
    if previous_base_chain_id is not None and previous_base_chain_id in cofactors_data:
        chain_cofactors = cofactors_data[previous_base_chain_id]
        plot_cofactors(axes[axis_idx], chain_cofactors['labels'], chain_cofactors['bfactors'],
                      cmap, norm, max_seq_length, debug=debug)
        axis_idx += 1
    
    return axis_idx


# ============================================================================
# Main Execution Loop
# ============================================================================

def create_argument_parser():
    """Create and return the command-line argument parser.
    
    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser with all available options.
    
    Examples
    --------
    >>> parser = create_argument_parser()
    >>> args = parser.parse_args(['-o', '/custom/path'])
    >>> args.output_dir
    '/custom/path'
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Generate protein sequence plots with B-factor coloring, cofactors, and helix annotations.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default configuration
  python plot_lifetimes_sequences.py
  
  # Use custom color configuration
  python plot_lifetimes_sequences.py -c /path/to/colors.yaml
  
  # Use custom output directory
  python plot_lifetimes_sequences.py -o /custom/output/path
  
  # Combine multiple options
  python plot_lifetimes_sequences.py \\
    -d /data/cifs_lifetimes/ \\
    -o /results/plots/ \\
    -c /configs/colors.yaml
        """
    )
    
    parser.add_argument(
        '-d', '--cifs-dir',
        dest='cifs_dir',
        default=CIFS_LIFETIMES_DIR,
        help=f'Directory containing CIF files (default: {CIFS_LIFETIMES_DIR})'
    )
    
    parser.add_argument(
        '-b', '--basenames-csv',
        dest='basenames_csv',
        default=BASENAMES_CSV,
        help=f'Path to basenames CSV file (default: {BASENAMES_CSV})'
    )
    
    parser.add_argument(
        '-l', '--chain-labels',
        dest='chain_labels',
        default=CHAIN_LABELS_YAML,
        help=f'Path to chain labels YAML file (default: {CHAIN_LABELS_YAML})'
    )
    
    parser.add_argument(
        '-c', '--color-config',
        dest='color_config',
        default=COLOR_DEFINITIONS_YAML,
        help=f'Path to color definitions YAML (default: {COLOR_DEFINITIONS_YAML})'
    )
    
    parser.add_argument(
        '-p', '--helix-psii',
        dest='helix_psii',
        default=HELIX_DEFINITIONS_YAML_PSII,
        help=f'Path to PSII helix definitions YAML (default: {HELIX_DEFINITIONS_YAML_PSII})'
    )
    
    parser.add_argument(
        '-s', '--helix-psbs',
        dest='helix_psbs',
        default=HELIX_DEFINITIONS_YAML_PSBS,
        help=f'Path to PsbS helix definitions YAML (default: {HELIX_DEFINITIONS_YAML_PSBS})'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        dest='output_dir',
        default=OUTPUT_FIGURES_DIR,
        help=f'Output directory for generated plots (default: {OUTPUT_FIGURES_DIR})'
    )
    
    parser.add_argument(
        '--split-sequences',
        dest='split_sequences',
        type=int,
        default=None,
        help='Split sequences longer than this value across multiple rows (e.g., 250). Default: None (no splitting)'
    )
    
    parser.add_argument(
        '--debug',
        dest='debug',
        action='store_true',
        default=False,
        help='Enable debug visualization (shows subplot bounds). Default: False'
    )
    
    return parser


def main():
    """Main function to process all cases and generate plots."""
    # Parse command-line arguments
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Update global paths with command-line arguments
    global CIFS_LIFETIMES_DIR, BASENAMES_CSV, CHAIN_LABELS_YAML
    global COLOR_DEFINITIONS_YAML, HELIX_DEFINITIONS_YAML_PSII
    global HELIX_DEFINITIONS_YAML_PSBS, OUTPUT_FIGURES_DIR
    
    CIFS_LIFETIMES_DIR = args.cifs_dir
    BASENAMES_CSV = args.basenames_csv
    CHAIN_LABELS_YAML = args.chain_labels
    COLOR_DEFINITIONS_YAML = args.color_config
    HELIX_DEFINITIONS_YAML_PSII = args.helix_psii
    HELIX_DEFINITIONS_YAML_PSBS = args.helix_psbs
    OUTPUT_FIGURES_DIR = args.output_dir
    
    # Reload configurations with updated paths
    global chain_labels, color_config, helix_config, all_tags
    global config_vars
    
    chain_labels = load_yaml_file(CHAIN_LABELS_YAML)
    color_config = load_yaml_file(COLOR_DEFINITIONS_YAML)
    helix_config = merge_helix_configs(HELIX_DEFINITIONS_YAML_PSII, HELIX_DEFINITIONS_YAML_PSBS)
    
    # Extract and set all configuration variables globally
    config_vars = extract_config_variables(color_config)
    globals().update(config_vars)
    
    # NOTE: Automatic rescaling of box_spacing disabled to allow direct control via config
    # If you want automatic rescaling based on sequence length, enable the code below:
    # max_expected_seq_length = 550  # Maximum expected sequence length
    # max_data_units = 100  # Target max data coordinate range
    # seq_box_width_expected = globals()['box_width']
    # expected_width = (max_expected_seq_length - 1) * globals()['box_spacing'] + seq_box_width_expected
    # if expected_width > max_data_units:
    #     scale_factor = max_data_units / expected_width
    #     globals()['box_spacing'] = globals()['box_spacing'] * scale_factor
    
    # Load tags from basenames CSV
    basenames_df = pd.read_csv(BASENAMES_CSV, header=0, sep=' ')
    try:
        all_tags = basenames_df['unique_basename'].values
    except KeyError:
        print(f"Warning: Column 'unique_basename' not found in {BASENAMES_CSV}")
        print(f"         Falling back to first column as basenames")
        basenames_df = pd.read_csv(BASENAMES_CSV, header=None, sep=' ')
        all_tags = basenames_df.iloc[:, 0].values
    
    # Create output directory if it doesn't exist
    Path(OUTPUT_FIGURES_DIR).mkdir(parents=True, exist_ok=True)
    
    # Setup colormap and normalization (shared across all plots)
    cmap_name = globals()['cmap_name']
    vmin = globals()['vmin']
    vmax = globals()['vmax']
    cmap = getattr(plt.cm, cmap_name)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    
    # Iterate over all tags/cases
    for current_case in all_tags:
        print(f"Processing case: {current_case}")
        
        # Build file paths for this case
        cif_protein = f"{CIFS_LIFETIMES_DIR}/{current_case}_protein.cif"
        cif_psbs = f"{CIFS_LIFETIMES_DIR}/{current_case}_psbs.cif"
        cif_cofactors = f"{CIFS_LIFETIMES_DIR}/{current_case}_cofactors.cif"
        output_figure = f"{OUTPUT_FIGURES_DIR}/{current_case}.png"
        
        # Check if files exist
        if not os.path.exists(cif_protein):
            print(f"  Warning: {cif_protein} not found, skipping")
            continue
        
        try:
            # Extract chain sequences from protein CIF
            chains_data = extract_chain_sequences_from_cif(cif_protein)
            chains_list = list(chains_data.keys())
            
            # Extract PsbS sequences and split
            psbs_chains = extract_chain_sequences_from_cif(cif_psbs)
            for chain_id, chain_data in psbs_chains.items():
                psbs_segments = split_psbs_sequence(chain_data['sequence'], chain_data['bfactors'])
                chains_data.update(psbs_segments)
            
            # Apply sequence splitting if requested
            if args.split_sequences:
                chains_data = split_long_sequences(chains_data, args.split_sequences)
            
            # Extract cofactors
            cofactors_data = extract_cofactors_by_chain(cif_cofactors)
            
            # Build chain labels (using current keys in case sequences were split)
            chain_id_to_label = build_chain_mapping(list(chains_data.keys()), chain_labels)
            
            # Validate chains
            num_chains = len(chains_data)
            if num_chains == 0:
                print(f"  Warning: No chains found in {cif_protein}")
                continue
            
            # Pad sequences and get max length
            max_seq_length = pad_sequences_to_max_length(chains_data)
            
            # Use fixed box dimensions (all in globals after config update)
            seq_box_width = globals()['box_width']
            seq_box_height = globals()['box_height']
            helix_box_height = globals()['helix_box_height']
            
            # Calculate figure width based on sequence length (boxes are adjacent with no spacing)
            total_width = (max_seq_length - 1) * seq_box_width + seq_box_width
            figure_width = max(20, total_width * 0.1)
            
            # Calculate subplot heights dynamically based on box sizes and font sizes
            height_ratios, num_axes_exact, figure_height = calculate_subplot_heights_and_count(
                chains_data, cofactors_data, 
                seq_box_height, helix_box_height,
                globals()['resid_label_fontsize'], globals()['cofactor_fontsize'],
                args.split_sequences)
            
            # Create figure and axes with height ratios
            fig, axes = plt.subplots(num_axes_exact, 1, figsize=(figure_width, figure_height), 
                                    sharex=False, gridspec_kw={'height_ratios': height_ratios})
            if num_axes_exact == 1:
                axes = [axes]
            
            # Plot all chains and get number of axes used
            axes_used = plot_all_chains(axes, chains_data, chain_id_to_label, cofactors_data,
                          helix_config, cmap, norm, max_seq_length, args.split_sequences, debug=args.debug)
            
            # Remove unused axes
            for unused_ax in axes[axes_used:]:
                fig.delaxes(unused_ax)
            
            # Apply tight layout with configurable parameters
            if globals()['tight_layout_enabled']:
                plt.tight_layout(pad=globals()['tight_layout_pad'], w_pad=globals()['tight_layout_w_pad'], h_pad=globals()['tight_layout_h_pad'])
            else:
                plt.tight_layout()
            
            # Save figure
            plt.savefig(output_figure, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            print(f"  ✓ Saved sequence plots for {num_chains} chains to {output_figure}")
            
        except Exception as e:
            print(f"  Error processing {current_case}: {e}")
            plt.close('all')
            continue


if __name__ == '__main__':
    main()

