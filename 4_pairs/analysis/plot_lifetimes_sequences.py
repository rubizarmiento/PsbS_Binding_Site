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
all_tags = basenames_df['unique_basename'].values

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
        Chain identifiers from CIF file.
    chain_labels : dict
        Dictionary of known chain labels from YAML config.
    
    Returns
    -------
    dict
        Mapping of chain ID to display label.
    
    Examples
    --------
    >>> mapping = build_chain_mapping(['n', 's'], {'n': 'PSII-n'})
    >>> mapping['n']
    'PSII-n'
    >>> mapping['s']
    'Psbs'
    """
    chain_id_to_label = {}
    for chain_id in chains:
        if chain_id in chain_labels:
            chain_id_to_label[chain_id] = chain_labels[chain_id]
        else:
            chain_id_to_label[chain_id] = f"Psb{chain_id}"
    
    # Add labels for PsbS segments
    chain_id_to_label['n_psbs_1'] = "PsbS (A)"
    chain_id_to_label['n_psbs_2'] = "PsbS (B)"
    
    return chain_id_to_label






# Function to plot a single sequence in a given axes
def plot_sequence(ax, one_letter_seq, lifetime_values, label, cmap, norm, helix_definitions=None):
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
    # Draw helices from definitions if available
    if helix_definitions:
        for helix_id, helix_info in helix_definitions.items():
            start_res = helix_info['start']
            end_res = helix_info['end']
            
            # Only draw if helix is within the sequence bounds
            if start_res < len(one_letter_seq) and end_res <= len(one_letter_seq):
                helix_left = (start_res - 1) * box_spacing - 0.8
                helix_right = (end_res - 1) * box_spacing + 0.8
                helix_width = helix_right - helix_left
                helix_center = (helix_left + helix_right) / 2
                
                ax.add_patch(plt.Rectangle((helix_left, helix_boxes_y), helix_width, 0.2,
                                          facecolor=helix_facecolor, edgecolor=helix_edgecolor, alpha=helix_alpha))
                ax.text(helix_center, helix_labels_y, helix_id, ha='center', va='center', 
                        fontsize=helix_label_fontsize, color=helix_label_color, fontweight=helix_label_fontweight)

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


# ============================================================================
# Plotting Support Functions
# ============================================================================

def get_helix_key_for_chain(chain_id):
    """Map chain ID to its helix definition key.
    
    Parameters
    ----------
    chain_id : str
        Chain identifier from CIF file.
    
    Returns
    -------
    str
        Key to use for looking up helix definitions.
        PsbS segments are mapped to chain_A and chain_B.
    
    Examples
    --------
    >>> get_helix_key_for_chain('n')
    'chain_n'
    >>> get_helix_key_for_chain('n_psbs_1')
    'chain_A'
    """
    if chain_id == 'n_psbs_1':
        return 'chain_A'
    elif chain_id == 'n_psbs_2':
        return 'chain_B'
    else:
        return f"chain_{chain_id}"


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


def create_figure_and_axes(num_chains, max_seq_length, box_spacing, box_width):
    """Create matplotlib figure and axes for multi-chain plotting.
    
    Parameters
    ----------
    num_chains : int
        Number of chains to plot.
    max_seq_length : int
        Maximum length of sequences.
    box_spacing : float
        Spacing between sequence boxes.
    box_width : float
        Width of each sequence box.
    
    Returns
    -------
    tuple
        (fig, axes) - matplotlib Figure and Axes objects.
    """
    total_width = (max_seq_length - 1) * box_spacing + box_width
    figure_width = max(20, total_width * 0.1)
    figure_height = 2.5 * num_chains
    
    fig, axes = plt.subplots(num_chains, 1, figsize=(figure_width, figure_height), sharex=False)
    if num_chains == 1:
        axes = [axes]
    
    return fig, axes


def plot_all_chains(axes, chains_data, chain_id_to_label, cofactors_data, 
                    helix_config, cmap, norm, max_seq_length):
    """Plot all chains on their respective axes.
    
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
    """
    for idx, (chain_id, chain_data) in enumerate(chains_data.items()):
        one_letter_seq = chain_data['sequence']
        lifetime_values = chain_data['bfactors']
        label = chain_id_to_label[chain_id]
        
        helix_key = get_helix_key_for_chain(chain_id)
        chain_helices = helix_config.get(helix_key, {})
        
        plot_sequence(axes[idx], one_letter_seq, lifetime_values, label, cmap, norm, 
                     helix_definitions=chain_helices)
        
        if chain_id in cofactors_data:
            chain_cofactors = cofactors_data[chain_id]
            plot_cofactors(axes[idx], chain_cofactors['labels'], chain_cofactors['bfactors'], 
                          cmap, norm, max_seq_length)
        else:
            axes[idx].set_ylim(-0.3, plot_title_y + 0.3)


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
    
    # Load tags from basenames CSV
    basenames_df = pd.read_csv(BASENAMES_CSV, header=0, sep=' ')
    all_tags = basenames_df['unique_basename'].values
    
    # Create output directory if it doesn't exist
    Path(OUTPUT_FIGURES_DIR).mkdir(parents=True, exist_ok=True)
    
    # Setup colormap and normalization (shared across all plots)
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
            
            # Extract cofactors
            cofactors_data = extract_cofactors_by_chain(cif_cofactors)
            
            # Build chain labels
            chain_id_to_label = build_chain_mapping(chains_list, chain_labels)
            
            # Validate chains
            num_chains = len(chains_data)
            if num_chains == 0:
                print(f"  Warning: No chains found in {cif_protein}")
                continue
            
            # Pad sequences and get max length
            max_seq_length = pad_sequences_to_max_length(chains_data)
            
            # Create figure and axes
            fig, axes = create_figure_and_axes(num_chains, max_seq_length, box_spacing, box_width)
            
            # Plot all chains
            plot_all_chains(axes, chains_data, chain_id_to_label, cofactors_data,
                          helix_config, cmap, norm, max_seq_length)
            
            # Save figure
            plt.tight_layout()
            plt.savefig(output_figure, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            print(f"  ✓ Saved sequence plots for {num_chains} chains to {output_figure}")
            
        except Exception as e:
            print(f"  Error processing {current_case}: {e}")
            plt.close('all')
            continue


if __name__ == '__main__':
    main()
    print("\nAll cases processed successfully!")

