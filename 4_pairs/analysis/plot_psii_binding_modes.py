"""
Render PSII structure with PsbS binding modes colored by occupancy.

Arguments:
-ref: Path to reference PSII structure PDB file.
-binding_dir: Directory containing binding site PDB structures.
-occupancy_csv: CSV file with occupancy data for binding sites.
-output_dir: Directory to save output figures.


"""

import os
import numpy as np
os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import warnings
warnings.filterwarnings("ignore")
import argparse

def parser():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Plot PSII structure with PsbS binding modes colored by occupancy.')
    parser.add_argument('-ref', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb',
                        help='Path to reference PSII structure PDB file.')
    parser.add_argument('-binding_dir', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned',
                        help='Directory containing binding site PDB structures.')
    parser.add_argument('-occupancy_csv', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped/occupancy.csv',
                        help='CSV file with occupancy data for binding sites.')
    parser.add_argument('-output_dir', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/figures',
                        help='Directory to save output figures.')
    return parser.parse_args()


def load_reference_structure(pdb_path):
    """Load reference PSII structure.
    
    Parameters
    ----------
    pdb_path : str
        Path to reference PDB file
        
    Returns
    -------
    mda.Universe
        MDAnalysis Universe object containing the reference structure
    """
    return mda.Universe(pdb_path)


def load_binding_site_structures(binding_dir):
    """Load all binding site PDB structures from directory.
    
    Parameters
    ----------
    binding_dir : str
        Directory containing PDB files
        
    Returns
    -------
    dict
        Dictionary mapping PDB file paths to Universe objects
    """
    pdb_files = [f for f in os.listdir(binding_dir) if f.endswith('.pdb')]
    pdb_files = [os.path.join(binding_dir, f) for f in pdb_files]
    
    universes = [mda.Universe(f) for f in pdb_files]
    
    return {pdb_files[i]: universes[i] for i in range(len(pdb_files))}


def load_occupancy_data(csv_path):
    """Load occupancy data from CSV file.
    
    Parameters
    ----------
    csv_path : str
        Path to occupancy CSV file
        
    Returns
    -------
    pd.DataFrame
        Occupancy data with columns 'trajectory' and 'occupancy_percent'
    """
    return pd.read_csv(csv_path)


def annotate_chain_dict_with_occupancy(chain_dict, df):
    """Add occupancy values to chain_dict from occupancy DataFrame.
    
    Parameters
    ----------
    chain_dict : dict
        Dictionary mapping PDB paths to Universe objects
    df : pd.DataFrame
        Occupancy data with columns 'trajectory' and 'occupancy_percent'
        
    Returns
    -------
    dict
        Modified chain_dict with occupancy attribute added to each Universe
    """
    for key in chain_dict:
        filename_no_ext = os.path.splitext(os.path.basename(key))[0]
        occupancy_row = df[df['trajectory'] == filename_no_ext]
        
        if not occupancy_row.empty:
            chain_dict[key].occupancy = occupancy_row['occupancy_percent'].values[0]
        else:
            chain_dict[key].occupancy = 0
    
    return chain_dict


def sort_chain_dict_by_occupancy(chain_dict):
    """Sort chain dictionary by occupancy in ascending order.
    
    Parameters
    ----------
    chain_dict : dict
        Dictionary with occupancy attribute on each Universe
        
    Returns
    -------
    dict
        Sorted dictionary (lowest to highest occupancy)
    """
    return dict(sorted(chain_dict.items(), 
                      key=lambda item: item[1].occupancy))


def setup_figure_and_axis(figsize=(12, 10)):
    """Create matplotlib figure and axis.
    
    Parameters
    ----------
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (12, 10)
        
    Returns
    -------
    tuple
        (fig, ax) - matplotlib figure and axis objects
    """
    return plt.subplots(figsize=figsize)


def plot_psii_background(ax, u0):
    """Plot PSII background atoms on axis.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    """
    selection = u0.select_atoms("all")
    ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
              c='lightgray', alpha=1, s=60)


def plot_reference_chains(ax, u0, cmap):
    """Plot reference protein chains (CP24, CP29, CP26, CP43).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    cmap : matplotlib.colors.Colormap
        Colormap object for coloring chains
    """
    colors_plot = [cmap(4), cmap(5), cmap(1), cmap(7), 
                   cmap(4), cmap(5), cmap(1), cmap(7)]
    labels_plot = ["CP24", "CP29", "CP26", "CP43", 
                   "CP24", "CP26", "CP29", "CP43"]
    chain_sel = ["8", "r", "s", "c", "4", "R", "S", "C"]
    
    for chain, color, label in zip(chain_sel, colors_plot, labels_plot):
        selection = u0.select_atoms(f"chainID {chain}")
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                  color=color, alpha=1, s=60)
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='grey', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0))


def plot_lhcii_complexes(ax, u0):
    """Plot LHCII antenna complexes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    """
    lhcii_chains = [["5", "6", "7"], ["1", "2", "3"], 
                    ["g", "y", "n"], ["G", "Y", "N"]]
    lhcii_labels = ["S-LHCII", "S-LHCII", "M-LHCII", "M-LHCII"]
    lhcii_colors = ["#DBE4C9", "#DBE4C9", "#87cd9d", "#87cd9d"]
    
    for chains, color, label in zip(lhcii_chains, lhcii_colors, lhcii_labels):
        selection = u0.select_atoms(f"chainID {' '.join(chains)}")
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                  c=color, alpha=1, s=60)
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='grey', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0))


def calculate_normalized_occupancy(value, max_occupancy):
    """Calculate normalized occupancy value for colormap.
    
    Parameters
    ----------
    value : float
        Occupancy value
    max_occupancy : float
        Maximum occupancy for normalization
        
    Returns
    -------
    float
        Normalized value between 0 and 1
    """
    return value / max_occupancy if max_occupancy > 0 else 0.5


def plot_binding_sites(ax, sorted_chain_dict, df, cmap_rdpu, max_occupancy):
    """Plot PsbS binding sites with occupancy-based coloring.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    sorted_chain_dict : dict
        Sorted dictionary of binding site structures
    df : pd.DataFrame
        Occupancy data
    cmap_rdpu : matplotlib.colors.Colormap
        PuRd colormap for occupancy visualization
    max_occupancy : float
        Maximum occupancy value for normalization
        
    Returns
    -------
    list
        List of site labels applied to each binding site
    """
    n_sites = len(sorted_chain_dict)
    site_labels = [f"S{n_sites - i}" for i in range(n_sites)]
    
    for (pdb_path, u), label in zip(sorted_chain_dict.items(), site_labels):
        basename = os.path.splitext(os.path.basename(pdb_path))[0]
        
        # Find matching occupancy value
        matching_rows = df[df['trajectory'] == basename]
        
        if matching_rows.empty:
            print(f"Warning: No matching trajectory found for {basename}")
            norm_value = 0.5
        else:
            value = matching_rows.iloc[0]['occupancy_percent']
            norm_value = calculate_normalized_occupancy(value, max_occupancy)
        
        # Plot binding site with occupancy color
        color = cmap_rdpu(norm_value)
        selection_str = ("chainID 9 and (not resname *MG* *HEME* *GG* *SQ* *PG* "
                        "W* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR CHL CLA CLB)")
        sel = u.select_atoms(selection_str)
        
        print(f"Plotting {label}: occupancy={norm_value:.2f}, atoms={sel.n_atoms}")
        
        ax.scatter(sel.positions[:, 0], sel.positions[:, 1], 
                  c='black', alpha=1, s=90)
        ax.scatter(sel.positions[:, 0], sel.positions[:, 1], 
                  c=color, alpha=1, s=60)
        
        # Add label at center of geometry
        cog = sel.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='black', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    return site_labels


def add_colorbar(ax, cmap_rdpu, vmin=0, vmax=35):
    """Add occupancy colorbar to plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to attach colorbar to
    cmap_rdpu : matplotlib.colors.Colormap
        PuRd colormap
    vmin : float, optional
        Minimum value for colorbar. Default is 0
    vmax : float, optional
        Maximum value for colorbar. Default is 35
    """
    sm = plt.cm.ScalarMappable(cmap=cmap_rdpu, 
                              norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Probability (%)', fontsize=20)
    cbar.ax.tick_params(labelsize=20)


def format_plot(ax):
    """Format plot appearance.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to format
    """
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')


def save_figure(fig, output_dir, filename_base='psii_binding_sites_overview'):
    """Save figure to disk in multiple formats.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save
    output_dir : str
        Output directory path
    filename_base : str, optional
        Base filename without extension. Default is 'psii_binding_sites_overview'
    """
    os.makedirs(output_dir, exist_ok=True)
    
    png_path = f"{output_dir}/{filename_base}.png"
    pdf_path = f"{output_dir}/{filename_base}.pdf"
    
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    
    print(f"Figure saved in {png_path}")


def main():
    """Main execution function."""
    # Setup paths
    args = parser()
    paths = {
        'pdb0': args.ref,
        'binding_dir': args.binding_dir,
        'occupancy_csv': args.occupancy_csv,
        'output_dir': args.output_dir
    }
        
    # Load data
    u0 = load_reference_structure(paths['pdb0'])
    chain_dict = load_binding_site_structures(paths['binding_dir'])
    df = load_occupancy_data(paths['occupancy_csv'])
    
    # Process occupancy data
    chain_dict = annotate_chain_dict_with_occupancy(chain_dict, df)
    sorted_chain_dict = sort_chain_dict_by_occupancy(chain_dict)
    
    # Setup figure
    cmap = plt.get_cmap("Set3")
    cmap_rdpu = plt.get_cmap("PuRd")
    fig, ax = setup_figure_and_axis()
    
    # Plot components
    plot_psii_background(ax, u0)
    plot_reference_chains(ax, u0, cmap)
    plot_lhcii_complexes(ax, u0)
    
    # Plot binding sites
    max_occupancy = df['occupancy_percent'].max()
    plot_binding_sites(ax, sorted_chain_dict, df, cmap_rdpu, max_occupancy)
    
    # Finalize plot
    add_colorbar(ax, cmap_rdpu)
    format_plot(ax)
    
    # Save figure
    save_figure(fig, paths['output_dir'])


if __name__ == "__main__":
    main()