"""
Render PSII structure with PsbS binding modes colored by lifetime.

Arguments:
-ref: Path to reference PSII structure PDB file.
-binding_dir: Directory containing binding site PDB structures.
-binding_modes_occupancy_csv: CSV file with occupancy data for binding sites (used for sorting/labeling).
-chains_occupancy_csv: CSV file with per-chain lifetime sums (columns: chain, sum_ns).
-output_dir: Directory to save output figures.
-vmin: Minimum value for color mapping.
-vmax: Maximum value for color mapping.
-log_transform: Apply log10 transformation to lifetime values and vmin/vmax.
-cmap: Colormap name for lifetime coloring (uses cmap library).
-move_site_label: Label to move site text away from center (default: None), e.g., "S10".
-move_offset: Offset for moving site label (default: (0,0)), e.g., (5, -5).

Workflow:
1. Data Loading:
   - Load reference PSII structure (reference PDB)
   - Load all PsbS binding site structures (from binding_dir)
   - Load occupancy/lifetime data from CSV files

2. Data Processing:
   - Annotate each binding site with its lifetime
   - Sort binding sites by lifetime (lowest to highest)
   - Apply log10 transformation if requested

3. Figure Setup:
   - Create matplotlib figure and axis
   - Initialize colormaps

4. Plotting Layers:
   - Plot PSII background atoms (light gray)
   - Plot reference protein chains (CP24, CP29, CP26, CP43) with distinct colors
   - Plot LHCII antenna complexes colored by per-chain lifetime
   - Plot PsbS binding sites colored by lifetime
   - Add site labels (S1, S2, S3, etc.) at center of geometry

5. Finalization:
   - Add colorbar showing lifetime scale
   - Format plot (equal aspect, remove axes)
   - Save figure as PNG and PDF in output directory

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
from cmap import Colormap

def parser():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Plot PSII structure with PsbS binding modes colored by lifetime.')
    parser.add_argument('-ref', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb',
                        help='Path to reference PSII structure PDB file.')
    parser.add_argument('-binding_dir', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster',
                        help='Directory containing binding site PDB structures.')
    parser.add_argument('-binding_modes_occupancy_csv', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped/binding_modes_lifetimes.csv',
                        help='CSV file with lifetime data for binding modes (columns: trajectory, frames, lifetimes).')
    parser.add_argument('-chains_occupancy_csv', type=str, required=False,
                        default=None,
                        help='CSV file with per-chain lifetime sums (columns: chain, sum_ns). Used to color chains by lifetime.')
    parser.add_argument('-output_dir', type=str, required=False,
                        default='/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures',
                        help='Directory to save output figures.')
    parser.add_argument('-vmin', type=float, default=None,
                        help='Minimum value for color mapping.')
    parser.add_argument('-vmax', type=float, default=None,
                        help='Maximum value for color mapping.')
    parser.add_argument('-log_transform', action='store_true',
                        help='Apply log10 transformation to lifetime values.')
    parser.add_argument('-cmap', type=str, default='colorcet:CET_L17',
                        help='Colormap name for lifetime coloring (uses cmap library).')
    parser.add_argument('-move_site_label', type=str, required=False,
                        default=None,
                        help='Label to move site text away from center (e.g., "S10").')
    parser.add_argument('-move_offset', type=str, required=False,
                        default="0 0",
                        help='Offset for moving site label (e.g., (5, -5)).')

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
        Example: {pdb_path: mda.Universe}
    """
    pdb_files = [f for f in os.listdir(binding_dir) if f.endswith('.pdb')]
    pdb_files = [os.path.join(binding_dir, f) for f in pdb_files]
    
    universes = [mda.Universe(f) for f in pdb_files]
    
    return {pdb_files[i]: universes[i] for i in range(len(pdb_files))}


def load_lifetime_data(csv_path):
    """Load lifetime data from CSV file.
    
    Parameters
    ----------
    csv_path : str
        Path to lifetime CSV file (columns: trajectory, frames, lifetimes)
        
    Returns
    -------
    pd.DataFrame
        Lifetime data
    """
    return pd.read_csv(csv_path)


def load_chains_lifetime_data(csv_path):
    """Load per-chain lifetime data from CSV file.
    
    Parameters
    ----------
    csv_path : str
        Path to chains lifetime CSV file (columns: chain, sum_ns)
        
    Returns
    -------
    dict
        Dictionary mapping chain ID to sum_ns value
    """
    if csv_path is None:
        return {}
    df = pd.read_csv(csv_path)
    return dict(zip(df['chain'].astype(str), df['sum_ns']))


def annotate_chain_dict_with_lifetime(chain_dict, df):
    """Add lifetime values to chain_dict from lifetime DataFrame.
    
    Parameters
    ----------
    chain_dict : dict
        Dictionary mapping PDB paths to Universe objects
    df : pd.DataFrame
        Lifetime data with columns 'trajectory' and 'lifetimes'
        
    Returns
    -------
    dict
        Modified chain_dict with lifetime attribute added to each Universe
    """
    for key in chain_dict:
        filename_no_ext = os.path.splitext(os.path.basename(key))[0]
        lifetime_row = df[df['trajectory'] == filename_no_ext]
        
        if not lifetime_row.empty:
            chain_dict[key].lifetime = lifetime_row['lifetimes'].values[0]
        else:
            chain_dict[key].lifetime = 0
    
    return chain_dict


def sort_chain_dict_by_lifetime(chain_dict):
    """Sort chain dictionary by lifetime in ascending order.
    
    Parameters
    ----------
    chain_dict : dict
        Dictionary with lifetime attribute on each Universe
        
    Returns
    -------
    dict
        Sorted dictionary (lowest to highest lifetime)
    """
    return dict(sorted(chain_dict.items(), 
                      key=lambda item: item[1].lifetime))


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


def plot_reference_chains(ax, u0, cmap, chains_lifetime_dict=None, lifetime_cmap=None, norm=None):
    """Plot reference protein chains (CP24, CP29, CP26, CP43).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    cmap : matplotlib.colors.Colormap
        Default colormap object for coloring chains (Set3)
    chains_lifetime_dict : dict or None
        Dictionary mapping chain ID to lifetime value. If provided, chains with
        data are colored by lifetime instead of default colors.
    lifetime_cmap : colormap or None
        Colormap for lifetime coloring
    norm : matplotlib.colors.Normalize or None
        Normalization for lifetime colormap
    """
    default_colors = [cmap(4), cmap(5), cmap(1), cmap(7), 
                      cmap(4), cmap(5), cmap(1), cmap(7)]
    labels_plot = ["CP24", "CP29", "CP26", "CP43", 
                   "CP24", "CP29", "CP26", "CP43"]
    chain_sel = ["8", "r", "s", "c", "4", "R", "S", "C"]
    
    for chain, default_color, label in zip(chain_sel, default_colors, labels_plot):
        selection = u0.select_atoms(f"chainID {chain}")
        
        # Use lifetime color if chain has lifetime data, otherwise default
        if chains_lifetime_dict and chain in chains_lifetime_dict and lifetime_cmap is not None and norm is not None:
            color = lifetime_cmap(norm(chains_lifetime_dict[chain]))
        else:
            color = default_color
        
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                  color=color, alpha=1, s=60)
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='black', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0))


def plot_lhcii_complexes(ax, u0, chains_lifetime_dict=None, lifetime_cmap=None, norm=None):
    """Plot LHCII antenna complexes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    chains_lifetime_dict : dict or None
        Dictionary mapping chain ID to lifetime value
    lifetime_cmap : colormap or None
        Colormap for lifetime coloring
    norm : matplotlib.colors.Normalize or None
        Normalization for lifetime colormap
    """
    lhcii_chains = [["5", "6", "7"], ["1", "2", "3"], 
                    ["g", "y", "n"], ["G", "Y", "N"]]
    lhcii_labels = ["S-LHCII", "S-LHCII", "M-LHCII", "M-LHCII"]
    lhcii_default_colors = ["#DBE4C9", "#DBE4C9", "#87cd9d", "#87cd9d"]
    
    for chains, default_color, label in zip(lhcii_chains, lhcii_default_colors, lhcii_labels):
        # Plot each chain individually so we can color by lifetime
        for chain in chains:
            selection = u0.select_atoms(f"chainID {chain}")
            if chains_lifetime_dict and chain in chains_lifetime_dict and lifetime_cmap is not None and norm is not None:
                color = lifetime_cmap(norm(chains_lifetime_dict[chain]))
            else:
                color = default_color
            ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                      c=[color], alpha=1, s=60)
        
        # Group label at center of geometry of all chains
        selection = u0.select_atoms(f"chainID {' '.join(chains)}")
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='black', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0))


def plot_binding_sites(ax, sorted_chain_dict, df, lifetime_cmap, norm, move_site_label=None, move_offset=(0,0)):
    """Plot PsbS binding sites with lifetime-based coloring.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    sorted_chain_dict : dict
        Sorted dictionary of binding site structures
    df : pd.DataFrame
        Lifetime data with columns 'trajectory' and 'lifetimes'
    lifetime_cmap : colormap
        Colormap for lifetime visualization
    norm : matplotlib.colors.Normalize
        Normalization object for color mapping
    move_site_label : str or None
        Label to move
    move_offset : str
        Offset for moving site label
        
    Returns
    -------
    list
        List of site labels applied to each binding site
    """
    n_sites = len(sorted_chain_dict)
    site_labels = [f"S{n_sites - i}" for i in range(n_sites)]
    
    for (pdb_path, u), label in zip(sorted_chain_dict.items(), site_labels):
        basename = os.path.splitext(os.path.basename(pdb_path))[0]
        
        # Find matching lifetime value
        matching_rows = df[df['trajectory'] == basename]
        
        if matching_rows.empty:
            print(f"Warning: No matching trajectory found for {basename}")
            lifetime_value = 0
        else:
            lifetime_value = matching_rows.iloc[0]['lifetimes']
        
        # Plot binding site with lifetime color
        color = lifetime_cmap(norm(lifetime_value))
        selection_str = ("chainID 9 and (not resname *MG* *HEME* *GG* *SQ* *PG* "
                        "W* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR CHL CLA CLB)")
        sel = u.select_atoms(selection_str)
        
        print(f"Plotting {label}: lifetime={lifetime_value:.1f}, atoms={sel.n_atoms}")
        
        ax.scatter(sel.positions[:, 0], sel.positions[:, 1], 
                  c='black', alpha=0, s=90)
        ax.scatter(sel.positions[:, 0], sel.positions[:, 1], 
                  c=[color], alpha=.2, s=60)
        
        # Add label at center of geometry
        cog = sel.center_of_geometry()
        cog_x = cog[0]
        cog_y = cog[1]

        if label == move_site_label:
            move_offset = np.array(move_offset.split(), dtype=float)
            cog_x += move_offset[0]
            cog_y += move_offset[1]

        ax.text(cog_x, cog_y, label, color='black', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=1))
    
    return site_labels


def add_colorbar(ax, lifetime_cmap, display_vmin, display_vmax, original_vmin=None, original_vmax=None, log_transform=False):
    """Add lifetime colorbar to plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to attach colorbar to
    lifetime_cmap : colormap
        Colormap for lifetime visualization
    display_vmin : float
        Minimum value for colorbar (in display scale)
    display_vmax : float
        Maximum value for colorbar (in display scale)
    original_vmin : float or None
        Original vmin before log transform
    original_vmax : float or None
        Original vmax before log transform
    log_transform : bool
        Whether log transformation was applied
    """
    # Convert cmap library colormap to matplotlib colormap
    mpl_cmap = lifetime_cmap.to_mpl()
    sm = plt.cm.ScalarMappable(cmap=mpl_cmap, 
                              norm=plt.Normalize(vmin=display_vmin, vmax=display_vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    
    if log_transform:
        # Create tick positions at log10(0+1)=0, log10(10+1), log10(100+1), ... up to vmax
        tick_positions = [0]  # log10(0 + 1) = 0
        max_power = int(np.ceil(display_vmax))
        for p in range(1, max_power + 1):
            if p <= display_vmax:
                tick_positions.append(p)  # log10(10^p) = p
        
        # Add exact vmax position if not already there
        if original_vmax is not None and display_vmax not in tick_positions:
            tick_positions.append(display_vmax)
        
        tick_positions = sorted(tick_positions)
        
        # Create labels showing original values
        tick_labels = []
        for t in tick_positions:
            if t == 0:
                tick_labels.append('0')
            elif t == display_vmax and original_vmax is not None:
                tick_labels.append(f'{int(original_vmax):,}')
            else:
                tick_labels.append(f'{int(10**t):,}')
        
        cbar.set_ticks(tick_positions)
        cbar.set_ticklabels(tick_labels)
    
    cbar.set_label('Lifetime (ns)', fontsize=20)
    cbar.ax.tick_params(labelsize=14)


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
    # Parse arguments
    args = parser()
        
    # Load data
    u0 = load_reference_structure(args.ref)
    chain_dict = load_binding_site_structures(args.binding_dir)
    df = load_lifetime_data(args.binding_modes_occupancy_csv)
    chains_lifetime_dict = load_chains_lifetime_data(args.chains_occupancy_csv)
    
    # Process lifetime data
    chain_dict = annotate_chain_dict_with_lifetime(chain_dict, df)
    sorted_chain_dict = sort_chain_dict_by_lifetime(chain_dict)
    
    # Get lifetime values for color scale
    lifetime_values = df['lifetimes'].values
    
    # Determine vmin/vmax
    original_vmin = args.vmin if args.vmin is not None else lifetime_values.min()
    original_vmax = args.vmax if args.vmax is not None else lifetime_values.max()
    
    # Apply log transform if requested
    if args.log_transform:
        # Transform lifetime values in df
        df['lifetimes'] = np.log10(df['lifetimes'] + 1)
        # Transform chain lifetime values
        chains_lifetime_dict = {k: np.log10(v + 1) for k, v in chains_lifetime_dict.items()}
        # Transform vmin/vmax
        display_vmin = np.log10(original_vmin + 1)
        display_vmax = np.log10(original_vmax + 1)
        print(f"Log transform applied: vmin {original_vmin} -> {display_vmin:.2f}, vmax {original_vmax} -> {display_vmax:.2f}")
    else:
        display_vmin = original_vmin
        display_vmax = original_vmax
    
    print(f"Color scale: display_vmin={display_vmin:.2f}, display_vmax={display_vmax:.2f}")
    
    # Setup colormaps
    cmap_ref = plt.get_cmap("Set3")
    lifetime_cmap = Colormap(args.cmap)
    norm = plt.Normalize(vmin=display_vmin, vmax=display_vmax)
    
    fig, ax = setup_figure_and_axis()
    
    # Plot components (binding sites first so proteins are plotted in front)
    plot_psii_background(ax, u0)
    plot_binding_sites(ax, sorted_chain_dict, df, lifetime_cmap, norm,
                       move_site_label=args.move_site_label,
                       move_offset=args.move_offset)
    plot_reference_chains(ax, u0, cmap_ref, 
                         chains_lifetime_dict=chains_lifetime_dict,
                         lifetime_cmap=lifetime_cmap, norm=norm)
    plot_lhcii_complexes(ax, u0,
                        chains_lifetime_dict=chains_lifetime_dict,
                        lifetime_cmap=lifetime_cmap, norm=norm)
    
    # Finalize plot
    add_colorbar(ax, lifetime_cmap, display_vmin, display_vmax,
                original_vmin=original_vmin, original_vmax=original_vmax,
                log_transform=args.log_transform)
    format_plot(ax)
    
    # Save figure
    save_figure(fig, args.output_dir)


if __name__ == "__main__":
    main()