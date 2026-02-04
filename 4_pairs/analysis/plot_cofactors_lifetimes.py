"""
Render PSII structure with residues plotted as circles and colored by lifetime.

Arguments:
-ref: Path to reference PSII structure PDB file.
-lifetimes: CSV file with  the headers: resid, resname, chain, lifetime
-vmin: Minimum lifetime value for color mapping
-vmax: Maximum lifetime value for color mapping
-cmap: Colormap name for lifetime coloring
-log_transform: Whether to apply log transformation to lifetime values (True/False)
-output: Directory to save output figure

Workflow:
1. Data Loading:
   - Load reference PSII structure (reference PDB)
   - Load CSV file with  the headers: resid, resname, chain, lifetime 

2. Data Processing:
   - Add the x, y positions of each residue to the dataframe from the reference structure
   - Calculate normalized occupancy values for color mapping

3. Figure Setup:
   - Create matplotlib figure and axis

4. Plotting Layers:
   - Plot PSII background atoms (light gray)
   - Plot reference protein chains (CP24, CP29, CP26, CP43) with distinct colors
   - Plot LHCII antenna complexes (S-LHCII and M-LHCII)
   - Plot residues colored by lifetime using specified colormap as circles with black borders and a size of 60

5. Finalization:
   - Format plot (equal aspect, remove axes)
   - Save figure as eps with a 600 dpi resolution

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
import cmap

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot PSII residues colored by lifetimes.")
    parser.add_argument('-ref', type=str, required=True, help='Path to reference PSII structure PDB file.')
    parser.add_argument('-lifetimes', type=str, required=True, help='CSV file with headers: resid, resname, chain, lifetime')
    parser.add_argument('-vmin', type=float, default=None, help='Minimum lifetime value for color mapping')
    parser.add_argument('-vmax', type=float, default=None, help='Maximum lifetime value for color mapping')
    parser.add_argument('-cmap', type=str, default='viridis', help='Colormap name for lifetime coloring')
    parser.add_argument('-log_transform', action='store_true', help='Apply log transformation to lifetime values')
    parser.add_argument('-output', type=str, required=True, help='Directory to save output figure')
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

def load_csv(csv_path):
    """Load CSV file with lifetimes.
    
    Parameters
    ----------
    csv_path : str
        Path to CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the lifetimes data
    """
    return pd.read_csv(csv_path)

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
                   "CP24", "CP29", "CP26", "CP43"]
    chain_sel = ["8", "r", "s", "c", "4", "R", "S", "C"]
    
    for chain, color, label in zip(chain_sel, colors_plot, labels_plot):
        selection = u0.select_atoms(f"chainID {chain}")
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                  color=color, alpha=1, s=60)
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='black', fontsize=12, 
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
        ax.text(cog[0], cog[1], label, color='black', fontsize=12, 
               ha='center', va='center',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0))


def add_positions_to_dataframe(df, u0):
    """Add x, y positions from reference structure to dataframe.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns: resid, resname, chain, lifetime
    u0 : mda.Universe
        Reference PSII Universe
        
    Returns
    -------
    pd.DataFrame
        DataFrame with added x, y position columns
    """
    positions = []
    
    for idx, row in df.iterrows():
        resid = row['resid']
        resname = row['resname']
        chain = row['chain']
        
        # Select the specific residue
        sel = u0.select_atoms(f"resid {resid} and resname {resname} and chainID {chain}")
        
        if sel.n_atoms > 0:
            cog = sel.center_of_geometry()
            positions.append({'x': cog[0], 'y': cog[1]})
        else:
            print(f"Warning: No atoms found for resid {resid}, resname {resname}, chain {chain}")
            positions.append({'x': np.nan, 'y': np.nan})
    
    # Add positions to dataframe
    pos_df = pd.DataFrame(positions)
    df['x'] = pos_df['x']
    df['y'] = pos_df['y']
    
    return df


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


def plot_residues_by_lifetime(ax, df, cmap_lifetime, vmin, vmax, log_transform=False):
    """Plot residues as circles colored by lifetime.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    df : pd.DataFrame
        DataFrame with columns: resid, resname, chain, lifetime, x, y
    cmap_lifetime : matplotlib.colors.Colormap
        Colormap for lifetime visualization
    vmin : float
        Minimum lifetime value for color mapping
    vmax : float
        Maximum lifetime value for color mapping
    log_transform : bool, optional
        Whether to apply log transformation to lifetime values
    """
    # Remove rows with NaN positions
    df_valid = df.dropna(subset=['x', 'y'])
    
    # Get lifetime values
    lifetimes = df_valid['lifetime'].values
    
    # Apply log transformation if requested
    if log_transform:
        lifetimes = np.log10(lifetimes + 1)  # Add 1 to avoid log(0)
        if vmin is not None:
            vmin = np.log10(vmin + 1)
        if vmax is not None:
            vmax = np.log10(vmax + 1)
    
    # Determine vmin and vmax if not provided
    if vmin is None:
        vmin = lifetimes.min()
    if vmax is None:
        vmax = lifetimes.max()
    
    # Normalize lifetime values
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    colors = cmap_lifetime(norm(lifetimes))
    
    # Plot residues with black border
    ax.scatter(df_valid['x'].values, df_valid['y'].values, 
              c='black', alpha=1, s=90, zorder=10)
    ax.scatter(df_valid['x'].values, df_valid['y'].values, 
              c=colors, alpha=1, s=60, zorder=11)
    
    print(f"Plotted {len(df_valid)} residues with lifetimes")


def add_colorbar(ax, cmap_lifetime, vmin, vmax, label='Lifetime', log_transform=False):
    """Add lifetime colorbar to plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to attach colorbar to
    cmap_lifetime : matplotlib.colors.Colormap
        Colormap for lifetime
    vmin : float
        Minimum value for colorbar
    vmax : float
        Maximum value for colorbar
    label : str, optional
        Label for colorbar. Default is 'Lifetime'
    log_transform : bool, optional
        Whether log transformation was applied
    """
    # Apply log transformation to vmin/vmax if needed
    if log_transform:
        display_vmin = np.log10(vmin + 1) if vmin is not None else 0
        display_vmax = np.log10(vmax + 1) if vmax is not None else 1
    else:
        display_vmin = vmin if vmin is not None else 0
        display_vmax = vmax if vmax is not None else 1
    
    sm = plt.cm.ScalarMappable(cmap=cmap_lifetime, 
                              norm=plt.Normalize(vmin=display_vmin, vmax=display_vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(label, fontsize=20)
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


def save_figure(fig, output_dir, filename_base='psii_residues_lifetimes'):
    """Save figure to disk as EPS.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save
    output_dir : str
        Output directory path
    filename_base : str, optional
        Base filename without extension. Default is 'psii_residues_lifetimes'
    """
    os.makedirs(output_dir, exist_ok=True)
    
    eps_path = f"{output_dir}/{filename_base}.eps"
    
    fig.savefig(eps_path, dpi=600, bbox_inches='tight', format='eps')
    
    print(f"Figure saved in {eps_path}")


def main():
    """Main execution function."""
    # Parse arguments
    args = parse_arguments()
        
    # Load data
    u0 = load_reference_structure(args.ref)
    df = load_csv(args.lifetimes)

    # Add positions from structure
    df = add_positions_to_dataframe(df, u0)
    
    # Setup figure
    cmap_ref_chains = plt.get_cmap("Set3")
    cmap_sites = plt.get_cmap(args.cmap)
    fig, ax = setup_figure_and_axis()
    
    # Plot components
    plot_psii_background(ax, u0)
    plot_reference_chains(ax, u0, cmap_ref_chains)
    plot_lhcii_complexes(ax, u0)
    
    # Plot residues colored by lifetime
    plot_residues_by_lifetime(ax, df, cmap_sites, args.vmin, args.vmax, args.log_transform)
    
    # Finalize plot
    add_colorbar(ax, cmap_sites, vmin=args.vmin, vmax=args.vmax, 
                label='Lifetime', log_transform=args.log_transform)
    format_plot(ax)
    
    # Save figure
    save_figure(fig, args.output)


if __name__ == "__main__":
    main()


