"""
Render PSII structure with residues plotted as circles and colored by lifetime.

Arguments:
-ref: Path to reference PSII structure PDB file.
-csv: CSV file with  the headers: resid, resname, chain, lifetime
-vmin: Minimum lifetime value for color mapping
-vmax: Maximum lifetime value for color mapping
-cmap: Colormap name for lifetime coloring
-log_transform: Whether to apply log transformation to lifetime values (True/False)
-output: Output file path (e.g., ./figures/plot.png)
-change_resnames_yaml: (Optional) YAML file with resid to resname mapping for renaming residues in the plot
-equivalent_chains_yaml: (Optional) YAML file with equivalent chains mapping.
-helix_labels_yaml: (Optional) YAML file with helix definitions for plotting helices B and C in green behind residues.
-change_chains_yaml: (Optional) YAML file with chain ID mapping for renaming chains.


Example usage:
python3 plot_cofactors_lifetimes.py -ref reference_psii.pdb -csv lifetimes.csv -vmin 0 -vmax 100 -cmap viridis -log_transform -output ./figures/plot.png


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
   - Save figure as PNG with a 600 dpi resolution

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
import yaml

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot PSII residues colored by lifetimes.")
    parser.add_argument('-ref', type=str, required=True, help='Path to reference PSII structure PDB file.')
    parser.add_argument('-csv', type=str, required=True, help='CSV file with headers: resid, resname, chain, lifetime')
    parser.add_argument('-vmin', type=float, default=None, help='Minimum lifetime value for color mapping')
    parser.add_argument('-vmax', type=float, default=None, help='Maximum lifetime value for color mapping')
    parser.add_argument('-cmap', type=str, default='viridis', help='Colormap name for lifetime coloring')
    parser.add_argument('-log_transform', action='store_true', help='Apply log transformation to lifetime values')
    parser.add_argument('-output', type=str, required=True, help='Output file path (e.g., ./figures/plot.png)')
    parser.add_argument('-change_resnames_yaml', type=str, default=None, help='YAML file with resid to resname mapping for renaming residues in the plot')
    parser.add_argument('-equivalent_chains_yaml', type=str, default=None, help='YAML file with equivalent chains mapping')
    parser.add_argument('-helix_labels_yaml', type=str, default=None, help='YAML file with helix definitions (for plotting helices B and C in green)')
    parser.add_argument('-change_chains_yaml', type=str, default=None, help='YAML file with chain ID mapping for renaming chains in the dataframe')
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
    """Load CSV file.
    
    Parameters
    ----------
    csv_path : str
        Path to CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the data
    """
    return pd.read_csv(csv_path)

def change_resnames(df, yaml_path):
    """Change residue names in dataframe based on YAML mapping.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'resname' column
    yaml_path : str
        Path to YAML file with old_resname: new_resname mapping
        
    Returns
    -------
    pd.DataFrame
        DataFrame with updated resname values
    """
    if yaml_path is None:
        return df
    
    # Load the YAML file
    with open(yaml_path, 'r') as f:
        resname_mapping = yaml.safe_load(f)
    
    # Create a mapping from old resname to new resname
    mapping = {old_resname: new_resname for old_resname, new_resname in resname_mapping.items()}
    
    # Update the resname column based on the mapping
    df['resname'] = df['resname'].map(lambda x: mapping.get(x, x))
    
    return df


def change_chains(df, yaml_path):
    """Change chain IDs in dataframe based on YAML mapping.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'chain' column
    yaml_path : str
        Path to YAML file with old_chain: new_chain mapping
        
    Returns
    -------
    pd.DataFrame
        DataFrame with updated chain values
    """
    if yaml_path is None:
        return df
    
    # Load the YAML file
    with open(yaml_path, 'r') as f:
        chain_mapping = yaml.safe_load(f)
    
    # Create a mapping from old chain to new chain
    mapping = {str(old_chain): str(new_chain) for old_chain, new_chain in chain_mapping.items()}
    
    # Update the chain column based on the mapping
    df['chain'] = df['chain'].map(lambda x: mapping.get(str(x), x))
    
    print(f"Changed chains: {mapping}")
    
    return df

def duplicate_entries_for_equivalent_chains(df, yaml_path):
    """Duplicate dataframe entries for equivalent chains.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'chain' column
    yaml_path : str
        Path to YAML file with equivalent chains mapping
        
    Returns
    -------
    pd.DataFrame
        DataFrame with duplicated entries for equivalent chains
    """
    if yaml_path is None:
        return df
    
    # Remove x,y columns if they exist (will be recalculated after duplication)
    columns_to_remove = ['x', 'y']
    df = df.drop(columns=[col for col in columns_to_remove if col in df.columns], errors='ignore')
    
    # Load the YAML file
    with open(yaml_path, 'r') as f:
        equivalent_chains = yaml.safe_load(f)
    
    # Create a list to store duplicated rows
    duplicated_rows = []
    
    # Iterate through each row in the dataframe
    for idx, row in df.iterrows():
        chain = row['chain']
        
        # Check if this chain has equivalents
        for key, equivalent_list in equivalent_chains.items():
            if chain in equivalent_list:
                # Duplicate this row for each equivalent chain
                for equiv_chain in equivalent_list:
                    new_row = row.copy()
                    new_row['chain'] = equiv_chain
                    duplicated_rows.append(new_row)
                break
        else:
            # No equivalent chains found, keep original row
            duplicated_rows.append(row)
    
    # Create new dataframe from duplicated rows
    df_duplicated = pd.DataFrame(duplicated_rows)
    
    print(f"Duplicated entries: {len(df)} -> {len(df_duplicated)} rows")
    
    return df_duplicated

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


def plot_psii_background(ax, u0, selected_chains=None, mode=None):
    """Plot PSII background atoms on axis.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    selected_chains : set or None
        If provided, only these chains are considered "selected"
    mode : str or None
        'white': non-selected chains in white, selected in lightgray
        'only_chains': only plot selected chains
        None: plot all in lightgray (default)
    """
    if mode == 'only_chains' and selected_chains is not None:
        sel_str = ' '.join(selected_chains)
        selection = u0.select_atoms(f"chainID {sel_str}")
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1],
                  c='lightgray', alpha=1, s=60)
    elif mode == 'white' and selected_chains is not None:
        # Plot non-selected chains in white
        all_chain_ids = set(u0.atoms.chainIDs)
        non_selected = all_chain_ids - selected_chains
        if non_selected:
            sel_str = ' '.join(non_selected)
            non_sel_atoms = u0.select_atoms(f"chainID {sel_str}")
            ax.scatter(non_sel_atoms.positions[:, 0], non_sel_atoms.positions[:, 1],
                      c='black', alpha=1, s=80)
            ax.scatter(non_sel_atoms.positions[:, 0], non_sel_atoms.positions[:, 1],
                      c='white', alpha=1, s=60)
        # Plot selected chains in lightgray
        sel_str = ' '.join(selected_chains)
        sel_atoms = u0.select_atoms(f"chainID {sel_str}")
        ax.scatter(sel_atoms.positions[:, 0], sel_atoms.positions[:, 1],
                  c='lightgray', alpha=1, s=60)
    else:
        selection = u0.select_atoms("all")
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1],
                  c='lightgray', alpha=1, s=60)


def plot_reference_chains(ax, u0, cmap, selected_chains=None, mode=None):
    """Plot reference protein chains (CP24, CP29, CP26, CP43).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    cmap : matplotlib.colors.Colormap
        Colormap object for coloring chains
    selected_chains : set or None
        If provided, only these chains are considered "selected"
    mode : str or None
        'white': non-selected chains in white, 'only_chains': skip non-selected
    """
    colors_plot = [cmap(4), cmap(5), cmap(1), cmap(7), 
                   cmap(4), cmap(5), cmap(1), cmap(7)]
    labels_plot = ["CP24", "CP29", "CP26", "CP43", 
                   "CP24", "CP29", "CP26", "CP43"]
    chain_sel = ["8", "r", "s", "c", "4", "R", "S", "C"]
    
    for chain, color, label in zip(chain_sel, colors_plot, labels_plot):
        is_selected = selected_chains is None or chain in selected_chains
        
        if mode == 'only_chains' and not is_selected:
            continue
        
        selection = u0.select_atoms(f"chainID {chain}")
        display_color = color if is_selected or mode is None else 'white'
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                  color="black", alpha=1, s=80)
        ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                  color=display_color, alpha=1, s=60)
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='black', fontsize=8, 
               ha='center', va='center', zorder=20,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="None", alpha=0))


def plot_lhcii_complexes(ax, u0, selected_chains=None, mode=None):
    """Plot LHCII antenna complexes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    selected_chains : set or None
        If provided, only these chains are considered "selected"
    mode : str or None
        'white': non-selected chains in white, 'only_chains': skip non-selected
    """
    lhcbm_color = "#DBE4C9"
    lhcb3_color = "#87cd9d"
    lhcii_chains = [["7", "6", "5"], ["3", "2", "1"], 
                    ["n", "y", "g"], ["N", "Y", "G"]]
    lhcii_labels = ["S-LHCII", "S-LHCII", "M-LHCII", "M-LHCII"]
    lhcii_colors = [[lhcb3_color, lhcbm_color, lhcbm_color], [lhcb3_color, lhcbm_color, lhcbm_color], 
                    [lhcbm_color, lhcbm_color, lhcbm_color], [lhcbm_color, lhcbm_color, lhcbm_color]]
    
    # Chains that will have lifetime labels plotted on them - use higher zorder
    labeled_chains = ["n", "g", "6", "7", "8", "s"]
    
    for chains, colors, label in zip(lhcii_chains, lhcii_colors, lhcii_labels):
        for chain, color in zip(chains, colors):
            is_selected = selected_chains is None or chain in selected_chains
            
            if mode == 'only_chains' and not is_selected:
                continue
            
            selection = u0.select_atoms(f"chainID {chain}")
            display_color = color if is_selected or mode is None else 'white'
            # Plot LHCII with higher zorder for labeled chains
            zorder = 5 if chain in labeled_chains else 3
            ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                      c="black", alpha=1, s=80, zorder=zorder)
            ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                      c=display_color, alpha=1, s=60, zorder=zorder)
        
        # Plot label at center of geometry of all chains in this group
        selection = u0.select_atoms(f"chainID {' '.join(chains)}")
        cog = selection.center_of_geometry()
        ax.text(cog[0], cog[1], label, color='black', fontsize=8, fontweight='bold',
               ha='center', va='center', zorder=20,
               bbox=dict(boxstyle="round,pad=0.3", facecolor=None, alpha=0))
    
    # Add individual chain labels for labeled chains
    # LHCBM labels at center of chains n, g, 6
    for chain in ["n", "g", "6"]:
        selection = u0.select_atoms(f"chainID {chain}")
        if selection.n_atoms > 0:
            cog = selection.center_of_geometry()
            ax.text(cog[0], cog[1], "LHCBM", color='black', fontsize=8,
                   ha='center', va='center', zorder=20,
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="None", alpha=0))
    
    # LHCB3 label at center of chain 7
    for chain in ["7"]:
        selection = u0.select_atoms(f"chainID {chain}")
        if selection.n_atoms > 0:
            cog = selection.center_of_geometry()
            ax.text(cog[0], cog[1], "LHCB3", color='black', fontsize=8,
                   ha='center', va='center', zorder=20,
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="None", alpha=0))


def plot_helices(ax, u0, yaml_path, chains_to_plot=None, helices_to_plot=None, color='green'):
    """Plot specific helices from selected chains.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to plot on
    u0 : mda.Universe
        Reference PSII Universe
    yaml_path : str
        Path to YAML file with helix definitions
    chains_to_plot : list, optional
        List of chain IDs to plot helices for. Default is ['n', 'g', '6', '7', '8', 's']
    helices_to_plot : list, optional
        List of helix names to plot. Default is ['B', 'C']
    color : str, optional
        Color for the helices. Default is 'green'
    """
    if yaml_path is None:
        return
    
    # Default chains and helices
    if chains_to_plot is None:
        chains_to_plot = ["1","2","5","6","g","n","y","Y","G","N", '7', '8', 's','6', '3', '4', 'S',"R","r"]
    if helices_to_plot is None:
        helices_to_plot = ['C','A']
        colors = ['green', 'blue']

    # Load the YAML file
    with open(yaml_path, 'r') as f:
        helix_definitions = yaml.safe_load(f)
    
    print(f"\n=== Plotting Helices ===")
    
    for chain in chains_to_plot:
        chain_key = f"chain_{chain}"
        
        if chain_key not in helix_definitions:
            print(f"Warning: No helix definitions found for chain {chain}")
            continue
        
        for helix_name in helices_to_plot:
            if helix_name not in helix_definitions[chain_key]:
                print(f"Warning: Helix {helix_name} not found for chain {chain}")
                continue
            
            helix_info = helix_definitions[chain_key][helix_name]
            start_resid = helix_info['start']
            end_resid = helix_info['end']
            color = colors[helices_to_plot.index(helix_name)] if helix_name in helices_to_plot else color
            
            # Select atoms in the helix
            selection = u0.select_atoms(f"chainID {chain} and resid {start_resid}:{end_resid}")
            
            if selection.n_atoms > 0:
                # Plot helix atoms with zorder lower than residue data points (zorder=11)
                ax.scatter(selection.positions[:, 0], selection.positions[:, 1], 
                          c=color, alpha=0.7, s=60, zorder=10, edgecolors=color, linewidths=0.3)
                print(f"Plotted helix {helix_name} of chain {chain}: residues {start_resid}-{end_resid} ({selection.n_atoms} atoms)")
            else:
                print(f"Warning: No atoms found for helix {helix_name} of chain {chain}")
    
    print("=" * 50 + "\n")


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
    # Reset index to ensure proper alignment
    df = df.reset_index(drop=True)
    df['x'] = pos_df['x'].values
    df['y'] = pos_df['y'].values
    
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


def plot_residues_by_lifetime(ax, df, cmap_lifetime, vmin, vmax, log_transform=False, show_labels=True):
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
        Minimum value for color mapping (in the final display scale - log scale if log_transform=True)
    vmax : float
        Maximum value for color mapping (in the final display scale - log scale if log_transform=True)
    log_transform : bool, optional
        Whether to apply log transformation to lifetime values
        
    Returns
    -------
    tuple
        (actual_vmin, actual_vmax) - the actual min/max values used for coloring
    """
    # Remove rows with NaN positions
    df_valid = df.dropna(subset=['x', 'y'])
    
    # Sort by sum_ns so that zero/low values are plotted first
    df_valid = df_valid.sort_values('sum_ns')
    
    # Get lifetime values
    lifetimes = df_valid['sum_ns'].values
    
    print(f"\n=== Plot Color Scale ===")
    print(f"Original data range: Min={lifetimes.min():.2f}, Max={lifetimes.max():.2f}")
    
    # Apply log transformation if requested
    if log_transform:
        lifetimes = np.log10(lifetimes + 1)  # log10(x + 1)
        print(f"Data range (after log transform): Min={lifetimes.min():.2f}, Max={lifetimes.max():.2f}")
        # Transform vmin/vmax to log scale as well (user provides values in original scale)
        display_vmin = np.log10(vmin + 1) if vmin is not None else lifetimes.min()
        display_vmax = np.log10(vmax + 1) if vmax is not None else lifetimes.max()
        print(f"vmin/vmax transformed: {vmin} -> {display_vmin:.2f}, {vmax} -> {display_vmax:.2f}")
    else:
        # Use provided vmin/vmax or data range
        display_vmin = vmin if vmin is not None else lifetimes.min()
        display_vmax = vmax if vmax is not None else lifetimes.max()
    
    print(f"Color scale being used (vmin, vmax): ({display_vmin:.2f}, {display_vmax:.2f})")
    print("=" * 50 + "\n")
    
    # Normalize lifetime values
    norm = plt.Normalize(vmin=display_vmin, vmax=display_vmax)
    colors = cmap_lifetime(norm(lifetimes))
    
    # Plot residues with colored fill and black border
    ax.scatter(df_valid['x'].values, df_valid['y'].values, 
              c=colors, edgecolors='black', linewidths=0.5, 
              alpha=1, s=40, zorder=11)
    
    # Add labels to the right of residues with lifetime > 0
    if show_labels:
        for idx, row in df_valid.iterrows():
                if row['sum_ns'] > 0:
                    label = f"{row['resid']}{row['resname']}"
                    ax.text(row['x'] +2, row['y'], label, fontsize=5, 
                           ha='left', va='center', zorder=12)
    
    print(f"Plotted {len(df_valid)} residues with lifetimes")
    
    return display_vmin, display_vmax


def add_colorbar(ax, cmap_lifetime, vmin, vmax, original_vmin=None, original_vmax=None, label='Lifetime', log_transform=False):
    """Add lifetime colorbar to plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to attach colorbar to
    cmap_lifetime : matplotlib.colors.Colormap
        Colormap for lifetime
    vmin : float
        Minimum value for colorbar (already log-transformed if log_transform=True)
    vmax : float
        Maximum value for colorbar (already log-transformed if log_transform=True)
    original_vmin : float, optional
        Original (non-transformed) minimum value provided by user
    original_vmax : float, optional
        Original (non-transformed) maximum value provided by user
    label : str, optional
        Label for colorbar. Default is 'Lifetime'
    log_transform : bool, optional
        Whether log transformation was applied (for label only)
    """
    # vmin and vmax are already in the correct scale (log-transformed if needed)
    # from plot_residues_by_lifetime, so use them directly
    display_vmin = vmin if vmin is not None else 0
    display_vmax = vmax if vmax is not None else 1
    
    # Convert cmap library Colormap to matplotlib colormap
    mpl_cmap = cmap_lifetime.to_mpl() if hasattr(cmap_lifetime, "to_mpl") else cmap_lifetime
    
    sm = plt.cm.ScalarMappable(cmap=mpl_cmap, 
                              norm=plt.Normalize(vmin=display_vmin, vmax=display_vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    
    # If log transform was applied, show original values on colorbar ticks
    if log_transform:
        # Generate tick positions at powers of 10 (1, 10, 100, 1000, ...)
        # Start from 0 (which represents original value 0) and go up in integer steps
        tick_positions = list(np.arange(np.floor(display_vmin), np.floor(display_vmax) + 1))
        
        # Always include the exact vmax position if user provided original_vmax
        if original_vmax is not None and display_vmax not in tick_positions:
            tick_positions.append(display_vmax)
        
        tick_positions = sorted(tick_positions)
        
        # Create labels: 0, 1, 10, 100, 1000, ... and user's vmax
        tick_labels = []
        for t in tick_positions:
            if t == 0:
                tick_labels.append('0')
            elif t == display_vmax and original_vmax is not None:
                # Use the exact user-provided vmax value
                tick_labels.append(f'{int(original_vmax):,}')
            else:
                # Powers of 10: 10^1=10, 10^2=100, etc.
                tick_labels.append(f'{int(10**t):,}')
        
        cbar.set_ticks(tick_positions)
        cbar.set_ticklabels(tick_labels)
        label = f'{label} (ns)'
    else:
        # Linear scale - use nice round tick values
        # Generate ~5-6 evenly spaced ticks
        n_ticks = 5
        tick_positions = np.linspace(display_vmin, display_vmax, n_ticks)
        tick_labels = [f'{int(t):,}' for t in tick_positions]
        cbar.set_ticks(tick_positions)
        cbar.set_ticklabels(tick_labels)
        label = f'{label} (ns)'
    
    cbar.set_label(label, fontsize=20)
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


def save_figure(fig, output_path):
    """Save figure to disk as PNG.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save
    output_path : str
        Output file path
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    fig.savefig(output_path, dpi=600, bbox_inches='tight', format='png')
    
    print(f"Figure saved in {output_path}")


def save_pdb_with_bfactors(u0, df, output_path, log_transform=False):
    """Save PDB with lifetime values as B-factors.
    
    Parameters
    ----------
    u0 : mda.Universe
        Reference PSII Universe
    df : pd.DataFrame
        DataFrame with columns: resid, resname, chain, sum_ns
    output_path : str
        Output file path (will change extension to .pdb)
    log_transform : bool
        Whether to apply log transformation to lifetime values
    """
    # Create PDB output path
    base, _ = os.path.splitext(output_path)
    pdb_output_path = f"{base}.pdb"
    
    # Set all B-factors to 0 initially
    u0.atoms.tempfactors = 0.0
    
    # Create a lookup dictionary from dataframe
    # Key: (resid, chain), Value: sum_ns
    lifetime_dict = {}
    for _, row in df.iterrows():
        key = (int(row['resid']), str(row['chain']))
        lifetime_dict[key] = row['sum_ns']
    
    # Assign B-factors based on lifetime values
    for atom in u0.atoms:
        key = (atom.resid, atom.chainID)
        if key in lifetime_dict:
            value = lifetime_dict[key]
            if log_transform:
                value = np.log10(value + 1)
            atom.tempfactor = value
    
    # Save PDB
    u0.atoms.write(pdb_output_path)
    print(f"PDB with B-factors saved in {pdb_output_path}")


def main():
    """Main execution function."""
    # Parse arguments
    args = parse_arguments()
        
    # Load data
    u0 = load_reference_structure(args.ref)
    df = load_csv(args.csv)

    # Change chain IDs BEFORE adding positions (so positions are looked up with new chain IDs)
    df = change_chains(df, args.change_chains_yaml)

    # Duplicate entries for equivalent chains BEFORE adding positions
    df = duplicate_entries_for_equivalent_chains(df, args.equivalent_chains_yaml)

    # Print data statistics before processing
    print(f"\n=== Data Statistics (Before Processing) ===")
    if 'sum_ns' in df.columns:
        print(f"sum_ns column - Min: {df['sum_ns'].min():.2f}, Max: {df['sum_ns'].max():.2f}")


    # Add positions from structure (will use the updated chain IDs)
    df = add_positions_to_dataframe(df, u0)

    # Change residue names AFTER adding positions (for display purposes only)
    df = change_resnames(df, args.change_resnames_yaml)
    
    # Print data statistics after adding positions
    print(f"\n=== Data Statistics (After Adding Positions) ===")
    df_valid = df.dropna(subset=['x', 'y'])
    if 'sum_ns' in df_valid.columns:
        print(f"sum_ns column - Min: {df_valid['sum_ns'].min():.2f}, Max: {df_valid['sum_ns'].max():.2f}")
    
    # Setup figure
    cmap_ref_chains = plt.get_cmap("Set3")
    cmap_sites = Colormap(args.cmap)
    
     # Create figure and axis
    
    fig, ax = setup_figure_and_axis()
    
    # Plot components
    plot_psii_background(ax, u0)
    plot_reference_chains(ax, u0, cmap_ref_chains)
    plot_lhcii_complexes(ax, u0)
    
    # Plot helices B and C behind residue data points
    plot_helices(ax, u0, args.helix_labels_yaml)

    
    # Plot residues colored by lifetime
    actual_vmin, actual_vmax = plot_residues_by_lifetime(ax, df, cmap_sites, args.vmin, args.vmax, args.log_transform)
    
    # Finalize plot
    add_colorbar(ax, cmap_sites, vmin=actual_vmin, vmax=actual_vmax,
                original_vmin=args.vmin, original_vmax=args.vmax,
                label='Lifetime', log_transform=args.log_transform)
    format_plot(ax)
    
    # Save figure
    save_figure(fig, args.output)

    # Save figure without labels
    base, ext = os.path.splitext(args.output)
    output_nolabels = f"{base}_no_labels{ext}"
    
    # Create new figure for no-labels version
    fig2, ax2 = setup_figure_and_axis()
    
    # Plot same components without labels
    plot_psii_background(ax2, u0)
    plot_reference_chains(ax2, u0, cmap_ref_chains)
    plot_lhcii_complexes(ax2, u0)
    plot_helices(ax2, u0, args.helix_labels_yaml)
    plot_residues_by_lifetime(ax2, df, cmap_sites, args.vmin, args.vmax, args.log_transform, show_labels=False)
    add_colorbar(ax2, cmap_sites, vmin=actual_vmin, vmax=actual_vmax,
                original_vmin=args.vmin, original_vmax=args.vmax,
                label='Lifetime', log_transform=args.log_transform)
    format_plot(ax2)
    save_figure(fig2, output_nolabels)
    
    # Determine selected chains (chains that have residue data)
    selected_chains = set(df.dropna(subset=['x', 'y'])['chain'].unique())
    print(f"Selected chains (with data): {selected_chains}")

    # Save figure with white non-selected chains (suffix _white)
    output_white = f"{base}_white{ext}"
    fig3, ax3 = setup_figure_and_axis()
    plot_psii_background(ax3, u0, selected_chains=selected_chains, mode='white')
    plot_reference_chains(ax3, u0, cmap_ref_chains, selected_chains=selected_chains, mode='white')
    plot_lhcii_complexes(ax3, u0, selected_chains=selected_chains, mode='white')
    plot_helices(ax3, u0, args.helix_labels_yaml)
    plot_residues_by_lifetime(ax3, df, cmap_sites, args.vmin, args.vmax, args.log_transform, show_labels=False)
    add_colorbar(ax3, cmap_sites, vmin=actual_vmin, vmax=actual_vmax,
                original_vmin=args.vmin, original_vmax=args.vmax,
                label='Lifetime', log_transform=args.log_transform)
    format_plot(ax3)
    save_figure(fig3, output_white)

    # Save figure with only selected chains (suffix _only_chains)
    output_only = f"{base}_only_chains{ext}"
    fig4, ax4 = setup_figure_and_axis()
    plot_psii_background(ax4, u0, selected_chains=selected_chains, mode='only_chains')
    plot_reference_chains(ax4, u0, cmap_ref_chains, selected_chains=selected_chains, mode='only_chains')
    plot_lhcii_complexes(ax4, u0, selected_chains=selected_chains, mode='only_chains')
    plot_helices(ax4, u0, args.helix_labels_yaml)
    plot_residues_by_lifetime(ax4, df, cmap_sites, args.vmin, args.vmax, args.log_transform, show_labels=False)
    add_colorbar(ax4, cmap_sites, vmin=actual_vmin, vmax=actual_vmax,
                original_vmin=args.vmin, original_vmax=args.vmax,
                label='Lifetime', log_transform=args.log_transform)
    format_plot(ax4)
    save_figure(fig4, output_only)

    # Save PDB with B-factors if log_transform is specified
    if args.log_transform:
        save_pdb_with_bfactors(u0, df, args.output, log_transform=args.log_transform)


if __name__ == "__main__":
    main()


