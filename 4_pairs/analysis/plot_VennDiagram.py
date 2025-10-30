"""
Binding Site Occupancy Visualization

Author: Rubi Zarmiento-Garcia

This module creates a visualization of protein binding site occupancy using circles 
and convex hulls. Each circle represents a protein chain, and convex hulls group 
chains that form binding sites together, with colors indicating occupancy levels.

Workflow:
    1. Load PDB structure with chain information
    2. Load occupancy data from CSV (trajectory, frames, normalized_occupancy, occupancy_percent)
    3. Calculate radius and center of geometry (x0, y0) per chain
    4. Group chains by trajectory (format: {id}_{chains})
    5. Map chains to colors and labels
    6. Plot circles for each chain using radius, position, and chain-specific colors
    7. Draw convex hulls around groups of chains, colored by occupancy level
"""

import csv
import math
import os
import warnings

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from matplotlib.patches import Circle
from scipy.spatial import ConvexHull

# Ignore MDAnalysis warnings
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')


def check_if_file_exists(file):
    """
    Check if a file exists at the specified path.
    
    Parameters
    ----------
    file : str
        Path to the file to check.
        
    Returns
    -------
    bool
        True if file exists.
        
    Raises
    ------
    SystemExit
        If file does not exist.
    """
    if not os.path.exists(file):
        print(f"Error: File {file} does not exist.")
        exit()
    return True


def get_universe(f):
    """
    Load a molecular structure using MDAnalysis.
    
    Parameters
    ----------
    f : str
        Path to structure file (PDB, GRO, etc.).
        
    Returns
    -------
    MDAnalysis.Universe
        Universe object containing the structure.
    """
    return mda.Universe(f)

def csv_to_dict(csv_file):
    """
    Read CSV file and convert to dictionary with columns as keys.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV file.
        
    Returns
    -------
    dict
        Dictionary with column headers as keys and lists of column values as values.
        
    Examples
    --------
    >>> data = csv_to_dict("occupancy.csv")
    >>> data.keys()
    dict_keys(['trajectory', 'frames', 'normalized_occupancy', 'occupancy_percent'])
    """
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        data = {}
        for row in reader:
            for key, value in row.items():
                if key not in data:
                    data[key] = []
                data[key].append(value)
    return data

def get_r_x0_y0(u):
    """
    Calculate radius and center of geometry for each chain.
    
    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe containing the molecular structure.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'chains' : list of str
            Chain IDs
        - 'radius' : list of float
            Radius for each chain (Angstroms)
        - 'x0' : list of float
            X coordinate of center of geometry (Angstroms)
        - 'y0' : list of float
            Y coordinate of center of geometry (Angstroms)
    """
    unique_chains = set(atom.chainID for atom in u.atoms)
    x0 = []
    y0 = []
    radius = []
    for chain in unique_chains:
        sel = u.select_atoms(f"chainID {chain} and name BB")
        cog = sel.center_of_geometry() # Returns x, y , z
        # Get radius 
        dx = max(sel.positions[:,0]) - min(sel.positions[:,0])
        dy = max(sel.positions[:,1]) - min(sel.positions[:,1])
        r = math.sqrt((dx/2)**2 + (dy/2)**2)
        x0.append(cog[0])
        y0.append(cog[1])
        radius.append(r)
    dict_chains_r_x0_y0 = {
        "chains": list(unique_chains),
        "radius": radius,
        "x0": x0,
        "y0": y0
    }

    return dict_chains_r_x0_y0

def get_grouped_chains_r_x0_y0(csv_dict, dict_chains_r_x0_y0):
    """
    Group chains by trajectory and retrieve their geometric properties.
    
    Parameters
    ----------
    csv_dict : dict
        Dictionary from CSV containing 'trajectory' key with values like "1_8_s_c".
    dict_chains_r_x0_y0 : dict
        Dictionary with keys 'chains', 'radius', 'x0', 'y0' containing geometric 
        properties for individual chains.
        
    Returns
    -------
    dict
        Dictionary with grouped chains and their properties:
        - 'chains' : list of list of str
            Grouped chain IDs (e.g., [['8', 's', 'c'], ['7', '6', 'g']])
        - 'radius' : list of list of float
            Radii for each chain in each group
        - 'x0' : list of list of float
            X coordinates for each chain in each group
        - 'y0' : list of list of float
            Y coordinates for each chain in each group
    """
    trajectory = csv_dict['trajectory']
    grouped_chains = [ trj.split('_')[1:] for trj in trajectory ] # e.g [['8', 's', 'c'], ['7', '6', 'g', 'n']]

    x0_arr_of_arr = []
    y0_arr_of_arr = []
    r_arr_of_arr = []
    for group in grouped_chains:
        temp_r = []
        temp_x0 = []
        temp_y0 = []
        for chain in group:
            index = dict_chains_r_x0_y0["chains"].index(chain)
            r = dict_chains_r_x0_y0["radius"][index]  
            x0 = dict_chains_r_x0_y0["x0"][index]             
            y0 = dict_chains_r_x0_y0["y0"][index]    
            temp_r.append(r)
            temp_x0.append(x0)
            temp_y0.append(y0)
        r_arr_of_arr.append(temp_r)
        x0_arr_of_arr.append(temp_x0)
        y0_arr_of_arr.append(temp_y0)

    dict_grouped_chains_r_x0_y0 = {
        "chains": list(grouped_chains),
        "radius": r_arr_of_arr,
        "x0": x0_arr_of_arr,
        "y0": y0_arr_of_arr
    }
    return dict_grouped_chains_r_x0_y0

def get_dict_chains_occupancy(csv_dict):
    """
    Extract occupancy data grouped by chain.
    
    Parameters
    ----------
    csv_dict : dict
        Dictionary from CSV with 'trajectory' and 'occupancy_percent' keys.
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'chains' : list of list of str
            Chain IDs extracted from trajectory names
        - 'occupancy' : list of str or float
            Occupancy percentages for each trajectory
    """
    occupancy = csv_dict["occupancy_percent"]
    trajectory = csv_dict["trajectory"]
    chains = []
    for chain in trajectory:
        arr = chain.split('_')
        arr = arr[1:]  # Skip trajectory ID, keep only chain names
        chains.append(arr)
    dict_chains_occupancy = {
        "chains": chains,
        "occupancy": occupancy
    }
    return dict_chains_occupancy

def generate_circle_points(x0, y0, radius, n_points=20):
    """
    Generate points along the circumference of a circle.
    
    Args:
        x0 (float): Center x coordinate
        y0 (float): Center y coordinate  
        radius (float): Circle radius
        n_points (int): Number of points to generate around circumference
        
    Returns:
        np.array: Array of shape (n_points, 2) with x, y coordinates
    """
    angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    x_points = x0 + radius * np.cos(angles)
    y_points = y0 + radius * np.sin(angles)
    return np.column_stack([x_points, y_points])

def plot_convex_hull_around_circles(x0_arr, y0_arr, r_arr, ax, n_points=20, 
                                   color='black', linestyle='--', linewidth=2, alpha=1.0):
    """
    Plot a convex hull that wraps around the circumference of a group of circles.
    
    Parameters
    ----------
    x0_arr : array-like
        X coordinates for circle centers.
    y0_arr : array-like
        Y coordinates for circle centers.
    r_arr : array-like
        Circle radii.
    ax : matplotlib.axes.Axes
        Matplotlib axis object to draw on.
    n_points : int, optional
        Number of points per circle circumference (default: 20).
    color : str, optional
        Color of the hull boundary (default: 'black').
    linestyle : str, optional
        Line style for hull boundary (default: '--').
    linewidth : float, optional
        Line width for hull boundary (default: 2).
    alpha : float, optional
        Transparency of hull boundary (default: 1.0).
        
    Returns
    -------
    scipy.spatial.ConvexHull or None
        ConvexHull object if successful, None if insufficient points.
    """
    # Collect all circumference points from all circles
    all_points = []
    for x0, y0, r in zip(x0_arr, y0_arr, r_arr):
        circle_points = generate_circle_points(x0, y0, r, n_points)
        all_points.append(circle_points)
    
    # Combine all points into a single array
    all_points = np.vstack(all_points)
    
    # Compute convex hull
    if len(all_points) >= 3:
        hull = ConvexHull(all_points)
        hull_points = all_points[hull.vertices]
        
        # Plot the hull as a polygon
        hull_polygon = plt.Polygon(hull_points, closed=True, fill=None, 
                                   edgecolor=color, linestyle=linestyle, linewidth=linewidth, alpha=alpha)
        ax.add_patch(hull_polygon)
        
        return hull
    else:
        print("Not enough points to create convex hull (need at least 3)")
        return None

def add_occupancy_to_dict(dict_grouped_chains_r_x0_y0, dict_chains_occupancy):
    """
    Add occupancy information to grouped chains dictionary.
    
    Parameters
    ----------
    dict_grouped_chains_r_x0_y0 : dict
        Dictionary with grouped chains and their geometric properties.
    dict_chains_occupancy : dict
        Dictionary with chains and their occupancy values.
        
    Returns
    -------
    dict
        Updated dictionary with 'occupancy' key added containing occupancy 
        values for each chain group.
    """
    occupancy_list = []
    for group in dict_grouped_chains_r_x0_y0["chains"]:
        index = dict_chains_occupancy["chains"].index(group)
        occupancy_value = dict_chains_occupancy["occupancy"][index]
        occupancy_list.append(occupancy_value)
    
    dict_grouped_chains_r_x0_y0["occupancy"] = occupancy_list
    return dict_grouped_chains_r_x0_y0

def plot_circle(dict_grouped_chains_r_x0_y0, dict_chains_color_labels):
    """
    Create visualization of binding sites with circles and convex hulls.
    
    Parameters
    ----------
    dict_grouped_chains_r_x0_y0 : dict
        Dictionary with grouped chains, their geometric properties, and occupancy.
    dict_chains_color_labels : dict
        Dictionary mapping chains to colors and labels.
        
    Returns
    -------
    tuple
        (fig, ax) - Matplotlib figure and axes objects.
    """
    # Get max occupancy value for normalization
    max_occupancy = max([float(val) for val in dict_grouped_chains_r_x0_y0["occupancy"]])

    # Sort by occupancy (ascending order)
    occupancy_list = [float(val) for val in dict_grouped_chains_r_x0_y0["occupancy"]]
    sorted_indices = sorted(range(len(occupancy_list)), key=lambda i: occupancy_list[i], reverse=False)
    for key in dict_grouped_chains_r_x0_y0:
        dict_grouped_chains_r_x0_y0[key] = [dict_grouped_chains_r_x0_y0[key][i] for i in sorted_indices]

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal', 'box')
    
    counter = 0  # Used to progressively increase hull linewidth
    cmap_rdpu = plt.get_cmap("RdPu")  # Sequential colormap for occupancy
    for group in dict_grouped_chains_r_x0_y0["chains"]:
        index_group = dict_grouped_chains_r_x0_y0["chains"].index(group)
        r_arr = dict_grouped_chains_r_x0_y0["radius"][index_group]
        x0_arr = dict_grouped_chains_r_x0_y0["x0"][index_group]
        y0_arr = dict_grouped_chains_r_x0_y0["y0"][index_group]
        # Special case: adjust chain e position to avoid overlap with chain f
        if 'e' in group:
            index_e = group.index('e')
            y0_arr[index_e] += 15
        
        occupancy_value = float(dict_grouped_chains_r_x0_y0["occupancy"][index_group])
        normalized_occupancy = occupancy_value / max_occupancy  # Normalize to [0, 1]
        color_occupancy = cmap_rdpu(normalized_occupancy)
        for chain in group:
            if chain not in dict_chains_color_labels["chains"]:
                color = "lightgray"
                uppercase_chain = chain.upper()
                label = f"PSB{uppercase_chain}"
            else:
                index_chain = dict_chains_color_labels["chains"].index(chain)
                color = dict_chains_color_labels["colors"][index_chain]
                label = dict_chains_color_labels["labels"][index_chain]
            index_r = group.index(chain)
            circle = Circle((x0_arr[index_r], y0_arr[index_r]), 15, 
                          color=color, alpha=1, label=label)
            ax.add_patch(circle)
            
            # Add label at circle center
            ax.text(x0_arr[index_r], y0_arr[index_r], label, 
                   fontsize=8, ha='center', va='center')

    for group in dict_grouped_chains_r_x0_y0["chains"]:
        index_group = dict_grouped_chains_r_x0_y0["chains"].index(group)
        r_arr = dict_grouped_chains_r_x0_y0["radius"][index_group]
        x0_arr = dict_grouped_chains_r_x0_y0["x0"][index_group]
        y0_arr = dict_grouped_chains_r_x0_y0["y0"][index_group]

        occupancy_value = float(dict_grouped_chains_r_x0_y0["occupancy"][index_group])
        normalized_occupancy = occupancy_value / max_occupancy
        color_occupancy = cmap_rdpu(normalized_occupancy)
        
        # Plot convex hull around the circle group
        n_chains = len(group)
        r_fixed = np.full(n_chains, 20)  # Fixed radius for hull calculation
        scale_factor = 0.5  # Linewidth increment per group
        linewidth = 2 + counter * scale_factor
        
        hull = plot_convex_hull_around_circles(
            x0_arr, y0_arr, r_fixed, ax, n_points=30, 
            color=color_occupancy, linestyle='-', linewidth=linewidth, alpha=0.7
        )
        counter += 1
    
    # Set plot limits and styling
    ax.set_xlim(150, 500)
    ax.set_ylim(100, 250)
    ax.axis('off')
    
    # Save figure
    output_path = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_site_occupancy.png"
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    return fig, ax

def main():
    """
    Main function to generate binding site occupancy visualization.
    
    Loads structure and occupancy data, processes chains, and creates
    a plot with circles representing chains and convex hulls showing
    binding sites colored by occupancy level.
    """
    # Input files
    f = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb"
    csv_file = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned/occupancy.csv"
    
    cmap = plt.get_cmap("Set3")
    colors_plot = [cmap(4), cmap(5),  cmap(7),"#05802c","#98E400","#98E400","#98E400", "#98E400"]  
    labels_plot = ["CP24","CP26", "CP43","LHCB3","LHCBM","LHCBM","LHCBM","LHCBM"]
    chain_sel = ["8", "s", "c","7","6","g","n","Y"]
    dict_chains_color_labels = {
        "chains": chain_sel,
        "colors": colors_plot,
        "labels": labels_plot
    }

    u = get_universe(f)
    csv_dict = csv_to_dict(csv_file)
    dict_chains_r_x0_y0 = get_r_x0_y0(u)

    # Save dict_chains_r_x0_y0 dict to a csv file
    with open("/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/dict_chains_r_x0_y0.csv", 'w', newline='') as csvfile:
        fieldnames = ['chain', 'radius', 'x0', 'y0']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for i, chain in enumerate(dict_chains_r_x0_y0["chains"]):
            writer.writerow({
                'chain': chain,
                'radius': dict_chains_r_x0_y0["radius"][i],
                'x0': dict_chains_r_x0_y0["x0"][i],
                'y0': dict_chains_r_x0_y0["y0"][i]
            })

    # Filter to only include chains that exist in both PDB and CSV
    pdb_chains = set(dict_chains_r_x0_y0["chains"])
    csv_chains = set()
    for trajectory in csv_dict["trajectory"]:
        chains = trajectory.split('_')[1:]
        csv_chains.update(chains)
    common_chains = sorted(list(pdb_chains.intersection(csv_chains)))
    print(f"Using common chains: {common_chains}")
    
    dict_chains_occupancy = get_dict_chains_occupancy(csv_dict)
    dict_grouped_chains_r_x0_y0 = get_grouped_chains_r_x0_y0(csv_dict,dict_chains_r_x0_y0)
    dict_grouped_chains_r_x0_y0 = add_occupancy_to_dict(dict_grouped_chains_r_x0_y0, dict_chains_occupancy)

    fig, ax = plot_circle(dict_grouped_chains_r_x0_y0, dict_chains_color_labels)



# Main function
main()



