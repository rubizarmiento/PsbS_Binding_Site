"""
Author: Rubi Zarmiento-Garcia

Workflow:
    - Load PDB with chain information.
    - Load occupancy.csv with the headers trajectory,frames,normalized_occupancy,occupancy_percent i    cmap = plt.get_cmap("Set3")
    # Include ALL chains found in the data
    chain_sel = ['6', '7', '8', 'c', 'e', 'f', 'g', 'j', 'k', 'n', 'p', 's', 'z']
    colors_plot = [cmap(i % 12) for i in range(len(chain_sel))]  # Generate colors for all chains
    labels_plot = [f"Chain_{chain}" for chain in chain_sel]  # Generate labels for all chainsa dict.
    - Get dimensions (radius) and cog (x0, y0) per chain into a dict.
    - The column "trajectory" has the format {id}_{chains}, load the chains and "occupancy" into a dict.
    - Define dict chains and color.
    - Plot circles using the radius, y0, x0 and the color per chain.
    - Plot convex hull around the group of chains per "trajectory"
"""
import csv
import MDAnalysis as mda
import os
# Igonore warnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')
import math
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from matplotlib.patches import Circle
import numpy as np
def check_if_file_exists(file):
    if not os.path.exists(file):
        print(f"Error: File {file} does not exist.")
        exit()
    return True

def get_universe(f):
    return mda.Universe(f)



def csv_to_dict(csv_file):
    """
    Returns a dictionary from a CSV file where each key maps to a list of values from that column.
    Args:
        csv_file (str): Path to the CSV file.
    Returns:
        dict: Dictionary with column headers as keys and lists of column values.

    Example:
        {"trajectory": ["A1_8_s_c", "A2_7_6_g_n"],
         "frames": [1000, 1200],
         "normalized_occupancy": [0.8, 0.9],
         "occupancy_percent": [80, 90]}
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
    """Get radius, x0, y0 per chain in the universe u
    Args:
        u (MDAnalysis.Universe): Universe with the structure
    Returns:
        radius (list): List of radius per chain
        x0 (list): List of x0 per chain
        y0 (list): List of y0 per chain

    Units are in Angstroms
    
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

def get_grouped_chains_r_x0_y0(csv_dict,dict_chains_r_x0_y0):
    """
    Groups chains from the CSV dictionary and retrieves their corresponding radius, x0, and y0 values from the provided dictionary.

    Args:
        csv_dict (dict): Dictionary containing CSV data with a 'trajectory' key.
        dict_chains_r_x0_y0 (dict): Dictionary containing 'chains', 'radius', 'x0', and 'y0' lists.

    Returns:
        dict: Dictionary with grouped chains and their corresponding radius, x0, and y0 values.
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
    occupancy = csv_dict["occupancy_percent"]
    trajectory = csv_dict["trajectory"]
    chains = []
    for chain in trajectory:
        arr = chain.split('_')
        arr = arr[1:]
        chains.append(arr)
    dict_chains_occupancy = {
        "chains": chains,
        "occupancy": occupancy
    }
    return dict_chains_occupancy


#cmap_rdpu = plt.get_cmap("PuRd")
def plot_circle(dict_grouped_chains_r_x0_y0, dict_chains_color_labels):
    # Plot circles
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_title("Binding Site Occupancy")
    ax.set_aspect('equal', 'box')


    for group in dict_grouped_chains_r_x0_y0["chains"]:
        index_group = dict_grouped_chains_r_x0_y0["chains"].index(group)
        r_arr = dict_grouped_chains_r_x0_y0["radius"][index_group]
        x0_arr = dict_grouped_chains_r_x0_y0["x0"][index_group]
        y0_arr = dict_grouped_chains_r_x0_y0["y0"][index_group]
        for chain in group:
            print(chain)
            index_chain = dict_chains_color_labels["chains"].index(chain)
            color = dict_chains_color_labels["colors"][index_chain]
            label = dict_chains_color_labels["labels"][index_chain]
            index_r = group.index(chain)
            circle = Circle((x0_arr[index_r], y0_arr[index_r]), r_arr[index_r], color=color, alpha=0.5, label=label)
            ax.add_patch(circle)
    # Set limits
    ax.set_xlim(200, 500)
    ax.set_ylim(0, 350)
    # Save figure
    plt.savefig("/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_site_occupancy.png", dpi=300)
    return fig, ax

def main():
    # Test files
    f="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb"
    csv_file="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned/occupancy.csv"
    
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
    
    # Filter to only include chains that exist in both PDB and CSV
    pdb_chains = set(dict_chains_r_x0_y0["chains"])
    csv_chains = set()
    for trajectory in csv_dict["trajectory"]:
        chains = trajectory.split('_')[1:]
        csv_chains.update(chains)
    common_chains = sorted(list(pdb_chains.intersection(csv_chains)))
    print(f"Using common chains: {common_chains}")
    
    # Update the color dictionary to only include common chains
    dict_chains_color_labels["chains"] = common_chains[:len(dict_chains_color_labels["chains"])]
    
    dict_chains_occupancy = get_dict_chains_occupancy(csv_dict)
    dict_grouped_chains_r_x0_y0 = get_grouped_chains_r_x0_y0(csv_dict,dict_chains_r_x0_y0)
    fig, ax = plot_circle(dict_chains_r_x0_y0, dict_grouped_chains_r_x0_y0)



# Main function
main()



