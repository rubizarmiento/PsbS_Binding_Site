"""
Author: Rubi Zarmiento-Garcia
Date: 20/08/2025

Generate PDB files with the binding poses of two atomsgroups.
Uses gmx cluster

Workflow:
1. Reads a trajectory and two selections. The first selection is the mobile group and the second selection is the binding site.
2. Filters the frames where the selections are in contact.
3. Extract clusters using the gromos clustering algorithm using different cut-off.
4. The user should identify the best cluster representative.

Arguments:
-f: Path to the structure file, compatible with MDAnalysis.
-trj: Path to the trajectory file, compatible with MDAnalysis.
-tpr: TPR file, compatible with MDAnalysis.
-sel1: Selection string for the first atom group.
-sel2: Selection string for the second atom group.
-filter: Distance filter for the selections (default: 0.8 nm).
-o: Output directory for the PDB files.
-cutoff: Cut-off distance for clustering in Angstroms (default: [0.15,0.30,0.45])

# Example usage:
python binding_pose.py -f chain_4/initial_fit.pdb -trj chain_4/test.xtc -tpr chain_4/protein.tpr -sel1 "chainID 4" -sel2 "chainID A B" -o binding_poses_main_test -filter 0.8 --cutoff 0.15 0.30 0.45
"""
import argparse
import os
import subprocess

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

def parse_args():
    parser = argparse.ArgumentParser(description="Generate PDB files with binding poses of two atom groups.")
    parser.add_argument('-f', '--structure', required=True, help='Path to the structure file (PDB, GRO, etc.)')
    parser.add_argument('-trj', '--trajectory', required=True, help='Path to the trajectory file (DCD, XTC, etc.)')
    parser.add_argument('-tpr', '--tpr', required=True, help='TPR file, compatible with MDAnalysis')
    parser.add_argument('-sel1', '--selection1', required=True, help='Selection string for the first atom group')
    parser.add_argument('-sel2', '--selection2', required=True, help='Selection string for the second atom group')
    parser.add_argument('-o', '--output_dir', default='binding_poses', help='Output directory for PDB files')
    parser.add_argument('-filter', '--distance_filter', type=float, default=0.8, help='Distance filter for the selections in nm')
    parser.add_argument('--cutoff', nargs='+', type=float, default=[0.15, 0.30, 0.45], help='Cut-off distance for clustering in nanometers')
    
    return parser.parse_args()

def filter_frames_by_distance(u, sel1, sel2, distance_filter):
    """
    Filter frames where the distance between selections is less than the distance_filter.
    
    Parameters:
        distance_filter: Distance filter in nanometers
    """
    filtered_frames = []
    distance_filter_angstrom = distance_filter * 10  # Convert nm to Angstrom
    
    for ts in u.trajectory:
        dist_matrix = distance_array(sel1.positions, sel2.positions)
        min_dist = np.min(dist_matrix)
        if min_dist < distance_filter_angstrom:
            filtered_frames.append(ts.frame)
    
    print(f"Distance filter: {distance_filter} nm ({distance_filter_angstrom} Å)")
    return filtered_frames
def write_ndx_file(ndx_file, selections):
    """
    Write a GROMACS index file with the selections.

    Parameters:
        ndx_file (str): Path to the output index file.
        selections (list): List of two MDAnalysis AtomGroup objects, where
            selections[0] is written as [ group1 ] and selections[1] as [ group2 ].
    """
    with open(ndx_file, 'w') as f:
        f.write("[ group1 ]\n")
        # GROMACS indices are 1-based, so add 1 to MDAnalysis 0-based indices
        f.write(" ".join(map(str, selections[0].indices + 1)) + "\n")
        f.write("[ group2 ]\n")
        f.write(" ".join(map(str, selections[1].indices + 1)) + "\n")
        f.write("[ combined ]\n")
        # Combine both groups for output
        combined_indices = np.concatenate([selections[0].indices, selections[1].indices])
        f.write(" ".join(map(str, combined_indices + 1)) + "\n")


def extract_clusters(trj, tpr, ndx, cutoff,out_path, group1_name="group1", group2_name="group2"):
    """
    Extract clusters from the filtered frames using the Gromos clustering algorithm.

    Side effects:
        - Writes clustering output files (e.g., cluster-id.xvg, sizes.xvg, centers.pdb, rmsd-clust.xpm, rmsd-dist.xvg)
          into subdirectories named according to the cutoff value.
        - Executes the 'gmx cluster' command via subprocess for each cutoff value.
    """
    # Store current working directory
    original_cwd = os.getcwd()
    
    for c in cutoff:
        outdir = f"{out_path}/clust_c{str(c).replace('.','')}"
        os.makedirs(outdir, exist_ok=True)

        # Convert relative paths to absolute paths
        abs_tpr = os.path.abspath(tpr)
        abs_trj = os.path.abspath(trj)
        abs_ndx = os.path.abspath(ndx)

        # Build the gmx command with relative paths (since we'll be in outdir)
        cmd = [
            "gmx", "cluster",
            "-s", abs_tpr,
            "-f", abs_trj,
            "-n", abs_ndx,
            "-method", "gromos",
            "-cutoff", str(c),
            "-clid", "cluster-id.xvg",
            "-sz", "sizes.xvg",
            "-cl", "centers.pdb",
            "-o", "rmsd-clust.xpm",
            "-dist", "rmsd-dist.xvg"
        ]

        # Change to output directory so log files are written there
        os.chdir(outdir)
        
        # Pipe group names (fit group and output group)
        # Use group1 for fitting, combined for output to include both selections
        result = subprocess.run(
            cmd,
            input=f"{group1_name}\ncombined\n",
            capture_output=True,
            text=True
        )
        
        # Change back to original directory
        os.chdir(original_cwd)
        
        # Print the command being executed
        print(f"Running command: {' '.join(cmd)}")
        print(f"Fit group: {group1_name}, Output group: combined")
        
        # Print stdout and stderr for debugging
        if result.stdout:
            print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
            
        if result.returncode != 0:
            print(f"Error running gmx cluster for cutoff {c}. Return code: {result.returncode}")
        else:
            print(f"Clustering completed for cutoff {c}. Output files are in {outdir}.")

        

def main():
    args = parse_args()

    # Load the universe
    u = mda.Universe(args.structure, args.trajectory)

    # Select the atom groups
    sel1 = u.select_atoms(args.selection1)
    sel2 = u.select_atoms(args.selection2)

    if len(sel1) == 0 or len(sel2) == 0:
        raise ValueError("One of the selections is empty. Please check your selection strings.")

    # Filter frames based on distance
    filtered_frames = filter_frames_by_distance(u, sel1, sel2, args.distance_filter)
    print(f"Filtered frames: {len(filtered_frames)}")

    # Save the filtered trajectory
    filtered_traj_path = os.path.join(args.output_dir, "filtered_frames.xtc")
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Writing filtered trajectory to: {filtered_traj_path}")
    print(f"Total atoms: {u.atoms.n_atoms}")
    
    with mda.Writer(filtered_traj_path, u.atoms.n_atoms) as W:
        frames_written = 0
        for ts in u.trajectory:
            if ts.frame in filtered_frames:
                W.write(u.atoms)
                frames_written += 1
        print(f"Frames written to trajectory: {frames_written}")
    
    # Check if file was created and has content
    if os.path.exists(filtered_traj_path):
        size = os.path.getsize(filtered_traj_path)
        print(f"Filtered trajectory file size: {size} bytes")
    else:
        print("ERROR: Filtered trajectory file was not created!")

    if not filtered_frames:
        print("No frames found where the selections are within the distance filter.")
        return
    else:
        print(f"Filtered {len(filtered_frames)} frames where the selections are within the distance filter.")

    # Write index file for GROMACS
    ndx_file = os.path.join(args.output_dir, 'selections.ndx')
    write_ndx_file(ndx_file, [sel1, sel2])
    print(f"Index file written to {ndx_file}")
    
    # Extract clusters using Gromos clustering algorithm on filtered trajectory
    extract_clusters(filtered_traj_path, args.tpr, ndx_file, args.cutoff, args.output_dir, group1_name="group1", group2_name="group2")


def test():
    """
    Test function to ensure the script runs without errors.
    """
    tpr = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/protein.tpr"
    structure = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/initial_fit.pdb"
    trajectory = '/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4/test.xtc'
    selection2 = 'chainID A B'
    selection1 = 'chainID 4'
    args = argparse.Namespace()
    args.structure = structure
    args.trajectory = trajectory
    args.tpr = tpr
    args.selection1 = selection1
    args.selection2 = selection2
    args.output_dir = 'binding_poses_test'
    args.distance_filter = 0.8
    args.cutoff = [0.15, 0.30, 0.45]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # Load the universe and selections
    u = mda.Universe(args.structure, args.trajectory)
    sel1 = u.select_atoms(args.selection1)
    sel2 = u.select_atoms(args.selection2)
    if len(sel1) == 0 or len(sel2) == 0:
        raise ValueError("One of the selections is empty. Please check your selection strings.")
    # Filter frames based on distance
    # Filter frames based on distance
    filtered_frames = filter_frames_by_distance(u, sel1, sel2, args.distance_filter)
    print(f"Filtered frames: {len(filtered_frames)}")
    # Save the filtered trajectory
    filtered_traj_path = os.path.join(args.output_dir, "filtered_frames.xtc")
    if not os.path.exists(filtered_traj_path):
        with mda.Writer(filtered_traj_path, u.atoms.n_atoms) as W:
            for ts in u.trajectory:
                if ts.frame in filtered_frames:
                    W.write(u.atoms)
    if not filtered_frames:
        print("No frames found where the selections are within the distance filter.")
        exit(0)
    else:
        print(f"Filtered {len(filtered_frames)} frames where the selections are within the distance filter.")
    # Write index file for GROMACS
    ndx_file = os.path.join(args.output_dir, 'selections.ndx')
    write_ndx_file(ndx_file, [sel1, sel2])
    print(f"Index file written to {ndx_file}")
    # Extract clusters using Gromos clustering algorithm on filtered trajectory
    extract_clusters(filtered_traj_path, args.tpr, ndx_file, args.cutoff, args.output_dir, group1_name="group1", group2_name="group2")



if __name__ == "__main__":
    main()

#test()
#main()