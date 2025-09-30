"""
Author: Rubi Zarmiento-Garcia
Date 30/08/2025 

Workflow for Clustering Binding Modes from Molecular Dynamics Trajectories

This script performs leader-based clustering of binding events from MD trajectories.
The workflow consists of three main steps:

1. PREPARE EVENT REPRESENTATIVES (prepare_event_reps function):
   - Load topology files for each binding event using MDAnalysis
- Read the lifetime for each event from a np.npy file
   - Align all frames to a common reference using receptor atoms
   - Extract ligand coordinates for RMSD comparisons
   - Return aligned ligand coordinates for all events

2. CALCULATE LIGAND RMSD (ligand_rmsd function):
   - Compute RMSD between pairs of ligand coordinate sets
   - Support for permutation-based RMSD (for symmetric ligands)
   - Uses MDAnalysis RMSD with centering and superposition

3. PERFORM LEADER CLUSTERING (leader_cluster function):
   - Sort events by duration (longest-lived first)
   - Start with longest event as first cluster leader
   - For each subsequent event:
     * Calculate RMSD to all existing cluster leaders
     * Assign to closest cluster if RMSD < cutoff
     * Create new cluster if no close match found
   - Calculate cluster probabilities based on duration weights
   - Return cluster assignments, members, and statistics

Arguments:
    - f: Structure files
    - npy: Lifetimes .npy file
    - receptor_sel: Selection string for receptor atoms
    - ligand_sel: Selection string for ligand atoms
    - odir: Output directory for results
    - cutoff: RMSD cutoff in nm (default 0.34 nm)

Outputs:
    leaders_idx.csv: Indices of cluster leaders
    cluster_members.csv: List of members for each cluster
    prob.csv: Cluster probabilities based on duration weights
    assigning.csv: Cluster assignments for each event

Example:
python3 cluster_binding_modes.py -f f1.pdb f2.pdb -npy lifetime_file -receptor_sel "protein and name CA" -ligand_sel "resname LIG"

"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import argparse
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Cluster binding modes from MD trajectories")
    parser.add_argument("-f", "--files", nargs="+", required=True, help="Topology files for each event")
    parser.add_argument("-npy", "--lifetime_npy", required=True, help="Path to .npy file with event lifetimes")
    parser.add_argument("-receptor_sel", "--receptor_selection", required=True, help="MDAnalysis selection string for receptor atoms")
    parser.add_argument("-ligand_sel", "--ligand_selection", required=True, help="MDAnalysis selection string for ligand atoms")
    parser.add_argument("-odir", "--output_dir", required=True, help="Output directory for results")
    parser.add_argument("-cutoff", "--rmsd_cutoff", type=float, default=0.34, help="RMSD cutoff in nm (default 0.34 nm)")
    return parser.parse_args()


def prepare_event_reps(topology_files, lifetime_npy, receptor_sel, ligand_sel, ref_frame=None):
    """
    Load topology files for each binding event, read lifetimes from .npy file,
    align frames to common reference, and extract ligand coordinates.

    Parameters
    ----------
    topology_files : list[str]
        List of topology file paths (one per binding event)
    lifetime_npy : str
        Path to .npy file containing event lifetimes
    receptor_sel, ligand_sel : str
        MDAnalysis selection strings for alignment and ligand RMSD atoms
    ref_frame : int or None
        Frame index to serve as alignment reference

    Returns
    -------
    coords : (N, M, 3) float64
        Ligand coordinates for each event (Å)
    durations_ns : (N,) float64
        Event lifetimes in nanoseconds
    """
    # Read lifetimes from .npy file
    durations_ns = np.load(lifetime_npy)
    n_events = len(durations_ns)

    if len(topology_files) != n_events:
        raise ValueError(f"Number of topology files ({len(topology_files)}) must match number of events ({n_events})")

    coords_list = []
    uref = None

    for i, (top_file, duration) in enumerate(zip(topology_files, durations_ns)):
        # Load universe for this event
        u = mda.Universe(top_file)

        # Use first event as reference if not specified
        if ref_frame is None and uref is None:
            uref = mda.Universe(top_file)
            ref = uref.select_atoms(receptor_sel)

        # Align to reference
        if uref is not None:
            aligner = align.AlignTraj(u, uref, select=receptor_sel, weights="mass", in_memory=True)
            aligner.run()

        # Extract combined receptor and ligand coordinates (assuming single frame per topology)
        rec = u.select_atoms(receptor_sel)
        lig = u.select_atoms(ligand_sel)
        if rec.n_atoms == 0:
            raise ValueError(f"Receptor selection \"{receptor_sel}\" matched 0 atoms in {top_file}")
        if lig.n_atoms == 0:
            raise ValueError(f"Ligand selection \"{ligand_sel}\" matched 0 atoms in {top_file}")
        combined_coords = np.concatenate([rec.positions, lig.positions])
        coords_list.append(combined_coords.astype(np.float64).copy())

    return np.asarray(coords_list), durations_ns

def ligand_rmsd(A, B, perms=None):
    """
    RMSD between two ligand coordinate sets (Å). If 'perms' provided, return the
    minimum RMSD over all permutations.
    """
    if perms is None:
        # Simple RMSD without permutations
        return rmsd(A, B, center=True, superposition=True)
    else:
        # Find minimum RMSD over permutations
        best = np.inf
        for p in perms:
            if len(p) != len(B) or np.any(np.array(p) >= len(B)):
                raise ValueError("Permutation in 'perms' is invalid for the shape of B.")
            best = min(best, rmsd(A, B[p, :], center=True, superposition=True))
        return best

def leader_cluster(coords_Ang, durations_ns, cutoff_nm=0.34, perms=None):
    """
    Greedy (leader) clustering of event representatives, seeded by the longest-lived event.

    Parameters
    ----------
    coords_Ang : (N, M, 3)
        Ligand coordinates per event in Å (output of prepare_event_reps).
    durations_ns : (N,)
        Event lifetimes; used for seeding order (descending) and weighting probabilities.
    cutoff_nm : float
        Ligand RMSD cutoff in nm (will be converted to Å internally).
    perms : list or None
        Optional permutations for RMSD calculation.

    Returns
    -------
    dict
        {
          "leaders_idx": [int],            # event indices that define each mode
          "members":     [list[int]],      # event indices per mode
          "assign":      (N,) int,         # mode id for each event
          "prob":        (K,) float,       # duration-weighted probabilities per mode
          "dur_sum":     (K,) float,       # total duration per mode (ns)
          "cutoff_A":    float             # cutoff in Å actually used
        }
    """
    N = coords_Ang.shape[0]
    # Use durations as weights for probabilities
    frames_per_event = durations_ns.copy()
    order = np.argsort(durations_ns)[::-1]  # longest first
    cutoff_A = cutoff_nm * 10.0

    leaders_idx = []
    members = []
    totals_frames = []
    totals_dur = []
    assign = -np.ones(N, dtype=int)

    for k in order:
        x = coords_Ang[k]
        if not leaders_idx:
            # First leader
            leaders_idx.append(k)
            members.append([k])
            totals_frames.append(frames_per_event[k])
            totals_dur.append(durations_ns[k])
            assign[k] = 0
            continue

        # Distance to existing leaders
        dists = [ligand_rmsd(x, coords_Ang[i0], perms) for i0 in leaders_idx]
        i_min = int(np.argmin(dists))
        if dists[i_min] < cutoff_A:
            # Assign to existing cluster
            members[i_min].append(k)
            totals_frames[i_min] += frames_per_event[k]
            totals_dur[i_min] += durations_ns[k]
            assign[k] = i_min
        else:
            # New cluster
            leaders_idx.append(k)
            members.append([k])
            totals_frames.append(frames_per_event[k])
            totals_dur.append(durations_ns[k])
            assign[k] = len(leaders_idx) - 1

    totals_frames = np.asarray(totals_frames, dtype=float)
    total_sum = totals_frames.sum()
    if total_sum == 0:
        prob = np.zeros_like(totals_frames)
    else:
        prob = totals_frames / total_sum

    return {
        "leaders_idx": leaders_idx,
        "members": members,
        "assign": assign,
        "prob": prob,
        "dur_sum": np.asarray(totals_dur, dtype=float),
        "cutoff_A": cutoff_A,
    }

def write_results(results, output_dir):
    # Write each result array to a separate CSV file
    for key, value in results.items():
        if key in ["cutoff_A", "assign"]:
            continue  # Skip cutoff_A and assign
        if key == "members":
            # Write each cluster's members as a row with cluster_id
            with open(f"{output_dir}/{key}.csv", "w") as f:
                f.write("cluster_id\tmember_indices\n")
                for i, sublist in enumerate(value):
                    f.write(f"{i}\t" + "\t".join(map(str, [x+1 for x in sublist])) + "\n")
        elif key == "leaders_idx":
            np.savetxt(f"{output_dir}/{key}.csv", [x + 1 for x in value], delimiter="\t", fmt='%d', header="leader_index")
        elif key == "prob":
            np.savetxt(f"{output_dir}/{key}.csv", value, delimiter="\t", fmt='%.2f', header="probability")
        elif key == "dur_sum":
            np.savetxt(f"{output_dir}/{key}.csv", value, delimiter="\t", fmt='%.2f', header="duration_ns")
        else:
            np.savetxt(f"{output_dir}/{key}.csv", value, delimiter="\t")
    
    # Write README.md
    with open(f"{output_dir}/README.md", "w") as f:
        f.write("This directory contains clustering results:\n")
        f.write("- leaders_idx.csv: indices of cluster leaders (with header 'leader_index')\n")
        f.write("- members.csv: cluster members with cluster_id (with headers 'cluster_id' and 'member_indices')\n")
        f.write("- prob.csv: probability of each cluster (with header 'probability')\n")
        f.write("- dur_sum.csv: total duration of each cluster in ns (with header 'duration_ns')\n")
    
    print(f"Results written to {output_dir}")

def main():
    args = parse_arguments()

    # Prepare coordinates and read lifetimes
    coords, durations_ns = prepare_event_reps(args.files, args.lifetime_npy, args.receptor_selection, args.ligand_selection)

    # Perform clustering (uses durations as weights automatically)
    result = leader_cluster(coords, durations_ns, cutoff_nm=args.rmsd_cutoff)

    # Write results to output directory
    write_results(result, args.output_dir)


if __name__ == "__main__":
    main()