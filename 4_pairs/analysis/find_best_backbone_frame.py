"""
Scan cluster member frames and find the one with the best backbone integrity
(smallest maximum BB-BB gap), then extract it as a replacement PDB.

The "middle" frame chosen by gmx cluster may have large backbone gaps (distorted CG
backbone), making cg2at fail. This script finds an alternative frame from the same
cluster that has the best backbone continuity.

Usage:
    python3 find_best_backbone_frame.py \
        -f structure.pdb -trj trajectory.xtc \
        -g cluster.log \
        -chains 7 8 \
        -o output.pdb \
        [-n 10] \
        [--max_gap_threshold 7.0] \
        [--exclude_chains 9]
"""

import MDAnalysis as mda
import numpy as np
import warnings
import argparse
import os
import sys

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find cluster frame with best backbone integrity"
    )
    parser.add_argument('-f', '--structure-file', type=str, required=True,
                        help='Path to the structure/topology file (.pdb)')
    parser.add_argument('-trj', '--trajectory-file', type=str, required=True,
                        help='Path to the trajectory file (.xtc)')
    parser.add_argument('-g', '--cluster-log', type=str, required=True,
                        help='Path to the cluster.log file from gmx cluster')
    parser.add_argument('-chains', nargs='+', required=True,
                        help='Protein chain IDs to check backbone integrity for')
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help='Path for the output PDB file')
    parser.add_argument('-n', '--n-frames', type=int, default=1,
                        help='Number of frames to extract (1=best only, >1 includes previous frames)')
    parser.add_argument('--max_gap_threshold', type=float, default=7.0,
                        help='Maximum allowed BB-BB distance (Angstrom) between consecutive residues. '
                             'Frames where all chains have max gap below this are considered "good". '
                             'Default: 7.0 A (normal CG BB spacing is ~3.5-4.5 A)')
    parser.add_argument('--exclude_chains', nargs='*', default=[],
                        help='Chain IDs to exclude from integrity check (e.g., PsbS dimer chain '
                             'that always has a large inter-protomer gap)')
    parser.add_argument('--cluster_id', type=int, default=1,
                        help='Cluster ID to scan (default: 1, the largest cluster)')
    return parser.parse_args()


def parse_cluster_log(log_file, cluster_id=1):
    """
    Parse gmx cluster log file and return member timestamps and middle timestamp
    for the specified cluster.
    
    Returns:
        middle_time (float): timestamp of the middle structure (ps)
        member_times (list of float): all member timestamps (ps)
    """
    with open(log_file, 'r') as f:
        lines = f.readlines()

    # Find header line
    header_idx = None
    for i, line in enumerate(lines):
        if 'cl. | #st  rmsd | middle rmsd | cluster members' in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find cluster log header")

    # Parse clusters
    clusters = {}
    current_cl = None
    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split('|')
        first = parts[0].strip()

        # New cluster line: starts with a digit and has >=4 pipe-separated fields
        if first and first[0].isdigit() and len(parts) >= 4:
            cl = int(first)
            # If we've already parsed our target cluster and hit a new one, stop
            if current_cl == cluster_id and cl != cluster_id:
                break
            current_cl = cl
            mid_parts = parts[2].strip().split()
            middle_time = float(mid_parts[0])
            member_str = parts[-1].strip()
            members = [float(x) for x in member_str.split() if x.strip()]
            clusters[cl] = {'middle': middle_time, 'members': members}
        elif current_cl is not None:
            # Continuation line with more member timestamps
            member_str = parts[-1].strip()
            more = [float(x) for x in member_str.split() if x.strip()]
            clusters[current_cl]['members'].extend(more)

    if cluster_id not in clusters:
        raise ValueError(f"Cluster {cluster_id} not found in log file")

    return clusters[cluster_id]['middle'], clusters[cluster_id]['members']


def calc_backbone_integrity(universe, frame_idx, chains):
    """
    Calculate backbone integrity metrics for a given frame.
    
    Returns:
        max_gap (float): maximum BB-BB gap across all specified chains
        total_gaps_gt10 (int): total number of gaps > 10 Å across all chains
        chain_details (dict): {chain_id: (max_gap, n_gaps_gt10)}
    """
    universe.trajectory[frame_idx]
    
    max_gap = 0.0
    total_gaps_gt10 = 0
    chain_details = {}

    for chain_id in chains:
        sel = universe.select_atoms(f'name BB and chainID {chain_id}')
        if len(sel) == 0:
            chain_details[chain_id] = (0.0, 0)
            continue
        positions = sel.positions
        distances = np.linalg.norm(np.diff(positions, axis=0), axis=1)
        chain_max = float(distances.max())
        n_gaps = int(np.sum(distances > 10))
        chain_details[chain_id] = (chain_max, n_gaps)

        if chain_max > max_gap:
            max_gap = chain_max
        total_gaps_gt10 += n_gaps

    return max_gap, total_gaps_gt10, chain_details


def extract_frames(universe, best_idx, output_file, n_frames=1):
    """Write the best frame (and optionally previous frames) to PDB files."""
    base_output = os.path.splitext(output_file)[0]

    # Write best frame
    universe.trajectory[best_idx]
    universe.atoms.write(output_file)
    print(f"  Wrote best frame (index {best_idx}, time {universe.trajectory.time:.0f} ps) -> {output_file}")

    if n_frames > 1:
        # Write previous frames
        start_frame = max(0, best_idx - n_frames + 1)
        counter = 1
        for frame_idx in range(start_frame, best_idx):
            universe.trajectory[frame_idx]
            out_path = f"{base_output}_prev_{counter}.pdb"
            universe.atoms.write(out_path)
            print(f"  Previous frame {counter}: index {frame_idx} (time {universe.trajectory.time:.0f} ps) -> {out_path}")
            counter += 1


def main():
    args = parse_args()

    # Parse cluster log
    middle_time, member_times = parse_cluster_log(args.cluster_log, args.cluster_id)
    print(f"Cluster {args.cluster_id}: {len(member_times)} member frames, middle at {middle_time:.0f} ps")

    # Load trajectory
    u = mda.Universe(args.structure_file, args.trajectory_file)
    timesteps = np.array([ts.time for ts in u.trajectory])
    print(f"Trajectory: {len(u.trajectory)} frames, {u.trajectory[0].time:.0f}-{u.trajectory[-1].time:.0f} ps")

    # Determine which chains to check
    check_chains = [c for c in args.chains if c not in args.exclude_chains]
    print(f"Checking backbone integrity for chains: {check_chains}")
    if args.exclude_chains:
        print(f"Excluding chains from integrity check: {args.exclude_chains}")

    # Check middle frame first
    middle_idx = int(np.argmin(np.abs(timesteps - middle_time)))
    mid_gap, mid_ngaps, mid_details = calc_backbone_integrity(u, middle_idx, check_chains)
    print(f"\nOriginal middle frame (time={middle_time:.0f} ps, idx={middle_idx}):")
    for ch, (mg, ng) in mid_details.items():
        print(f"  Chain {ch}: max_gap={mg:.1f} A, gaps>10A={ng}")
    print(f"  Overall: max_gap={mid_gap:.1f} A, total_gaps>10A={mid_ngaps}")

    # Check if middle frame is already good
    if mid_gap <= args.max_gap_threshold:
        print(f"\nMiddle frame backbone is already good (max_gap={mid_gap:.1f} <= {args.max_gap_threshold} A).")
        print("No replacement needed.")
        sys.exit(0)

    # Scan all cluster member frames
    print(f"\nScanning all {len(member_times)} cluster member frames...")
    results = []
    for t in member_times:
        idx = int(np.argmin(np.abs(timesteps - t)))
        max_gap, n_gaps, chain_details = calc_backbone_integrity(u, idx, check_chains)
        results.append({
            'time': t,
            'frame_idx': idx,
            'max_gap': max_gap,
            'n_gaps_gt10': n_gaps,
            'chain_details': chain_details,
        })

    # Sort by: (1) total gaps > 10 Å, (2) max gap
    results.sort(key=lambda x: (x['n_gaps_gt10'], x['max_gap']))

    # Report top candidates
    print(f"\nTop 10 frames with best backbone integrity:")
    print(f"{'Time (ps)':>12} {'Frame':>6} {'MaxGap':>8} {'Gaps>10':>8} ", end="")
    for ch in check_chains:
        print(f" Ch{ch}_max", end="")
    print()

    for r in results[:10]:
        print(f"{r['time']:12.0f} {r['frame_idx']:6d} {r['max_gap']:8.1f} {r['n_gaps_gt10']:8d} ", end="")
        for ch in check_chains:
            mg = r['chain_details'].get(ch, (0, 0))[0]
            print(f" {mg:7.1f}", end="")
        print()

    # Count good frames
    good = [r for r in results if r['max_gap'] <= args.max_gap_threshold]
    print(f"\nFrames with max_gap <= {args.max_gap_threshold} A: {len(good)}/{len(results)}")

    best = results[0]
    print(f"\nBest frame: time={best['time']:.0f} ps (idx={best['frame_idx']})")
    print(f"  max_gap={best['max_gap']:.1f} A, total_gaps>10A={best['n_gaps_gt10']}")
    for ch, (mg, ng) in best['chain_details'].items():
        print(f"  Chain {ch}: max_gap={mg:.1f} A, gaps>10A={ng}")

    if best['max_gap'] >= mid_gap:
        print(f"\nWARNING: Best frame ({best['max_gap']:.1f} A) is not better than middle ({mid_gap:.1f} A).")
        print("The backbone distortion may be present in all cluster frames.")

    # Extract best frame
    print(f"\nExtracting best frame...")
    extract_frames(u, best['frame_idx'], args.output_file, args.n_frames)
    print("Done.")


if __name__ == "__main__":
    main()
