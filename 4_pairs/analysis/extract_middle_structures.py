#!/usr/bin/env python3
"""
Script to extract middle structures from clustering results.

This script parses cluster.log files from GROMACS clustering analysis,
extracts the middle structure frame numbers     parser.add_argument('--clustering-dir', required=True,
                       help='Directory containing clustering results')
    parser.add_argument('--sim-base-dir', required=True,
                       help='Base directory containing trajectory files')
    parser.add_argument('--output-dir', default='extracted_structures',
                       help='Output directory for extracted structures')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without actually extracting frames')
    parser.add_argument('--biggest-only', action='store_true',
                       help='Only extract the biggest cluster from clust_c075 (default: extract all clusters)')
    parser.add_argument('--clust-c075-only', action='store_true',
                       help='Only process clust_c075 directories (default: process all clust_* directories)') cluster,
and uses gmx trjconv to extract those frames from the corresponding trajectory files.
"""

import os
import re
import subprocess
import argparse
from pathlib import Path
from typing import Dict, List, Tuple


def parse_cluster_log(log_file: str) -> Dict[int, Dict[str, any]]:
    """
    Parse a cluster.log file to extract cluster information.

    Args:
        log_file: Path to the cluster.log file

    Returns:
        Dictionary with cluster ID as key and cluster info as value
    """
    clusters = {}

    try:
        with open(log_file, 'r') as f:
            content = f.read()

        # Find the cluster table section
        # Look for lines that match the pattern: "cl. | #st  rmsd | middle rmsd | cluster members"
        lines = content.split('\n')

        # Find the header line
        header_idx = -1
        for i, line in enumerate(lines):
            if 'cl.' in line and '#st' in line and 'rmsd' in line:
                header_idx = i
                break

        if header_idx == -1:
            print(f"Warning: Could not find cluster table header in {log_file}")
            return clusters

        # Parse cluster lines
        for line in lines[header_idx + 1:]:
            line = line.strip()
            if not line or line.startswith('Found') or line.startswith('Writing'):
                continue

            # Match pattern like: "  1 | 1784  0.422 | 944000 .351 | 121000 123000 ..."
            # Also handle scientific notation like "1.35e+06"
            match = re.match(r'\s*(\d+)\s*\|\s*(\d+)\s+([\d.]+)\s*\|\s*([\d.eE+-]+)\s+[\d.]+\s*\|\s*(.+)', line)
            if match:
                cluster_id = int(match.group(1))
                n_structures = int(match.group(2))
                rmsd = float(match.group(3))
                middle_frame_str = match.group(4)
                
                # Handle scientific notation
                try:
                    middle_frame = int(float(middle_frame_str))
                except ValueError:
                    print(f"Warning: Could not parse middle frame '{middle_frame_str}' in line: {line}")
                    continue
                
                # members = match.group(5).strip()

                clusters[cluster_id] = {
                    'n_structures': n_structures,
                    'rmsd': rmsd,
                    'middle_frame': middle_frame,
                    'log_file': log_file
                }

    except Exception as e:
        print(f"Error parsing {log_file}: {e}")

    return clusters


def find_trajectory_file(cluster_dir: str, base_dir: str) -> Tuple[str, str]:
    """
    Map cluster directory name to corresponding trajectory file.

    Args:
        cluster_dir: Name of the cluster directory (e.g., '1_n_s' or '1_sim_1_A2_s')
        base_dir: Base directory containing trajectory files

    Returns:
        Tuple of (trajectory_path, tpr_path) or empty strings if not found
    """
    # Try new format first: direct mapping like '1_n_s' -> 'base_dir/1_n_s.xtc'
    traj_file = f"{cluster_dir}.xtc"
    tpr_file = f"{cluster_dir}.tpr"
    pdb_file = f"{cluster_dir}.pdb"
    
    traj_path = os.path.join(base_dir, traj_file)
    tpr_path = os.path.join(base_dir, tpr_file)
    pdb_path = os.path.join(base_dir, pdb_file)
    
    if os.path.exists(traj_path):
        # Try to find structure file (tpr or pdb)
        if os.path.exists(tpr_path):
            return traj_path, tpr_path
        elif os.path.exists(pdb_path):
            return traj_path, pdb_path
        else:
            print(f"Warning: Structure file not found for {cluster_dir} (looked for .tpr and .pdb)")
            return traj_path, ""
    
    print(f"Warning: Trajectory file not found: {traj_path}")

    # Fall back to old format: {number}_sim_{sim_num}_{chains}
    match = re.match(r'(\d+)_sim_(\d+)_(.+)', cluster_dir)
    if not match:
        print(f"Warning: Could not parse cluster directory name: {cluster_dir}")
        return "", ""

    cluster_num = match.group(1)
    sim_num = int(match.group(2))
    chains = match.group(3)

    # Look for trajectory file in the corresponding sim directory
    sim_dir = os.path.join(base_dir, f'sim_{sim_num}')
    if not os.path.exists(sim_dir):
        print(f"Warning: Simulation directory not found: {sim_dir}")
        return "", ""

    # Look for trajectory files - prefer 3-run.xtc if it exists
    trajectory_files = ['3-run.xtc', 'analysis.xtc', 'proteins.xtc', 'aligned.xtc']
    tpr_files = ['3-run.tpr', 'topol.tpr', 'analysis.tpr', '2-eq.tpr']

    traj_path = ""
    for traj_file in trajectory_files:
        traj_path_candidate = os.path.join(sim_dir, traj_file)
        if os.path.exists(traj_path_candidate):
            traj_path = traj_path_candidate
            break
    
    if not traj_path:
        print(f"Warning: No trajectory file found in {sim_dir}")
        return "", ""
    
    # Find corresponding TPR file
    tpr_path = ""
    for tpr_file in tpr_files:
        tpr_path_candidate = os.path.join(sim_dir, tpr_file)
        if os.path.exists(tpr_path_candidate):
            tpr_path = tpr_path_candidate
            break
    
    return traj_path, tpr_path


def extract_frame(trajectory_file: str, frame_number: int, output_file: str, structure_file: str = None) -> bool:
    """
    Extract a specific frame from a trajectory file using gmx trjconv.

    Args:
        trajectory_file: Path to the trajectory file
        frame_number: Frame number to extract (0-based)
        output_file: Path for the output PDB file
        structure_file: Path to the structure file (.pdb or .tpr, optional)

    Returns:
        True if successful, False otherwise
    """
    try:
        # Find structure file if not provided - prefer PDB to preserve chain IDs
        if not structure_file:
            # First try: PDB file with same base name as trajectory file
            traj_basename = os.path.basename(trajectory_file)
            pdb_basename = traj_basename.replace('.xtc', '.pdb')
            pdb_path = os.path.join(os.path.dirname(trajectory_file), pdb_basename)
            
            if os.path.exists(pdb_path):
                structure_file = pdb_path
            else:
                # Try TPR file as fallback
                tpr_basename = traj_basename.replace('.xtc', '.tpr')
                tpr_path = os.path.join(os.path.dirname(trajectory_file), tpr_basename)
                if os.path.exists(tpr_path):
                    structure_file = tpr_path
                else:
                    # Fall back to old method: look for standard files in sim directory
                    sim_dir = os.path.dirname(trajectory_file)
                    structure_files = ['3-run.pdb', 'analysis.pdb', '3-run.tpr', 'topol.tpr', 'analysis.tpr', '2-eq.tpr']

                    for struct_file in structure_files:
                        struct_path = os.path.join(sim_dir, struct_file)
                        if os.path.exists(struct_path):
                            structure_file = struct_path
                            break

        if not structure_file or not os.path.exists(structure_file):
            print(f"Warning: Structure file not found for {trajectory_file}")
            print(f"  Looked in: {os.path.dirname(trajectory_file)}")
            return False

        # Run gmx trjconv to extract the frame
        cmd = [
            'gmx', 'trjconv',
            '-f', trajectory_file,
            '-s', structure_file,
            '-o', output_file
        ]

        # Use -b and -e to select specific frame
        # Frame numbers in cluster.log correspond to time in ps
        # (since nstxout-compressed = 1000 and dt = 0.001 ps/step)
        time_ps = frame_number

        cmd.extend(['-b', str(time_ps), '-e', str(time_ps)])

        print(f"Running: {' '.join(cmd)}")
        print(f"Input: 0")

        # Run the command with input
        result = subprocess.run(cmd, input='0\n', text=True, capture_output=True)

        if result.returncode == 0:
            print(f"Successfully extracted frame {frame_number} to {output_file}")
            return True
        else:
            print(f"Error extracting frame {frame_number}: {result.stderr}")
            print(f"Command output: {result.stdout}")
            return False

    except Exception as e:
        print(f"Error extracting frame {frame_number}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description='Extract middle structures from clustering results')
    parser.add_argument('--clustering-dir', required=True,
                       help='Directory containing clustering results')
    parser.add_argument('--sim-base-dir', required=True,
                       help='Base directory containing trajectory files')
    parser.add_argument('--output-dir', default='extracted_structures',
                       help='Output directory for extracted structures')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without actually extracting frames')
    parser.add_argument('--biggest-only', action='store_true',
                       help='Only extract the biggest cluster from clust_c075 (default: extract all clusters)')
    parser.add_argument('--clust-c075-only', action='store_true',
                       help='Only process clust_c075 directories (default: process all clust_* directories)')

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Find all cluster directories (main directories containing clust_c* subdirs)
    cluster_dirs = []
    for item in os.listdir(args.clustering_dir):
        item_path = os.path.join(args.clustering_dir, item)
        if os.path.isdir(item_path):
            # Check if this directory contains clust_* subdirectories
            has_clust_subdirs = any(subdir.startswith('clust_') for subdir in os.listdir(item_path)
                                   if os.path.isdir(os.path.join(item_path, subdir)))
            if has_clust_subdirs:
                cluster_dirs.append(item_path)

    print(f"Found {len(cluster_dirs)} cluster directories")

    # Process each cluster directory
    all_clusters = {}
    for cluster_dir_path in cluster_dirs:
        cluster_dir_name = os.path.basename(cluster_dir_path)
        print(f"\nProcessing cluster directory: {cluster_dir_name}")

        # Find trajectory file for this cluster
        trajectory_file, structure_file = find_trajectory_file(cluster_dir_name, args.sim_base_dir)

        if not trajectory_file:
            print(f"No trajectory file found for {cluster_dir_name}")
            continue

        print(f"Using trajectory: {trajectory_file}")
        print(f"Using structure: {structure_file}")

        # Process clust_* subdirectories based on flags
        if args.clust_c075_only:
            # Only process clust_c075
            clust_c075_path = os.path.join(cluster_dir_path, 'clust_c075')
            if os.path.isdir(clust_c075_path):
                cluster_log = os.path.join(clust_c075_path, 'cluster.log')
                if os.path.exists(cluster_log):
                    print(f"  Processing clust_c075/cluster.log")
                    clusters = parse_cluster_log(cluster_log)

                    if clusters:
                        print(f"  Found {len(clusters)} clusters in clust_c075")

                        if args.biggest_only:
                            # Find the biggest cluster (most structures)
                            biggest_cluster_id = None
                            max_structures = 0

                            for cluster_id, cluster_info in clusters.items():
                                if cluster_info['n_structures'] > max_structures:
                                    max_structures = cluster_info['n_structures']
                                    biggest_cluster_id = cluster_id

                            if biggest_cluster_id is not None:
                                print(f"  Biggest cluster: {biggest_cluster_id} with {max_structures} structures")
                                clusters_to_process = {biggest_cluster_id: clusters[biggest_cluster_id]}
                            else:
                                print(f"  No clusters found in clust_c075/cluster.log")
                                clusters_to_process = {}
                        else:
                            clusters_to_process = clusters

                        # Process the selected clusters
                        for cluster_id, cluster_info in clusters_to_process.items():
                            frame_number = cluster_info['middle_frame']
                            suffix = "_biggest" if args.biggest_only else f"_cluster_{cluster_id}"
                            output_file = os.path.join(
                                args.output_dir,
                                f"{cluster_dir_name}_clust_c075{suffix}_frame_{frame_number}.pdb"
                            )

                            if args.dry_run:
                                print(f"    Would extract frame {frame_number} from {trajectory_file} to {output_file}")
                            else:
                                success = extract_frame(trajectory_file, frame_number, output_file, structure_file)
                                if success:
                                    cluster_info['extracted_file'] = output_file

                            all_clusters[f"{cluster_dir_name}_clust_c075{suffix}"] = cluster_info
                    else:
                        print(f"  No clusters found in clust_c075/cluster.log")
                else:
                    print(f"  cluster.log not found in clust_c075")
            else:
                print(f"  clust_c075 directory not found in {cluster_dir_name}")
        else:
            # Original behavior: process all clust_* subdirectories
            for subdir in os.listdir(cluster_dir_path):
                subdir_path = os.path.join(cluster_dir_path, subdir)
                if os.path.isdir(subdir_path) and subdir.startswith('clust_'):
                    cluster_log = os.path.join(subdir_path, 'cluster.log')
                    if os.path.exists(cluster_log):
                        print(f"  Processing {subdir}/cluster.log")
                        clusters = parse_cluster_log(cluster_log)

                        if clusters:
                            print(f"  Found {len(clusters)} clusters in {subdir}")

                            for cluster_id, cluster_info in clusters.items():
                                frame_number = cluster_info['middle_frame']
                                output_file = os.path.join(
                                    args.output_dir,
                                    f"{cluster_dir_name}_{subdir}_cluster_{cluster_id}_frame_{frame_number}.pdb"
                                )

                                if args.dry_run:
                                    print(f"    Would extract frame {frame_number} from {trajectory_file} to {output_file}")
                                else:
                                    success = extract_frame(trajectory_file, frame_number, output_file, structure_file)
                                    if success:
                                        cluster_info['extracted_file'] = output_file

                                all_clusters[f"{cluster_dir_name}_{subdir}_cluster_{cluster_id}"] = cluster_info
                        else:
                            print(f"  No clusters found in {subdir}/cluster.log")    # Print summary
    print("\n=== Summary ===")
    print(f"Processed {len(cluster_dirs)} cluster directories")
    print(f"Found {len(all_clusters)} clusters")

    if not args.dry_run:
        extracted_count = sum(1 for c in all_clusters.values() if 'extracted_file' in c)
        print(f"Successfully extracted {extracted_count} structures")

        # Save summary to file
        summary_file = os.path.join(args.output_dir, 'extraction_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("Cluster Extraction Summary\n")
            f.write("=" * 50 + "\n\n")

            for cluster_name, info in all_clusters.items():
                f.write(f"Cluster: {cluster_name}\n")
                f.write(f"  Structures: {info['n_structures']}\n")
                f.write(f"  RMSD: {info['rmsd']:.3f}\n")
                f.write(f"  Middle frame: {info['middle_frame']}\n")
                if 'extracted_file' in info:
                    f.write(f"  Extracted file: {info['extracted_file']}\n")
                f.write("\n")

        print(f"Summary saved to {summary_file}")


if __name__ == '__main__':
    main()