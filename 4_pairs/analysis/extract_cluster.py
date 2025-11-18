"""


python3 ${script}/extract_cluster.py -f structure.pdb -tr trajectory.xtc -g cluster.log -o output.pdb

Arguments:
    -f, --structure-file: Path to the structure file (.pdb, .gro, etc.)
    -trj, --trajectory-file: Path to the trajectory file (.xtc, .dcd, etc.)
    -g, --cluster-log: Path to the cluster log file
    -o, --output-file: Path for the output PDB file

"""


import MDAnalysis as mda
# Ignore warnings from MDAnalysis
import warnings
warnings.filterwarnings("ignore")
import argparse
import numpy as np  

def parser():

    parser = argparse.ArgumentParser(description="Extract middle structures from clusters")
    parser.add_argument('-f', '--structure-file', type=str, required=True, help='Path to the structure file (.pdb, .gro, etc.)')
    parser.add_argument('-trj', '--trajectory-file', type=str, required=True, help='Path to the trajectory file (.xtc, .dcd, etc.)')
    parser.add_argument('-g', '--cluster-log', type=str, required=True, help='Path to the cluster log file')
    parser.add_argument('-o', '--output-file', type=str, required=True, help='Path for the output PDB file')
    parser.add_argument('-n', '--n-frames', type=int, required=False, default=1, help='Number of frames to extract before the middle cluster time (default: 1, extracts only the middle frame)')
    parser.add_argument('--first-member', action='store_true', help='Extract the first cluster member instead of the middle structure')
    #parser.add_argument('-id', '--id-cluster', type=int, required=False, help='Cluster ID to extract (optional, if not provided, extracts all clusters). For the biggest cluster, use 1, as it is the first one in the log file.)')

    return parser.parse_args()




def extract_frame(trajectory_file: str, structure_file: str, timestamp: float, output_file: str, n_frames: int = 1) -> bool:
    import os
    u = mda.Universe(structure_file, trajectory_file) 

    total_frames = len(u.trajectory)
    trj_length = u.trajectory.totaltime
    trj_dt = u.trajectory.dt
    t0 = u.trajectory[0].time
    tf = u.trajectory[-1].time

    timestep_array = np.arange(t0, tf + trj_dt, trj_dt, dtype=float)

    # Find the frame closest to the target timestamp
    closest_value = timestep_array[np.argmin(np.abs(timestep_array - timestamp))]
    middle_frame = np.where(timestep_array == closest_value)[0][0]
    print(f"Middle frame: {middle_frame} at time {closest_value} ps")

    # Calculate the range of frames to extract
    start_frame = max(0, middle_frame - n_frames + 1)
    end_frame = middle_frame + 1

    print(f"Extracting {n_frames} frame(s) from frame {start_frame} to {middle_frame}")

    # If only extracting 1 frame, use the original behavior
    if n_frames == 1:
        u.trajectory[middle_frame]
        u.atoms.write(output_file)
        print(f"Wrote to {output_file}")
    else:
        # Save middle frame (without enumeration)
        base_output = os.path.splitext(output_file)[0]
        u.trajectory[middle_frame]
        u.atoms.write(output_file)
        print(f"Middle frame (frame {middle_frame}) -> {output_file}")

        # Save previous frames with enumeration
        frame_counter = 1
        for frame_idx in range(start_frame, middle_frame):
            u.trajectory[frame_idx]
            frame_time = u.trajectory.time
            output_path = f"{base_output}_prev_{frame_counter}.pdb"
            u.atoms.write(output_path)
            print(f"Previous frame {frame_counter}: frame index {frame_idx} (time {frame_time} ps) -> {output_path}")
            frame_counter += 1

def log_cluster_to_clusterid_and_timestamps(log_file: str, first_member: bool = False) -> float:
    # Get the line
    # cl. | #st  rmsd | middle rmsd | cluster members
    # And the next one, return the 4th numerical value (middle rmsd) or 5th (first member)
    #1 | 3949  0.348 | 3.521e+06 .263 |   2000   3000   4000   5000   6000   7000   8000
    with open(log_file, 'r') as f:
        lines = f.readlines()
    n_lines = len(lines)
    next_line = None
    for i in range(n_lines):
        line = lines[i].strip()
        if line.startswith("cl. | #st  rmsd | middle rmsd | cluster members"):
            next_line = lines[i+1].strip()
            break
    else:
        raise ValueError("Cluster log format not recognized")
    # Split by spaces
    parts = next_line.split()

    # Remove empty parts and |
    parts = [p for p in parts if p and p != '|']

    if first_member:
        # The first member is the 5th value (after cl, #st, rmsd, middle rmsd)
        # Index 4 is the first cluster member timestamp
        timestamp = parts[4]
    else:
        # Keep only the fourth value (middle rmsd)
        timestamp = parts[3]

    timestamp = float(timestamp)

    return timestamp


def main():
    args = parser()
    timestamp = log_cluster_to_clusterid_and_timestamps(args.cluster_log, first_member=args.first_member)

    extract_frame(args.trajectory_file, args.structure_file, timestamp, args.output_file, args.n_frames)

if __name__ == "__main__":
     main()

