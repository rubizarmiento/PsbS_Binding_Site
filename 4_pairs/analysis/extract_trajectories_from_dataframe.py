"""
Reads DataFrame with the format: with the columns start_frame,end_frame, the rest of the columns are ignored,
and writes the trajectories between start_frame and end_frame in a directory.

Arguments:
    input (str): CSV file containing the DataFrame with columns start_frame and end_frame.
    sort: Sort the dataframe by the column "sort", e.g. lifetime_ns.
    n_trj: Save only the top N trajectories by the "sort" column, e.g. 3.
    f: Structure file compatible with MDAnalysis.
    trj: Trajectory file compatible with MDAnalysis.
    sel: Selection string for MDAnalysis, default "all"
    output_dir (str): Directory to write the extracted trajectories.
    preffix: Prefix for the output trajectory files.
    join: Write a single trajectory instead of multiple files.

Example:
python3 extract_trajectories_from_lifetime.py -i input.csv -f structure.pdb -trj trajectory.xtc -sort lifetime_ns -n_trj 3 -o output_dir
"""

import MDAnalysis as mda
import pandas as pd
import argparse
import os 

def parse_args():
    parser = argparse.ArgumentParser(description="Extract trajectories from a lifetime DataFrame.")
    parser.add_argument("-i", "--input_file", required=True, help="CSV file containing the DataFrame with columns start_frame and end_frame.")
    parser.add_argument("-f", "--structure", required=True, help="Structure file compatible with MDAnalysis.")
    parser.add_argument("-trj", "--trajectory", required=True, help="Trajectory file compatible with MDAnalysis.")
    parser.add_argument("-sel", "--selection", default="all", help="Selection string for MDAnalysis.")
    parser.add_argument("-sort", "--sort", required=True, help="Column to sort the DataFrame by.")
    parser.add_argument("-n_trj", "--n_trj", type=int, help="Save only the top N trajectories by the sort column.")
    parser.add_argument("-prefix", "--prefix", required=True, help="Prefix for the output trajectory files.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to write the extracted trajectories.")
    parser.add_argument("-join", "--join", action="store_true", help="Write a single trajectory instead of multiple files.")
    return parser.parse_args()  

def extract_trajectories(u, start_frame_arr, end_frame_arr, preffix, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    for i, (start, end) in enumerate(zip(start_frame_arr, end_frame_arr)):
        output_file = os.path.join(output_dir, f"{preffix}_{i+1}_{start}_{end}.xtc")
        
        # Create a writer for this trajectory segment
        with mda.Writer(output_file, n_atoms=u.atoms.n_atoms) as writer:
            # Iterate through the specific frame range
            for frame_idx in range(start, end + 1):  # +1 to include end frame
                try:
                    u.trajectory[frame_idx]  # Go to specific frame
                    writer.write(u.atoms)   # Write current frame
                except IndexError:
                    print(f"Warning: Frame {frame_idx} not found in trajectory")
                    break
        
        print(f"Extracted trajectory {i+1}/{len(start_frame_arr)}: frames {start}-{end} -> {output_file}")

def join_trajectories(u, start_frame_arr, end_frame_arr, prefix, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, f"{prefix}_joined.xtc")
    total_frames = 0
    
    # Create a single writer for all trajectory segments
    with mda.Writer(output_file, n_atoms=u.atoms.n_atoms) as writer:
        for i, (start, end) in enumerate(zip(start_frame_arr, end_frame_arr)):
            frames_written = 0
            
            # Iterate through the specific frame range
            for frame_idx in range(start, end + 1):  # +1 to include end frame
                try:
                    u.trajectory[frame_idx]  # Go to specific frame
                    writer.write(u.atoms)   # Write current frame
                    frames_written += 1
                    total_frames += 1
                except IndexError:
                    print(f"Warning: Frame {frame_idx} not found in trajectory")
                    break
            
            print(f"Segment {i+1}/{len(start_frame_arr)}: frames {start}-{end} ({frames_written} frames written)")
    
    print(f"All segments joined into single trajectory: {output_file}")
    print(f"Total frames written: {total_frames}")

def main():
    args = parse_args()
    df = pd.read_csv(args.input_file)

    if args.n_trj:
        df = df.sort_values(by=args.sort, ascending=False)
        if args.n_trj:
            df = df.head(args.n_trj)
            df.to_csv(os.path.join(args.output_dir, f"{args.prefix}_top_{args.n_trj}.csv"), index=False)
        else:
            df.to_csv(os.path.join(args.output_dir, f"{args.prefix}_sorted.csv"), index=False)

    start_frame_arr = df['start_frame'].values
    end_frame_arr = df['end_frame'].values
    u = mda.Universe(args.structure, args.trajectory)
    u = u.select_atoms(args.selection)
    print(f"Universe has {u.atoms.n_residues} residues.")
    if args.join:
        join_trajectories(u, start_frame_arr, end_frame_arr, args.prefix, args.output_dir)
    else:
        extract_trajectories(u, start_frame_arr, end_frame_arr,args.prefix, args.output_dir)

if __name__ == "__main__":
    main()
