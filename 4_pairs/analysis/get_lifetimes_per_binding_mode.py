"""
Arguments:
    -d: Directories containing the trajectory files (e.g., /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/)
    -dt: Time step between frames in nanoseconds (e.g., 1)
    -o: Output CSV file path (e.g., /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes.csv)

Workflow:

    -Get the number of frames of the trajectories in dirs /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/*xtc. 
    The script counts all the files with the extension .xtc in the specified directories, and for each file, it looks for a corresponding PDB file with the same basename. 
    If the PDB file is found, it loads the trajectory using MDAnalysis and counts the number of frames. 
    the PDB files have the same basename.
    -The max number of frames is 32000, normalize the number of frames and write the lifetimes in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes.csv

"""

import os
os.environ['OMP_NUM_THREADS'] = '1'

import glob
import MDAnalysis as mda
import pandas as pd
import sys
import warnings
import argparse
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')

def parser_args():
    parser = argparse.ArgumentParser(description="Calculate lifetimes of trajectories based on frame counts.")
    parser.add_argument("-d", "--dir", type=str, nargs='+', required=True, help="Directories containing the trajectory files (e.g., /path/to/trajectories/)")
    parser.add_argument("-dt", "--dt", type=float, default=1.0, help="Time step between frames in nanoseconds (default: 1.0)")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file path (e.g., /path/to/output/lifetimes.csv)")
    return parser.parse_args()

def get_trajectory_frames(trajectory_dirs):
    """Get the number of frames for each trajectory file in the directories"""
    
    frame_counts = {}
    
    for trajectory_dir in trajectory_dirs:
        trajectory_files = glob.glob(os.path.join(trajectory_dir, "*.xtc"))
        print(f"Found {len(trajectory_files)} trajectory files in {trajectory_dir}:")
        
        for traj_file in sorted(trajectory_files):
            basename = os.path.basename(traj_file)
            name_without_ext = os.path.splitext(basename)[0]
            
            # Find corresponding PDB file
            pdb_file = os.path.join(trajectory_dir, f"{name_without_ext}.pdb")
            
            if not os.path.exists(pdb_file):
                print(f"Warning: PDB file not found for {basename} in {trajectory_dir}")
                continue
                
            try:
                # Load trajectory
                u = mda.Universe(pdb_file, traj_file)
                n_frames = len(u.trajectory)
                frame_counts[name_without_ext] = n_frames
                print(f"  {basename}: {n_frames} frames")
                
            except Exception as e:
                print(f"Error reading {basename} in {trajectory_dir}: {e}")
                continue
    
    return frame_counts

def main():
    args = parser_args()
    trajectory_dirs = args.dir
    dt = args.dt
    output_csv = args.output
    
    print("Analyzing trajectory files...")
    
    # Get frame counts
    frame_counts = get_trajectory_frames(trajectory_dirs)

    lifetimes = {traj: frames * dt for traj, frames in frame_counts.items()}    
    
    if not frame_counts:
        print("No trajectory files found or processed successfully.")
        return
    
    # Datagrame with basenames, frames, lifetimes
    df = pd.DataFrame({
        'trajectory': list(frame_counts.keys()),
        'frames': list(frame_counts.values()),
        'lifetimes': [lifetimes[traj] for traj in frame_counts.keys()]
    })
    
    # Sort by trajectory name for consistent output
    df = df.sort_values('trajectory')
    
    # Save to CSV
    df.to_csv(output_csv, index=False)
    
    print(f"\Lifetimes data saved to: {output_csv}")
    print("\nSummary:")
    print(df.to_string(index=False))
    
    # Print some statistics
    total_frames = df['frames'].sum()
    avg_lifetime = df['lifetimes'].mean()
    max_frames_actual = df['frames'].max()
    
    print("\nStatistics:")
    print(f"  Total frames across all trajectories: {total_frames}")
    print(f"  Average lifetime: {avg_lifetime:.2f}")
    print(f"  Maximum frames in any trajectory: {max_frames_actual}")
    print(f"  Number of trajectory groups: {len(df)}")

if __name__ == "__main__":
    main()