"""
Workflow:

    -Get the number of frames of the trajectories in dir /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/*xtc. 
    the PDB files have the same basename.
    -The max number of frames is 32000, normalize the number of frames and write the occupancy in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/occupancy.csv

"""

import os
os.environ['OMP_NUM_THREADS'] = '1'

import glob
import MDAnalysis as mda
import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')



def get_trajectory_frames(trajectory_dir):
    """Get the number of frames for each trajectory file in the directory"""
    
    trajectory_files = glob.glob(os.path.join(trajectory_dir, "*.xtc"))
    frame_counts = {}
    
    print(f"Found {len(trajectory_files)} trajectory files:")
    
    for traj_file in sorted(trajectory_files):
        basename = os.path.basename(traj_file)
        name_without_ext = os.path.splitext(basename)[0]
        
        # Find corresponding PDB file
        pdb_file = os.path.join(trajectory_dir, f"{name_without_ext}.pdb")
        
        if not os.path.exists(pdb_file):
            print(f"Warning: PDB file not found for {basename}")
            continue
            
        try:
            # Load trajectory
            u = mda.Universe(pdb_file, traj_file)
            n_frames = len(u.trajectory)
            frame_counts[name_without_ext] = n_frames
            print(f"  {basename}: {n_frames} frames")
            
        except Exception as e:
            print(f"Error reading {basename}: {e}")
            continue
    
    return frame_counts

def calculate_normalized_occupancy(frame_counts, max_frames=32000):
    """Calculate normalized occupancy based on frame counts"""
    
    occupancy_data = []
    
    for trajectory_name, n_frames in frame_counts.items():
        # Normalize by max_frames
        normalized_occupancy = round(n_frames / max_frames, 2)
        occupancy_percent = round(normalized_occupancy * 100, 2)
        
        occupancy_data.append({
            'trajectory': trajectory_name,
            'frames': n_frames,
            'normalized_occupancy': normalized_occupancy,
            'occupancy_percent': occupancy_percent
        })
    
    return occupancy_data

def main():
    # Directory containing trajectories
    trajectory_dir=sys.argv[1]
    max_frames=int(sys.argv[2]) 

    # Output CSV file
    output_csv = os.path.join(trajectory_dir, "occupancy.csv")
    
    print("Analyzing trajectory files...")
    
    # Get frame counts
    frame_counts = get_trajectory_frames(trajectory_dir)
    
    if not frame_counts:
        print("No trajectory files found or processed successfully.")
        return
    
    # Calculate normalized occupancy
    occupancy_data = calculate_normalized_occupancy(frame_counts, max_frames=max_frames)
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(occupancy_data)
    
    # Sort by trajectory name for consistent output
    df = df.sort_values('trajectory')
    
    # Save to CSV
    df.to_csv(output_csv, index=False)
    
    print(f"\nOccupancy data saved to: {output_csv}")
    print("\nSummary:")
    print(df.to_string(index=False))
    
    # Print some statistics
    total_frames = df['frames'].sum()
    avg_occupancy = df['normalized_occupancy'].mean()
    max_frames_actual = df['frames'].max()
    
    print("\nStatistics:")
    print(f"  Total frames across all trajectories: {total_frames}")
    print(f"  Average normalized occupancy: {avg_occupancy:.2f}")
    print(f"  Maximum frames in any trajectory: {max_frames_actual}")
    print(f"  Number of trajectory groups: {len(df)}")

if __name__ == "__main__":
    main()