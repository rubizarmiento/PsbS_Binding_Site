"""
Workflow:
    -Read /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv with the content:
    tag original new count tag_number old_chains
    n_s sim_4_A4_N_S sim_4_A4_n_s 6 1 N_S
    -Join all the trajectories with the same tag_number in the files /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/{original}_aligned.xtc
    and save them as /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/{tag_number}_{tag}.xtc
    -Also copy the corresponding .pdb, .top, and .tpr files to the same output directory
"""

import os
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
import MDAnalysis as mda
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')

def concatenate_trajectories(trajectory_files, output_file, output_dir):
    """Concatenate multiple trajectories into a single output file and copy topology/PDB/TPR files"""
    
    if not trajectory_files:
        print(f"No trajectory files found for concatenation")
        return False
    
    print(f"Concatenating {len(trajectory_files)} trajectories:")
    for traj_file in trajectory_files:
        print(f"  - {traj_file}")
    
    # Load first trajectory to get structure info
    first_traj = trajectory_files[0]
    if not os.path.exists(first_traj):
        print(f"First trajectory file not found: {first_traj}")
        return False
    
    # Try to find corresponding PDB file
    pdb_file = first_traj.replace('_aligned.xtc', '_aligned.pdb')
    if not os.path.exists(pdb_file):
        # Try grouped PDB file
        pdb_file = first_traj.replace('_aligned.xtc', '_grouped.pdb')
        if not os.path.exists(pdb_file):
            # Try regular PDB file
            pdb_file = first_traj.replace('_aligned.xtc', '.pdb')
            if not os.path.exists(pdb_file):
                # Try without _aligned suffix
                pdb_file = first_traj.replace('.xtc', '.pdb')
                if not os.path.exists(pdb_file):
                    print(f"No corresponding PDB file found for {first_traj}")
                    print(f"Tried: {first_traj.replace('_aligned.xtc', '_aligned.pdb')}")
                    print(f"Tried: {first_traj.replace('_aligned.xtc', '_grouped.pdb')}")
                    print(f"Tried: {first_traj.replace('_aligned.xtc', '.pdb')}")
                    print(f"Tried: {first_traj.replace('.xtc', '.pdb')}")
                    return False
    
    # Try to find corresponding TOP file
    top_file = first_traj.replace('_aligned.xtc', '.top')
    if not os.path.exists(top_file):
        # Try without _aligned suffix
        top_file = first_traj.replace('.xtc', '.top')
        if not os.path.exists(top_file):
            print(f"No corresponding TOP file found for {first_traj}")
            print(f"Tried: {first_traj.replace('_aligned.xtc', '.top')}")
            print(f"Tried: {first_traj.replace('.xtc', '.top')}")
            top_file = None
    
    # Try to find corresponding TPR file
    tpr_file = first_traj.replace('_aligned.xtc', '.tpr')
    if not os.path.exists(tpr_file):
        # Try without _aligned suffix
        tpr_file = first_traj.replace('.xtc', '.tpr')
        if not os.path.exists(tpr_file):
            print(f"No corresponding TPR file found for {first_traj}")
            print(f"Tried: {first_traj.replace('_aligned.xtc', '.tpr')}")
            print(f"Tried: {first_traj.replace('.xtc', '.tpr')}")
            tpr_file = None
    
    print(f"Using structure file: {pdb_file}")
    
    # Copy PDB file to output directory
    output_pdb = str(output_file).replace('.xtc', '.pdb')
    import shutil
    shutil.copy2(pdb_file, output_pdb)
    print(f"Copied PDB file to: {output_pdb}")
    
    if top_file:
        print(f"Using topology file: {top_file}")
        # Copy topology file to output directory
        output_top = str(output_file).replace('.xtc', '.top')
        shutil.copy2(top_file, output_top)
        print(f"Copied topology file to: {output_top}")
    
    if tpr_file:
        print(f"Using TPR file: {tpr_file}")
        # Copy TPR file to output directory
        output_tpr = str(output_file).replace('.xtc', '.tpr')
        shutil.copy2(tpr_file, output_tpr)
        print(f"Copied TPR file to: {output_tpr}")
    
    # Load first universe
    u_first = mda.Universe(pdb_file, first_traj)
    n_atoms = len(u_first.atoms)
    
    print(f"First trajectory has {len(u_first.trajectory)} frames and {n_atoms} atoms")
    
    # Create output writer
    with mda.Writer(str(output_file), n_atoms) as writer:
        # Write all frames from first trajectory
        for ts in u_first.trajectory:
            writer.write(u_first.atoms)
        
        # Write frames from remaining trajectories
        for traj_file in trajectory_files[1:]:
            if not os.path.exists(traj_file):
                print(f"Warning: Trajectory file not found: {traj_file}")
                continue
                
            print(f"Processing {traj_file}...")
            
            # Load trajectory
            u = mda.Universe(pdb_file, traj_file)
            
            # Check if atom count matches
            if len(u.atoms) != n_atoms:
                print(f"Warning: Atom count mismatch in {traj_file} ({len(u.atoms)} vs {n_atoms})")
                continue
            
            print(f"  Adding {len(u.trajectory)} frames")
            
            # Write all frames
            for ts in u.trajectory:
                writer.write(u.atoms)
    
    print(f"Concatenated trajectory saved as: {output_file}")
    return True

def main():
    # File paths
    csv_file = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv"
    trj_dir = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj"
    output_dir = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read CSV file
    print(f"Reading CSV file: {csv_file}")
    df = pd.read_csv(csv_file, sep='\s+')
    
    print(f"Found {len(df)} entries in CSV")
    print("Grouping by tag_number...")
    
    # Group by tag_number
    grouped = df.groupby('tag_number')
    
    for tag_number, group in grouped:
        print(f"\nProcessing tag_number {tag_number}:")
        
        # Get the tag for this group (should be the same for all rows in group)
        tag = group['tag'].iloc[0]
        print(f"  Tag: {tag}")
        
        # Collect all trajectory files for this tag_number
        trajectory_files = []
        
        for idx, row in group.iterrows():
            original = row['original']
            traj_file = os.path.join(trj_dir, f"{original}_aligned.xtc")
            
            if os.path.exists(traj_file):
                trajectory_files.append(traj_file)
                print(f"  Found trajectory: {original}_aligned.xtc")
            else:
                print(f"  Warning: Trajectory not found: {original}_aligned.xtc")
        
        if not trajectory_files:
            print(f"  No trajectory files found for tag_number {tag_number}")
            continue
        
        # Create output filename
        output_file = os.path.join(output_dir, f"{tag_number}_{tag}.xtc")
        
        print(f"  Concatenating {len(trajectory_files)} trajectories...")
        
        # Concatenate trajectories
        success = concatenate_trajectories(trajectory_files, output_file, output_dir)
        
        if success:
            print(f"  Successfully created: {tag_number}_{tag}.xtc")
        else:
            print(f"  Failed to create: {tag_number}_{tag}.xtc")
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()