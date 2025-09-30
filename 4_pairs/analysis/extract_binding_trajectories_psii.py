"""
Reads DataFrame with the format: resid_i,resid_j,start_frame,end_frame,frames,lifetime_ns
resid_1 are segids and resid_2 are chainIDs
Extracts sub    # If topology file is provided, modify it and run grompp
    if top_file and mdp_file and chain_id and seg_id:
        print(f"Topology file provided: {top_file}")
        print(f"MDP file provided: {mdp_file}")
        print(f"Chain ID: {chain_id}, Seg ID: {seg_id}")
        
        output_top = output_pdb.with_suffix('.top')
        output_tpr = output_pdb.with_suffix('.tpr')
        
        try:
            modify_topology_file(top_file, chain_id, seg_id, output_top)
            success = run_grompp(mdp_file, str(output_pdb), str(output_top), str(output_tpr))
            if not success:
                print(f"Warning: grompp failed for {output_pdb} - topology/coordinates mismatch is expected for subset extractions")
                print(f"Trajectory files are still valid: {output_xtc}")
        except Exception as e:
            print(f"Warning: Error processing topology/TPR for {output_pdb}: {e}")
            print(f"Trajectory files are still valid: {output_xtc}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Skipping topology processing - missing arguments: top_file={top_file}, mdp_file={mdp_file}, chain_id={chain_id}, seg_id={seg_id}")aj between start_frame and end_frame for each pair. Also extracts the first frame as a PDB file.

Arguments:
-f: Structure file compatible with MDAnalysis.
-trj: Trajectory file compatible with MDAnalysis.
-df: DataFrame file with the format: resid_i,resid_j,start_frame,end_frame,frames,lifetime_ns
-odir: Directory to write the extracted trajectories.
-preffix: Prefix for the output trajectory files.
-p: Topology file (.top) to modify for each trajectory.
-mdp: MDP file for grompp.

python extract_binding_trajectories_psii.py \
  -f structure.pdb \
  -trj trajectory.xtc \
  -df binding_events.csv \
  -odir output_dir \
  -p original.top \
  -mdp /path/to/em.mdp


"""

import pandas as pd
import MDAnalysis as mda
import argparse
import os
import subprocess
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')
def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract binding trajectories from DataFrame")
    parser.add_argument("-f", "--structure", required=True, help="Structure file compatible with MDAnalysis")
    parser.add_argument("-trj", "--trajectory", required=True, help="Trajectory file compatible with MDAnalysis")
    parser.add_argument("-df", "--dataframe", required=True, help="DataFrame file with binding events")
    parser.add_argument("-odir", "--output_dir", required=True, help="Directory to write extracted trajectories")
    parser.add_argument("-preffix", "--prefix", default="binding", help="Prefix for output trajectory files")
    parser.add_argument("-p", "--topology", help="Topology file (.top) to modify for each trajectory")
    parser.add_argument("-mdp", "--mdp_file", help="MDP file for grompp")
    
    return parser.parse_args()

def modify_topology_file(top_file, chain_ids, seg_id, output_top):
    """Modify topology file to comment out molecules except those containing chainIDs and segid"""
    
    print(f"Processing topology file: {top_file}")
    print(f"Looking for chain_ids: {chain_ids}, seg_id: {seg_id}")
    
    try:
        with open(str(top_file), 'r') as f:
            lines = f.readlines()
        
        print(f"Read {len(lines)} lines from topology file")
        
        modified_lines = []
        in_molecules_section = False
        
        for line in lines:
            stripped = line.strip()
            
            # Check if we're entering the molecules section
            if stripped == "[ molecules ]":
                in_molecules_section = True
                modified_lines.append(line)
                continue
            
            # If we're in molecules section and line is not empty
            if in_molecules_section and stripped and not stripped.startswith(';'):
                # Check if line contains any of the chain_ids or seg_id
                seg_id_str = str(seg_id)
                if any(f"chain_{cid}" in stripped or f"psbs_chain_{cid}" in stripped for cid in chain_ids) or seg_id_str in stripped:
                    modified_lines.append(line)
                    print(f"Keeping molecule line: {stripped}")
                else:
                    # Comment out the line
                    modified_lines.append(f";{line}")
            elif in_molecules_section and stripped.startswith(';') and any(f"; chainID {cid}" in stripped for cid in chain_ids):
                # Uncomment cofactor lines with matching chainID comment
                uncommented_line = line[1:]  # Remove the leading ;
                modified_lines.append(uncommented_line)
                print(f"Uncommenting cofactor line: {stripped[1:]}")
            else:
                modified_lines.append(line)
        
        # Write modified topology
        with open(str(output_top), 'w') as f:
            f.writelines(modified_lines)
        
        print(f"Modified topology written to: {output_top}")
        
    except Exception as e:
        print(f"Error in modify_topology_file: {e}")
        raise

def run_grompp(mdp_file, pdb_file, top_file, tpr_file):
    """Run gmx grompp to generate tpr file"""
    
    cmd = [
        "gmx", "grompp",
        "-f", str(mdp_file),
        "-c", str(pdb_file),
        "-p", str(top_file),
        "-o", str(tpr_file),
        "-maxwarn", "1000"
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"Successfully generated TPR file: {tpr_file}")
    else:
        print(f"Error running grompp: {result.stderr}")
        return False
    
    return True

def extract_binding_trajectory(u, start_frame, end_frame, output_xtc, output_pdb, 
                              top_file=None, mdp_file=None, chain_id=None, seg_id=None, prefix="binding", idx=0):
    """Extract subtrajectory between start_frame and end_frame"""
    
    # For our modified structure, we need to select atoms from BOTH the segid and the chainID
    # The DataFrame contains pairs like "A2,s" meaning segid A2 interacting with chainID s
    
    # Select atoms from the segid (our A1-A4 groups)
    seg_atoms = u.select_atoms(f"segid {seg_id}")
    
    # Select atoms from the chainID (the interacting partner)
    # For our A1-A4 segids, the chainID in the DataFrame needs to be mapped to the actual chainID
    actual_chain_id = chain_id
    chain_atoms = u.select_atoms(f"chainID {actual_chain_id}")
    
    # Combine both selections
    selected_atoms = seg_atoms + chain_atoms
    
    print(f"Selected {len(seg_atoms)} atoms from segid {seg_id} + {len(chain_atoms)} atoms from chainID {chain_id} = {len(selected_atoms)} total atoms")
    
    # First, go to first frame and write PDB
    u.trajectory[start_frame]
    # Set chainID to seg_id for all selected atoms
    for atom in selected_atoms:
        atom.chainID = seg_id
    with mda.Writer(str(output_pdb), selected_atoms.n_atoms) as pdb_writer:
        pdb_writer.write(selected_atoms)
    
    print(f"First frame PDB: {output_pdb}")
    
    # Then extract trajectory frames
    with mda.Writer(str(output_xtc), selected_atoms.n_atoms) as xtc_writer:
        for ts in u.trajectory[start_frame:end_frame+1]:
            xtc_writer.write(selected_atoms)
    
    print(f"Extracted trajectory: {start_frame}-{end_frame} -> {output_xtc}")
    
    # If topology file is provided, modify it and run grompp
    if top_file and mdp_file and chain_id and seg_id:
        print(f"Topology file provided: {top_file}")
        print(f"MDP file provided: {mdp_file}")
        print(f"Chain ID: {chain_id}, Seg ID: {seg_id}")
        
        output_top = output_pdb.with_suffix('.top')
        output_tpr = output_pdb.with_suffix('.tpr')
        
        try:
            modify_topology_file(top_file, [actual_chain_id], seg_id, output_top)
            success = run_grompp(mdp_file, str(output_pdb), str(output_top), str(output_tpr))
            if success:
                print(f"Generated TPR file: {output_tpr}")
            else:
                print(f"Warning: grompp failed for {output_pdb} - topology/coordinates mismatch is expected for subset extractions")
                print(f"Trajectory files are still valid: {output_xtc}")
        except Exception as e:
            print(f"Warning: Error processing topology/TPR for {output_pdb}: {e}")
            print(f"Trajectory files are still valid: {output_xtc}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Skipping topology processing - missing arguments: top_file={top_file}, mdp_file={mdp_file}, chain_id={chain_id}, seg_id={seg_id}")

def main():
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load structure and trajectory
    print(f"Loading structure: {args.structure}")
    print(f"Loading trajectory: {args.trajectory}")
    u = mda.Universe(args.structure, args.trajectory)
    
    # Load DataFrame
    print(f"Loading DataFrame: {args.dataframe}")
    df = pd.read_csv(args.dataframe)
    
    print(f"Found {len(df)} binding events")
    print(f"Trajectory has {len(u.trajectory)} frames")
    
    # Process each segid
    for resid_i, group in df.groupby('resid_i'):
        all_chain_ids = group['resid_j'].unique().tolist()
        actual_chain_ids = all_chain_ids
        
        first_row = group.iloc[0]
        start_frame_first = int(first_row['start_frame'])
        
        print(f"\nProcessing segid {resid_i} with chain_ids {actual_chain_ids}")
        
        # Select atoms
        seg_atoms = u.select_atoms(f"segid {resid_i}")
        seg_chain_ids = list(set(seg_atoms.chainIDs))
        all_chain_ids_to_keep = actual_chain_ids
        chain_suffix = '_'.join(sorted(all_chain_ids_to_keep))
        print(f"All chain_ids to keep: {all_chain_ids_to_keep}")
        
        # Create output filenames
        output_xtc = output_dir / f"{args.prefix}_{resid_i}_{chain_suffix}.xtc"
        output_pdb = output_dir / f"{args.prefix}_{resid_i}_{chain_suffix}.pdb"
        
        chain_atoms = sum(u.select_atoms(f"chainID {cid}") for cid in all_chain_ids_to_keep)
        selected_atoms = seg_atoms | chain_atoms
        
        print(f"Selected {len(seg_atoms)} atoms from segid {resid_i} + {len(chain_atoms)} atoms from chainIDs {all_chain_ids_to_keep} = {len(selected_atoms)} total atoms")
        
        # Go to first frame and write PDB
        u.trajectory[start_frame_first]
        # Set chainID to segid for all selected atoms
        for atom in selected_atoms:
            atom.chainID = resid_i
        with mda.Writer(str(output_pdb), selected_atoms.n_atoms) as pdb_writer:
            pdb_writer.write(selected_atoms)
        
        print(f"First frame PDB for {resid_i}: {output_pdb}")
        
        # Collect all unique frames across all binding events for this segid
        all_frames = set()
        for idx, row in group.iterrows():
            start_frame = int(row['start_frame'])
            end_frame = int(row['end_frame'])
            lifetime_ns = row['lifetime_ns']
            
            print(f"Event {idx}: {resid_i}-{row['resid_j']} (frames {start_frame}-{end_frame}, {lifetime_ns} ns)")
            
            # Add all frames from this event to the set
            for frame in range(start_frame, end_frame + 1):
                all_frames.add(frame)
        
        # Sort the unique frames
        unique_frames = sorted(list(all_frames))
        print(f"Total unique frames to extract: {len(unique_frames)} (from {len(group)} events)")
        
        # Then write each unique frame only once
        with mda.Writer(str(output_xtc), selected_atoms.n_atoms) as xtc_writer:
            for frame in unique_frames:
                # Validate frame range
                if frame >= len(u.trajectory):
                    print(f"WARNING: Frame {frame} exceeds trajectory length {len(u.trajectory)}")
                    continue
                
                u.trajectory[frame]
                xtc_writer.write(selected_atoms)
        
        print(f"Extracted trajectory for {resid_i}: {output_xtc}")
        
        # If topology file is provided, modify it and run grompp
        if args.topology and args.mdp_file:
            output_top = output_pdb.with_suffix('.top')
            output_tpr = output_pdb.with_suffix('.tpr')
            
            try:
                modify_topology_file(args.topology, all_chain_ids_to_keep, resid_i, output_top)
                success = run_grompp(args.mdp_file, str(output_pdb), str(output_top), str(output_tpr))
                if success:
                    print(f"Generated TPR for {resid_i}: {output_tpr}")
                else:
                    print(f"Warning: grompp failed for {output_pdb} - topology/coordinates mismatch is expected for subset extractions")
                    print(f"Trajectory files are still valid: {output_xtc}")
            except Exception as e:
                print(f"Warning: Error processing topology/TPR for {output_pdb}: {e}")
                print(f"Trajectory files are still valid: {output_xtc}")
                import traceback
                traceback.print_exc()
    
    print(f"\nCompleted processing {len(df['resid_i'].unique())} segids")
    print(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()
