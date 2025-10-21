"""
Reads DataFrame with the format: resid_i,resid_j,start_frame,end_frame,frames,lifetime_ns
resid_1 are segids and resid_2 are chainIDs

Extract subtrajectories (.pdb, .xtc .tpr) per binding event

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
        cofactors_lines = []
        in_molecules_section = False

        for line in lines:
            stripped = line.strip()
            # Check if we're entering the molecules section
            if stripped == "[ molecules ]":
                in_molecules_section = True
                modified_lines.append(line)
                cofactors_lines.append(line)  # Append to cofactors too
                continue

            # If we're in molecules section and line is not empty
            elif in_molecules_section and stripped and not stripped.startswith(';'):
                # Check if line contains any of the chain_ids or seg_id
                seg_id_str = str(seg_id)
                if seg_id_str in stripped or any(stripped.startswith(f"chain_{cid}") for cid in chain_ids):
                    modified_lines.append(line)
                    print(f"Keeping protein molecule line: {stripped}")
                # If the line is already commented and contains a specific chain ID, keep it as is
                elif any(f"; chain_{cid}" in stripped for cid in chain_ids):
                    modified_lines.append(line)
                    cofactors_lines.append(f"{line}")
                else:
                    # Comment out the line
                    modified_lines.append(f";{line}")
                    cofactors_lines.append(f";{line}")
            else:
                modified_lines.append(line)
                cofactors_lines.append(line)  # Append all other lines to cofactors

        # Write modified topology
        with open(str(output_top), 'w') as f:
            f.writelines(modified_lines)
        
        output_top_cofactors = output_top.with_stem(output_top.stem + '_cofactors')
        with open(str(output_top_cofactors), 'w') as f:
            f.writelines(cofactors_lines)
        
        print(f"Modified topology written to: {output_top}")
        print(f"Modified topology written to: {output_top_cofactors}")

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
    cofactors_atoms = u.select_atoms(f"chainID {actual_chain_id} and resname CLA CLB CHL *MG* *HEM* *GG* DGD *SQ* *PG* W* *HG* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR")

    
    # Combine both selections
    selected_atoms = seg_atoms + chain_atoms
    print(f"Selected {len(seg_atoms)} atoms from segid {seg_id} + {len(chain_atoms)} atoms from chainID {chain_id} = {len(selected_atoms)} total atoms")
    print(f"The system contains {cofactors_atoms.n_residues} cofactors")
    
    # First, go to first frame and write PDB
    u.trajectory[start_frame]
    # Set chainID to seg_id for all selected atoms
    #DELETE#
    #for atom in selected_atoms:
    #    atom.chainID = seg_id

    with mda.Writer(str(output_pdb), selected_atoms.n_atoms) as pdb_writer:
        pdb_writer.write(selected_atoms)
    print(f"First frame PDB: {output_pdb}")
    # Then extract trajectory frames
    with mda.Writer(str(output_xtc), selected_atoms.n_atoms) as xtc_writer:
        for ts in u.trajectory[start_frame:end_frame+1]:
            xtc_writer.write(selected_atoms)
    print(f"Extracted trajectory: {start_frame}-{end_frame} -> {output_xtc}")

    if cofactors_atoms.n_atoms !=0:
        ouput_pdb_cofactors =  output_pdb.with_stem(output_pdb.stem + '_cofactors')
        with mda.Writer(str(ouput_pdb_cofactors), cofactors_atoms.n_atoms) as pdb_writer:
            pdb_writer.write(cofactors_atoms)
        print(f"First frame PDB: {ouput_pdb_cofactors}")

        ouput_xtc_cofactors =  output_xtc.with_stem(output_xtc.stem + '_cofactors')
        with mda.Writer(str(ouput_xtc_cofactors), cofactors_atoms.n_atoms) as xtc_writer:
            for ts in u.trajectory[start_frame:end_frame+1]:
                xtc_writer.write(cofactors_atoms)
        print(f"Extracted trajectory: {start_frame}-{end_frame} -> {ouput_xtc_cofactors}")
    else:
        print(f"Warning: The {ouput_pdb_cofactors} does not contain cofactors")

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
        all_chain_ids_to_keep = [cid.upper() for cid in actual_chain_ids]
        chain_suffix = '_'.join(sorted(all_chain_ids_to_keep))
        print(f"All chain_ids to keep: {all_chain_ids_to_keep}")
        
        # Create output filenames
        output_xtc = output_dir / f"{args.prefix}_{resid_i}_{chain_suffix}.xtc"
        output_pdb = output_dir / f"{args.prefix}_{resid_i}_{chain_suffix}.pdb"
        
        chain_atoms = sum(u.select_atoms(f"chainID {cid}") for cid in all_chain_ids_to_keep)
        cofactors_atoms = chain_atoms.select_atoms(f"resname CLA CLB CHL *MG* *HEM* *GG* DGD *SQ* *PG* W* *HG* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR")
        selected_atoms = seg_atoms | chain_atoms 
        
        print(f"Selected {len(seg_atoms)} atoms from segid {resid_i} + {len(chain_atoms)} atoms from chainIDs {all_chain_ids_to_keep} = {len(selected_atoms)} total atoms")
        print(f"Found {cofactors_atoms.n_residues} cofactors")
        # Go to first frame and write PDB
        u.trajectory[start_frame_first]

        with mda.Writer(str(output_pdb), selected_atoms.n_atoms) as pdb_writer:
            pdb_writer.write(selected_atoms)
        print(f"First frame PDB for {resid_i}: {output_pdb}")
        
        if cofactors_atoms.n_atoms != 0:
            ouput_pdb_cofactors =  output_pdb.with_stem(output_pdb.stem + '_cofactors')
            with mda.Writer(str(ouput_pdb_cofactors), cofactors_atoms.n_atoms) as pdb_writer:
                pdb_writer.write(cofactors_atoms)
            print(f"First frame PDB for {resid_i}: {ouput_pdb_cofactors}")


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

        if cofactors_atoms.n_atoms != 0:
            ouput_xtc_cofactors =  output_xtc.with_stem(output_xtc.stem + '_cofactors')

            with mda.Writer(str(ouput_xtc_cofactors), cofactors_atoms.n_atoms) as xtc_writer:
                for frame in unique_frames:
                    # Validate frame range
                    if frame >= len(u.trajectory):
                        print(f"WARNING: Frame {frame} exceeds trajectory length {len(u.trajectory)}")
                        continue
                    
                    u.trajectory[frame]
                    xtc_writer.write(cofactors_atoms)
            
            print(f"Extracted trajectory for {resid_i}: {ouput_xtc_cofactors}")

        
        # If topology file is provided, modify it and run grompp
        if args.topology and args.mdp_file:
            output_top = output_pdb.with_suffix('.top')
            output_tpr = output_pdb.with_suffix('.tpr')
            
            try:
                modify_topology_file(args.topology, all_chain_ids_to_keep, resid_i, output_top)
                
                # Determine cases based on cofactors
                selected_cofactors = seg_atoms
                if selected_cofactors.n_atoms != 0:
                    cases = ["", "_cofactors"]
                else:
                    cases = [""]
                
                for case in cases:
                    pdb = output_pdb.with_stem(output_pdb.stem + case)
                    top = output_top.with_stem(output_top.stem + case)
                    tpr = output_tpr.with_stem(output_tpr.stem + case)
                    
                    success = run_grompp(args.mdp_file, str(pdb), str(top), str(tpr))
                    if success:
                        print(f"Generated TPR for {resid_i} ({case}): {tpr}")
                    else:
                        print(f"Warning: grompp failed for {pdb} - topology/coordinates mismatch is expected for subset extractions")
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
