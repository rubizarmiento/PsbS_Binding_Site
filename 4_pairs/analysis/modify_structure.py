#!/usr/bin/env python3
"""
Fixed modify_structure.py that correctly assigns segids and/or chainIDs to only selected atoms
"""

import MDAnalysis as mda
import argparse
import sys
import warnings
warnings.filterwarnings("ignore")

def parse_args():
    parser = argparse.ArgumentParser(description='Modify structure segids and/or chainIDs for selected atoms only')
    parser.add_argument('-f', '--file', type=str, required=True, help='Input PDB file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output PDB file')
    parser.add_argument('-sel', '--selection', type=str, required=True, help='MDAnalysis selection')
    parser.add_argument('-segid', '--segid', type=str, help='New segid for selected atoms (4 characters max)')
    parser.add_argument('-chain', '--chain', type=str, help='New chainID for selected atoms (1 character)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Validate that at least one modification option is provided
    if not args.segid and not args.chain:
        print("Error: Must specify either -segid or -chain (or both)")
        sys.exit(1)
    
    # Validate chain length
    if args.chain and len(args.chain) != 1:
        print("Error: Chain ID must be exactly 1 character")
        sys.exit(1)
    
    # Validate segid length
    if args.segid and len(args.segid) > 4:
        print("Error: Segid must be 4 characters or less")
        sys.exit(1)
    
    # Load universe
    u = mda.Universe(args.file)
    
    # Select atoms to modify
    selection = u.select_atoms(args.selection)
    print(f"Number of atoms selected: {len(selection)}")
    
    if len(selection) == 0:
        print("Warning: No atoms selected with the given criteria")
    
    # Get atom IDs that need to be changed
    atom_ids_to_change = set(selection.ids)
    
    # Read the input file and modify it line by line
    with open(args.file, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM')):
                # Extract atom ID from PDB line
                atom_id = int(line[6:11].strip())
                
                if atom_id in atom_ids_to_change:
                    new_line = line
                    
                    # Modify chainID if specified (column 22)
                    if args.chain:
                        new_line = new_line[:21] + args.chain + new_line[22:]
                    
                    # Modify segid if specified (columns 72-75)
                    if args.segid:
                        new_line = new_line[:72] + f"{args.segid:<4}" + new_line[76:]
                    
                    outfile.write(new_line)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
    
    print(f"Output file: {args.output}")
    if args.chain:
        print(f"Modified chainID to: {args.chain}")
    if args.segid:
        print(f"Modified segid to: {args.segid}")

if __name__ == "__main__":
    main()
