#!/usr/bin/env python3
"""
Merge per-chain cg2at de_novo PDB files into a single PDB.

Usage:
    python merge_cg2at_pdbs.py -pdbs chain7.pdb chain8.pdb chain9.pdb -o merged.pdb

Each input PDB is a cg2at de_novo output for a single chain.
The merged PDB concatenates all ATOM records and adds proper TER/END records.

Author: Rubi Zarmiento-Garcia
"""

import argparse
import os


def merge_pdbs(pdb_files, output):
    """Merge multiple PDB files into one, keeping only ATOM/HETATM records."""
    all_lines = []
    for pdb_file in pdb_files:
        if not os.path.isfile(pdb_file):
            print(f"WARNING: {pdb_file} not found, skipping.")
            continue
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    all_lines.append(line)
        # Add TER between chains
        if all_lines and not all_lines[-1].startswith('TER'):
            all_lines.append('TER\n')

    # Write merged PDB
    with open(output, 'w') as f:
        f.writelines(all_lines)
        if not all_lines[-1].startswith('END'):
            f.write('END\n')

    print(f"Merged {len(pdb_files)} PDBs -> {output} ({sum(1 for l in all_lines if l.startswith(('ATOM', 'HETATM')))} atoms)")


def main():
    parser = argparse.ArgumentParser(description="Merge per-chain cg2at PDBs into one.")
    parser.add_argument('-pdbs', nargs='+', required=True, help='Input PDB files to merge')
    parser.add_argument('-o', required=True, help='Output merged PDB file')
    args = parser.parse_args()
    merge_pdbs(args.pdbs, args.o)


if __name__ == "__main__":
    main()
