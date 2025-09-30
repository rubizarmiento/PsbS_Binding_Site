#!/usr/bin/env python3
"""
Extract DSSP secondary structure information for each chain in a PDB file.

Usage:
    python extract_dssp.py --pdb /path/to/pdb/file.pdb --output_dir /path/to/output/dir

This script uses MDAnalysis to:
1. Load the PDB file
2. For each chain, run DSSP analysis
3. Save the DSSP sequence as a simple text file
"""

import argparse
import os
from pathlib import Path
import sys

import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP


def extract_dssp_per_chain(pdb_file, output_dir):
    """
    Extract DSSP secondary structure for each chain in the PDB file.

    Parameters:
    pdb_file: Path to the PDB file
    output_dir: Directory to save DSSP files
    """
    print(f"Loading PDB file: {pdb_file}")

    # Load the universe
    u = mda.Universe(pdb_file)

    # Get unique chains
    chains = set(u.atoms.chainIDs)
    chains = sorted([c for c in chains if c.strip()])  # Remove empty chains

    print(f"Found {len(chains)} chains: {chains}")

    for chain_id in chains:
        print(f"Processing chain {chain_id}...")

        # Select atoms for this chain (protein atoms only)
        chain_atoms = u.select_atoms(f"chainID {chain_id} and protein")

        if len(chain_atoms) == 0:
            print(f"Warning: No protein atoms found for chain {chain_id}")
            continue

        # Run DSSP analysis
        try:
            dssp = DSSP(chain_atoms)
            dssp.run()

            # Extract the DSSP sequence
            dssp_sequence = ''.join(dssp.results.dssp[0])

            # Save to file
            output_file = os.path.join(output_dir, f"5XNL_chain_{chain_id}.dssp")
            with open(output_file, 'w') as f:
                f.write(f">5XNL_chain_{chain_id}\n")
                f.write(dssp_sequence + "\n")

            print(f"Saved DSSP for chain {chain_id} to {output_file}")
            print(f"  Length: {len(dssp_sequence)} residues")
            print(f"  Sequence: {dssp_sequence[:50]}{'...' if len(dssp_sequence) > 50 else ''}")

        except Exception as e:
            print(f"Error processing chain {chain_id}: {e}")
            continue

    print("DSSP extraction completed!")


def main():
    parser = argparse.ArgumentParser(description="Extract DSSP secondary structure per chain from PDB file")
    parser.add_argument('--pdb', required=True, help='Path to PDB file')
    parser.add_argument('--output_dir', required=True, help='Output directory for DSSP files')

    args = parser.parse_args()

    # Check if PDB file exists
    if not os.path.exists(args.pdb):
        print(f"Error: PDB file {args.pdb} does not exist")
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Extract DSSP
    extract_dssp_per_chain(args.pdb, args.output_dir)


if __name__ == "__main__":
    main()