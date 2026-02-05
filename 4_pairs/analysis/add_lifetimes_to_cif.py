#!/usr/bin/env python3
"""
Add lifetime B-factors to a CIF file from CSV data.

This script reads a structure file (PDB or mmCIF), loads lifetime data from a CSV file,
and assigns B-factors based on residue lifetimes. The output is written in mmCIF format
to support B-factors > 999.

Usage:
------
python add_lifetimes_to_cif.py -f <structure_file> -sel <selection> -csv <lifetime_csv> -o <output_cif>

Arguments:
----------
-f, --file : str
    Input structure file (PDB or mmCIF format).
-sel, --selection : str
    Atom selection string (MDAnalysis syntax). Examples:
    - "protein"
    - "segid A B C"
    - "resid 1:100"
    - "name CA CB"
-csv, --csv : str
    Path to lifetime CSV file containing 'resid' and 'sum_ns' columns.
-o, --output : str
    Output mmCIF file path.

CSV Format:
-----------
The CSV file must contain at least two columns:
- resid : int - Residue ID
- sum_ns : float - Lifetime value in nanoseconds (used as B-factor)

If no matching residue is found in the CSV, B-factor defaults to 0.

Examples:
---------
# Basic protein structure
python add_lifetimes_to_cif.py -f structure.pdb -sel "protein" \\
    -csv lifetimes.csv -o output.cif

# Specific segid
python add_lifetimes_to_cif.py -f structure.pdb -sel "segid A1 A2" \\
    -csv lifetimes.csv -o output.cif

# With chain selection
python add_lifetimes_to_cif.py -f structure.pdb -sel "segid PSII and protein" \\
    -csv lifetimes.csv -o output.cif

Author: Rubi Zarmiento-Garcia
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import MDAnalysis as mda
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO

# Suppress MDAnalysis warnings
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')


def load_lifetime_data(csv_file, log_transform=False):
    """
    Load lifetime data from CSV file.
    
    Parameters
    ----------
    csv_file : str
        Path to CSV file with 'resid' and lifetime columns.
    log_transform : bool
        If True, apply log10 transformation to lifetime values.
        
    Returns
    -------
    dict
        Mapping of residue ID to B-factor (lifetime value).
        
    Raises
    ------
    FileNotFoundError
        If CSV file doesn't exist.
    ValueError
        If required columns are missing.
    """
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    
    df = pd.read_csv(csv_file)
    
    # Check for required columns
    if 'resid' not in df.columns:
        raise ValueError(f"CSV file missing 'resid' column. Available columns: {list(df.columns)}")
    
    # Try different possible column names for lifetime data
    lifetime_col = None
    for col in ['sum_ns', 'lifetime_ns', 'total_ns', 'occupancy', 'sum_time']:
        if col in df.columns:
            lifetime_col = col
            break
    
    if lifetime_col is None:
        raise ValueError(
            f"CSV file missing lifetime column. Expected one of: "
            f"['sum_ns', 'lifetime_ns', 'total_ns', 'occupancy', 'sum_time']. "
            f"Available columns: {list(df.columns)}"
        )
    
    # Create mapping
    lifetime_values = df[lifetime_col].values
    
    # Print before transformation stats
    print(f"Loaded {len(df)} residue-bfactor mappings from {csv_file}")
    print(f"  Using column: {lifetime_col}")
    print(f"  Before transformation: min={lifetime_values.min():.2f}, max={lifetime_values.max():.2f}")
    
    # Apply log transformation if requested
    if log_transform:
        lifetime_values = np.log10(lifetime_values + 1)
        print(f"  After log10 transformation: min={lifetime_values.min():.2f}, max={lifetime_values.max():.2f}")
    
    resid_to_bfactor = dict(zip(df['resid'], lifetime_values))
    
    return resid_to_bfactor


def assign_bfactors(structure, resid_to_bfactor):
    """
    Assign B-factors to atoms in structure based on residue mapping.
    
    Parameters
    ----------
    structure : Bio.PDB.Structure
        Structure to modify.
    resid_to_bfactor : dict
        Mapping of residue ID to B-factor value.
        
    Returns
    -------
    Bio.PDB.Structure
        Modified structure with B-factors assigned.
    """
    atoms_updated = 0
    for residue in structure.get_residues():
        residue_id = residue.get_id()[1]
        bfactor = resid_to_bfactor.get(residue_id, 0.0)
        
        for atom in residue.get_atoms():
            atom.set_bfactor(bfactor)
            atoms_updated += 1
            
            # Fix unknown elements
            if atom.element == "X":
                atom.element = "C"
    
    print(f"Assigned B-factors to {atoms_updated} atoms")
    return structure


def save_to_mmcif(structure, output_file):
    """
    Save structure to mmCIF format.
    
    Parameters
    ----------
    structure : Bio.PDB.Structure
        Structure to save.
    output_file : str
        Output file path.
    """
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)
    print(f"Saved structure to: {output_file}")


def main():
    """
    Main function: parse arguments and process structure.
    """
    parser = argparse.ArgumentParser(
        description="Add lifetime B-factors to a CIF file from CSV data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f structure.pdb -sel "protein" -csv lifetimes.csv -o output.cif
  %(prog)s -f structure.pdb -sel "segid A1 A2" -csv lifetimes.csv -o output.cif
        """
    )
    
    parser.add_argument('-f', '--file', required=True,
                        help='Input structure file (PDB or mmCIF)')
    parser.add_argument('-sel', '--selection', required=True,
                        help='Atom selection string (MDAnalysis syntax)')
    parser.add_argument('-csv', '--csv', required=True,
                        help='Lifetime CSV file with resid and sum_ns columns')
    parser.add_argument('-o', '--output', required=True,
                        help='Output mmCIF file path')
    parser.add_argument('--log_transform', action='store_true',
                        help='Apply log10 transformation to lifetime values')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isfile(args.file):
        raise FileNotFoundError(f"Structure file not found: {args.file}")
    
    print(f"Input structure: {args.file}")
    print(f"Selection: {args.selection}")
    print(f"Lifetime CSV: {args.csv}")
    print(f"Output: {args.output}")
    print()
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}\n")
    
    # Load lifetime data
    resid_to_bfactor = load_lifetime_data(args.csv, log_transform=args.log_transform)
    
    # Load structure with MDAnalysis and apply selection
    print(f"Loading structure with MDAnalysis...")
    u = mda.Universe(args.file)
    selected_atoms = u.select_atoms(args.selection)
    print(f"Selected {len(selected_atoms)} atoms")
    
    # Write selected atoms to temporary PDB
    temp_pdb = args.output.replace('.cif', '_temp.pdb')
    selected_atoms.write(temp_pdb)
    
    # Load with Bio.PDB for B-factor assignment
    print(f"Loading structure with Bio.PDB...")
    parser_bio = PDBParser(QUIET=True)
    structure = parser_bio.get_structure("struct", temp_pdb)
    
    # Assign B-factors
    structure = assign_bfactors(structure, resid_to_bfactor)
    
    # Save to mmCIF
    save_to_mmcif(structure, args.output)
    
    # Clean up temporary file
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
