"""
Reads an structure file and returns a list of non-protein residues.

Args:
-f: str: Path to the structure file.
-sel: str: Selection string to filter residues.

Returns:
-list: List[str]: List of non-protein residue names.
"""

from typing import List
import MDAnalysis as mda
from MDAnalysis.core.groups import ResidueGroup
from MDAnalysis.core.universe import Universe
import argparse
# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

# Standard protein residue names
PROTEIN_RESNAMES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

def parser():
    parser = argparse.ArgumentParser(description="Return non-protein residues from a structure file based on a selection string.")
    parser.add_argument("-f", type=str, required=True, help="Path to the structure file.")
    parser.add_argument("-sel", type=str, required=True, help="Selection string to filter residues.")
    return parser.parse_args()


def return_non_protein_residues(f: str, sel: str) -> List[str]:
    """
    Reads an structure file and returns a list of non-protein residues.

    Args:
    - f: str: Path to the structure file.
    - sel: str: Selection string to filter residues.

    Returns:
    - list: List[str]: List of non-protein residue names.
    """
    u: Universe = mda.Universe(f)
    selected_residues: ResidueGroup = u.select_atoms(sel).residues
    non_protein_residues = [
        res.resname for res in selected_residues
        if res.resname not in PROTEIN_RESNAMES
    ]
    return list(set(non_protein_residues))  # Return unique residue names

def check_if_file_exists(f: str) -> bool:
    """Check if the given file exists."""
    import os
    return os.path.isfile(f)

def main():
    args = parser()
    if not check_if_file_exists(args.f):
        print(f"File not found: {args.f}")
        exit(1)

    non_protein_residues = return_non_protein_residues(args.f, args.sel)
    print("Non-protein residues found:", non_protein_residues)

if __name__ == "__main__":
    main()

