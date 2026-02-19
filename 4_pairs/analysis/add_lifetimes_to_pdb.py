#!/usr/bin/env python3
"""
Add lifetime B-factors to a PDB file and append CONECT records.

This script:
1. Loads a structure PDB (with correct chain IDs/resids) using MDAnalysis
2. Assigns B-factors from CSV lifetime data
3. Writes a PDB via MDAnalysis
4. Appends CONECT records from a separate PDB (e.g. from gmx trjconv -conect)

Usage:
------
python add_lifetimes_to_pdb.py -f <structure.pdb> -conect <conect.pdb> -csv <lifetime.csv> [-csv ...] -o <output.pdb>

Author: Rubi Zarmiento-Garcia
"""

import os
import argparse
import warnings
import numpy as np
import pandas as pd
import MDAnalysis as mda

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')


def load_lifetime_data(csv_files, log_transform=False):
    """
    Load lifetime data from one or more CSV files.

    Parameters
    ----------
    csv_files : list of str
        Path(s) to CSV file(s) with 'resid' and lifetime columns.
        When multiple files are provided, each must have a 'chainID_i' column.
    log_transform : bool
        If True, apply log10 transformation to lifetime values.

    Returns
    -------
    dict
        Mapping of resid (int) or (chainID, resid) tuple to B-factor value.
    """
    use_chain_key = len(csv_files) > 1
    resid_to_bfactor = {}

    for csv_file in csv_files:
        if not os.path.isfile(csv_file):
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        df = pd.read_csv(csv_file)

        if 'resid' not in df.columns:
            raise ValueError(f"CSV file missing 'resid' column. Available columns: {list(df.columns)}")

        if use_chain_key and 'chainID_i' not in df.columns:
            raise ValueError(
                f"Multiple CSV files provided but '{csv_file}' is missing 'chainID_i' column. "
                f"Available columns: {list(df.columns)}"
            )

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

        lifetime_values = df[lifetime_col].values.copy()

        print(f"Loaded {len(df)} residue-bfactor mappings from {csv_file}")
        print(f"  Using column: {lifetime_col}")
        print(f"  Before transformation: min={lifetime_values.min():.2f}, max={lifetime_values.max():.2f}")

        if log_transform:
            lifetime_values = np.log10(lifetime_values + 1)
            print(f"  After log10 transformation: min={lifetime_values.min():.2f}, max={lifetime_values.max():.2f}")

        if use_chain_key:
            chain_id = df['chainID_i'].iloc[0]
            print(f"  Chain: {chain_id}")
            for resid, bfactor in zip(df['resid'], lifetime_values):
                resid_to_bfactor[(chain_id, resid)] = bfactor
        else:
            resid_to_bfactor = dict(zip(df['resid'], lifetime_values))

    if use_chain_key:
        print(f"\nTotal: {len(resid_to_bfactor)} (chain, resid) mappings from {len(csv_files)} CSV files")

    return resid_to_bfactor


def assign_bfactors_mda(universe, resid_to_bfactor):
    """
    Assign B-factors to atoms using MDAnalysis, keyed by (chainID, resid) or resid.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis Universe with the structure.
    resid_to_bfactor : dict
        Mapping of resid or (chainID, resid) to B-factor value.

    Returns
    -------
    int
        Number of atoms updated.
    """
    use_chain_key = any(isinstance(k, tuple) for k in resid_to_bfactor)
    atoms_updated = 0

    for residue in universe.residues:
        resid = residue.resid
        chain_id = residue.atoms[0].chainID if hasattr(residue.atoms[0], 'chainID') else ''

        if use_chain_key:
            bfactor = resid_to_bfactor.get((chain_id, resid),
                                            resid_to_bfactor.get(resid, 0.0))
        else:
            bfactor = resid_to_bfactor.get(resid, 0.0)

        # Clamp to PDB max
        bfactor = min(bfactor, 999.99)

        for atom in residue.atoms:
            atom.tempfactor = bfactor
            atoms_updated += 1

    print(f"Assigned B-factors to {atoms_updated} atoms")
    return atoms_updated


def extract_conect_lines(conect_pdb):
    """
    Extract CONECT lines from a PDB file.

    Parameters
    ----------
    conect_pdb : str
        Path to PDB file with CONECT records (e.g. from gmx trjconv -conect).

    Returns
    -------
    list of str
        CONECT lines.
    """
    conect_lines = []
    with open(conect_pdb, 'r') as f:
        for line in f:
            if line.startswith('CONECT'):
                conect_lines.append(line)
    print(f"Extracted {len(conect_lines)} CONECT records from {conect_pdb}")
    return conect_lines


def append_conect_to_pdb(pdb_file, conect_lines):
    """
    Append CONECT lines to a PDB file (before END if present).

    Parameters
    ----------
    pdb_file : str
        Path to PDB file to modify in-place.
    conect_lines : list of str
        CONECT lines to append.
    """
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Find the last END line and insert CONECT before it
    end_idx = None
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip() == 'END':
            end_idx = i
            break

    if end_idx is not None:
        new_lines = lines[:end_idx] + conect_lines + lines[end_idx:]
    else:
        new_lines = lines + conect_lines

    with open(pdb_file, 'w') as f:
        f.writelines(new_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Add lifetime B-factors to a PDB file and append CONECT records.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-f', '--file', required=True,
                        help='Input structure PDB (with correct chain IDs and resids)')
    parser.add_argument('-conect', '--conect', default=None,
                        help='PDB file with CONECT records (e.g. from gmx trjconv -conect)')
    parser.add_argument('-csv', '--csv', required=True, nargs='+',
                        help='Lifetime CSV file(s) with resid and sum_ns columns.')
    parser.add_argument('-o', '--output', required=True,
                        help='Output PDB file path')
    parser.add_argument('--log_transform', action='store_true',
                        help='Apply log10 transformation to lifetime values')

    args = parser.parse_args()

    if not os.path.isfile(args.file):
        raise FileNotFoundError(f"Input PDB file not found: {args.file}")

    print(f"Input structure: {args.file}")
    print(f"CONECT source:   {args.conect}")
    print(f"Lifetime CSV(s): {args.csv}")
    print(f"Output:          {args.output}")
    print()

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # 1. Load lifetime data
    resid_to_bfactor = load_lifetime_data(args.csv, log_transform=args.log_transform)

    # 2. Load structure with MDAnalysis and assign B-factors
    print(f"Loading structure with MDAnalysis...")
    u = mda.Universe(args.file)
    print(f"Loaded {len(u.atoms)} atoms")
    assign_bfactors_mda(u, resid_to_bfactor)

    # 3. Write PDB with MDAnalysis (correct chains, resids, B-factors)
    u.atoms.write(args.output)
    print(f"Wrote PDB: {args.output}")

    # 4. Append CONECT records from trjconv output
    if args.conect and os.path.isfile(args.conect):
        conect_lines = extract_conect_lines(args.conect)
        if conect_lines:
            append_conect_to_pdb(args.output, conect_lines)
            print(f"Appended {len(conect_lines)} CONECT records")
    elif args.conect:
        print(f"WARNING: CONECT file not found: {args.conect}")

    print("\nDone!")


if __name__ == "__main__":
    main()
