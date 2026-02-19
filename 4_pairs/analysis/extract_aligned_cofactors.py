#!/usr/bin/env python3
"""
Extract cofactors from a CG structure after aligning it to an AT reference
using protein backbone atoms.

This script:
1. Loads a CG structure (protein + cofactors) from middle_cluster
2. Aligns it to the backmapped AT protein (final_aligned.pdb) using
   CG BB ↔ AT CA atoms as alignment anchors
3. Extracts cofactors from the aligned CG structure
4. Writes the aligned cofactors PDB

Usage:
    python extract_aligned_cofactors.py \
        -cg middle_cluster/1_n_s.pdb \
        -ref cg2at/1_n_s/FINAL/final_aligned.pdb \
        -sel_cofactors "resname CLA CLB ..." \
        -chains n s \
        -o cg2at/1_n_s_cofactors_cg.pdb

Author: Rubi Zarmiento Garcia
"""

import argparse
import warnings
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract aligned cofactors from CG structure."
    )
    parser.add_argument("-cg", required=True,
                        help="CG structure PDB (protein + cofactors, e.g. middle_cluster PDB).")
    parser.add_argument("-ref", required=True,
                        help="AT reference PDB (e.g. final_aligned.pdb).")
    parser.add_argument("-sel_cofactors", required=True,
                        help="MDAnalysis selection string for cofactors.")
    parser.add_argument("-chains", nargs="+", required=True,
                        help="Chain IDs to use for backbone alignment.")
    parser.add_argument("-o", required=True,
                        help="Output aligned cofactors PDB.")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load structures
    cg = mda.Universe(args.cg)
    ref = mda.Universe(args.ref)

    chain_str = " ".join(args.chains)
    sel_mobile = f"name BB and chainID {chain_str}"
    sel_ref = f"name CA and chainID {chain_str}"

    mobile_sel = cg.select_atoms(sel_mobile)
    ref_sel = ref.select_atoms(sel_ref)

    print(f"CG protein BB:  {len(mobile_sel)} atoms (chainID {chain_str})")
    print(f"AT protein CA:  {len(ref_sel)} atoms (chainID {chain_str})")

    # Handle atom count mismatch by using common resids
    if len(mobile_sel) != len(ref_sel):
        mobile_resids = set(mobile_sel.resids)
        ref_resids = set(ref_sel.resids)
        common = sorted(mobile_resids & ref_resids)
        print(f"Atom count mismatch, using {len(common)} common residues")

        # Filter to residues with matching atom counts
        matching = []
        for resid in common:
            n_mob = len(cg.select_atoms(f"{sel_mobile} and resid {resid}"))
            n_ref = len(ref.select_atoms(f"{sel_ref} and resid {resid}"))
            if n_mob == n_ref:
                matching.append(resid)

        resid_str = " ".join(map(str, matching))
        mobile_sel = cg.select_atoms(f"{sel_mobile} and resid {resid_str}")
        ref_sel = ref.select_atoms(f"{sel_ref} and resid {resid_str}")
        print(f"After filtering: mobile={len(mobile_sel)}, ref={len(ref_sel)}")

    if len(mobile_sel) != len(ref_sel):
        raise ValueError(
            f"Cannot align: {len(mobile_sel)} mobile vs {len(ref_sel)} ref atoms"
        )

    # Compute transformation
    mobile_com = mobile_sel.center_of_mass()
    ref_com = ref_sel.center_of_mass()
    R, rmsd = rotation_matrix(
        mobile_sel.positions - mobile_com,
        ref_sel.positions - ref_com
    )
    print(f"Alignment RMSD: {rmsd:.3f} Å")

    # Apply transformation to ALL atoms in the CG structure
    cg.atoms.translate(-mobile_com)
    cg.atoms.rotate(R)
    cg.atoms.translate(ref_com)

    # Extract and write cofactors
    cofactors = cg.select_atoms(args.sel_cofactors)
    print(f"Cofactors selected: {len(cofactors)} atoms")
    if len(cofactors) == 0:
        print("No cofactors found, skipping output.")
        return
    cofactors.write(args.o)
    print(f"Wrote aligned cofactors to {args.o}")


if __name__ == "__main__":
    main()
