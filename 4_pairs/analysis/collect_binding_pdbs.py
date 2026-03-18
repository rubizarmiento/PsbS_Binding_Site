"""
Collect and align binding site PDBs from cg2at output directories.

For each binding mode basename:
  1. Loads final_aligned.pdb (AA backmapped protein).
  2. Aligns it to a reference PDB via CA atoms (chain IDs parsed from basename).
  3. Applies the SAME rigid-body transform to the matching cofactors_cg.pdb.
  4. Writes protein -> binding_pdbs/<basename>.pdb
            cofactors -> binding_pdbs/<basename>_cofactors.pdb

Usage:
    python collect_binding_pdbs.py \
        -cg2at_dir  .../cg2at \
        -ref        .../psii_with_cofactors_aa.pdb \
        -odir       .../binding_pdbs \
        [-override  4_7_8=7  9_c_s_z=c]
"""

import argparse
import glob
import os
import warnings

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.align import rotation_matrix

warnings.filterwarnings("ignore")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-cg2at_dir", required=True,
                   help="Directory containing cg2at/<basename>/FINAL/final_aligned.pdb trees")
    p.add_argument("-ref", required=True,
                   help="AA reference PDB (e.g. psii_with_cofactors_aa.pdb)")
    p.add_argument("-odir", required=True,
                   help="Output directory for aligned PDBs")
    p.add_argument("-override", nargs="*", default=[],
                   metavar="BASENAME=CHAINS",
                   help="Per-basename chain override for alignment, e.g. 4_7_8=7")
    return p.parse_args()


def parse_overrides(override_list):
    """Parse ['4_7_8=7', '9_c_s_z=c'] -> {'4_7_8': '7', '9_c_s_z': 'c'}"""
    d = {}
    for entry in override_list:
        basename, chains = entry.split("=", 1)
        d[basename.strip()] = chains.strip()
    return d


def chains_from_basename(basename):
    """'4_7_8' -> ['7', '8'],  '1_n_s' -> ['n', 's']"""
    parts = basename.split("_")
    return parts[1:]   # drop the leading numeric ID


def align_and_collect(final_pdb, ref, align_chains, odir, basename):
    """
    Align final_pdb to ref on CA of align_chains.
    Returns (translation_vector, rotation_matrix, mobile_old_com, mobile_new_com).
    Writes <odir>/<basename>.pdb.
    """
    mobile = mda.Universe(final_pdb)

    chains_sel = " ".join(align_chains)
    mob_sel = mobile.select_atoms(f"name CA and chainID {chains_sel}")
    ref_sel = ref.select_atoms(f"name CA and chainID {chains_sel}")

    if len(mob_sel) == 0:
        raise ValueError(f"No CA atoms found in mobile for chainID {chains_sel}")
    if len(ref_sel) == 0:
        raise ValueError(f"No CA atoms found in reference for chainID {chains_sel}")
    if len(mob_sel) != len(ref_sel):
        raise ValueError(
            f"Atom count mismatch: mobile CA={len(mob_sel)}, ref CA={len(ref_sel)} "
            f"for chains {chains_sel}"
        )

    old_com = mob_sel.center_of_mass()
    ref_com = ref_sel.center_of_mass()

    mob0 = mob_sel.positions - old_com
    ref0 = ref_sel.positions - ref_com
    R, rmsd = rotation_matrix(mob0, ref0)

    # Apply to whole mobile universe
    mobile.atoms.translate(-old_com)
    mobile.atoms.rotate(R)
    mobile.atoms.translate(ref_com)

    out_pdb = os.path.join(odir, f"{basename}.pdb")
    mobile.atoms.write(out_pdb)

    new_com = mobile.select_atoms(f"name CA and chainID {chains_sel}").center_of_mass()
    print(f"  protein RMSD={rmsd:.2f} A  center: {np.round(old_com,1)} -> {np.round(ref_com,1)}")
    return R, old_com, ref_com


def transform_cofactors(cofactors_src, R, old_com, new_com, odir, basename):
    """Apply the same rigid transform to cofactors and write to odir."""
    cof = mda.Universe(cofactors_src)
    before = cof.atoms.center_of_mass()
    cof.atoms.translate(-old_com)
    cof.atoms.rotate(R)
    cof.atoms.translate(new_com)
    after = cof.atoms.center_of_mass()
    out = os.path.join(odir, f"{basename}_cofactors.pdb")
    cof.atoms.write(out)
    print(f"  cofactors center: {np.round(before,1)} -> {np.round(after,1)}")
    return out


def main():
    args = parse_args()
    overrides = parse_overrides(args.override)
    os.makedirs(args.odir, exist_ok=True)

    # Remove old PDBs
    for f in glob.glob(os.path.join(args.odir, "*.pdb")):
        os.remove(f)

    ref = mda.Universe(args.ref)

    pattern = os.path.join(args.cg2at_dir, "*/FINAL/final_aligned.pdb")
    final_pdbs = sorted(glob.glob(pattern))
    if not final_pdbs:
        raise FileNotFoundError(f"No final_aligned.pdb files found in {args.cg2at_dir}")

    print(f"Found {len(final_pdbs)} binding modes\n")

    for final_pdb in final_pdbs:
        basename = os.path.basename(os.path.dirname(os.path.dirname(final_pdb)))
        all_chains = chains_from_basename(basename)

        if basename in overrides:
            align_chains = overrides[basename].split()
            print(f"{basename}: aligning on chain(s) {align_chains} "
                  f"(override; all chains: {all_chains})")
        else:
            align_chains = all_chains
            print(f"{basename}: aligning on chain(s) {align_chains}")

        try:
            R, old_com, new_com = align_and_collect(
                final_pdb, ref, align_chains, args.odir, basename
            )
        except Exception as e:
            print(f"  ERROR during protein alignment: {e}")
            print(f"  Copying {final_pdb} without alignment")
            import shutil
            shutil.copy(final_pdb, os.path.join(args.odir, f"{basename}.pdb"))
            R, old_com, new_com = np.eye(3), np.zeros(3), np.zeros(3)

        cofactors_src = os.path.join(args.cg2at_dir, f"{basename}_cofactors_cg.pdb")
        if os.path.isfile(cofactors_src):
            try:
                transform_cofactors(cofactors_src, R, old_com, new_com, args.odir, basename)
            except Exception as e:
                print(f"  ERROR during cofactor transform: {e}")
                import shutil
                shutil.copy(cofactors_src,
                            os.path.join(args.odir, f"{basename}_cofactors.pdb"))
        else:
            print(f"  no cofactors_cg.pdb found, skipping")

        print()

    print("Done.")


if __name__ == "__main__":
    main()
