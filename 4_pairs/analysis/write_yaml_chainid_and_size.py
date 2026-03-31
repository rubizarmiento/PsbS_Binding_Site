"""
Write a YAML file mapping each chain ID in a PDB to its size (number of residues).

The output YAML has the format:
  chain_id:
    n_residues: <int>

Usage:
  python3 write_yaml_chainid_and_size.py -f reference.pdb -o chain_sizes.yaml
"""

import argparse
import yaml
import MDAnalysis as mda
import warnings

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write YAML mapping chain IDs to their number of residues."
    )
    parser.add_argument(
        "-f", "--pdb", required=True,
        help="Input PDB file (reference structure)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output YAML file."
    )
    parser.add_argument(
        "--selection", default="name CA",
        help="MDAnalysis selection to count residues. Default: 'name CA' "
             "(one atom per residue for proteins). Use 'all' for coarse-grained."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    u = mda.Universe(args.pdb)

    chain_ids = sorted(set(u.atoms.chainIDs))

    chain_sizes = {}
    for cid in chain_ids:
        sel = u.select_atoms(f"chainID {cid} and {args.selection}")
        n = sel.n_atoms
        if n > 0:
            chain_sizes[cid] = {"n_residues": int(n)}

    with open(args.output, "w") as f:
        yaml.dump(chain_sizes, f, default_flow_style=False, sort_keys=False)

    print(f"Wrote {len(chain_sizes)} chains to {args.output}")


if __name__ == "__main__":
    main()
