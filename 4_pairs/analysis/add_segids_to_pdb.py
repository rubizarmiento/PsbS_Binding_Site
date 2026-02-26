"""
Adds segids to a PDB file based on helix definitions from a YAML file.

Reads the PDB file with MDAnalysis and a helix definitions YAML file.
For each chain defined in the YAML, residues whose resid falls within
any helix range get segid "H" (helix), all other residues get segid "L" (loop).
Chains not present in the YAML are left unchanged.

helix definitions yaml format:
  chain_ID:
    helix_name:
      start: residue_number  # int
      end: residue_number    # int
      ...

Usage:
  python3 add_segids_to_pdb.py \
    --input_pdb path/to/input.pdb \
    --helix_yaml path/to/helix_definitions.yaml \
    --output_pdb path/to/output.pdb
"""

import argparse
import yaml
import MDAnalysis as mda


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Add segids (H=helix, L=loop) to a PDB based on helix definitions."
    )
    parser.add_argument(
        "-i", "--input_pdb", required=True, help="Path to the input PDB file."
    )
    parser.add_argument(
        "-y", "--helix_yaml", required=True,
        help="Path to the helix definitions YAML file."
    )
    parser.add_argument(
        "-o", "--output_pdb", required=True, help="Path to the output PDB file."
    )
    return parser.parse_args()


def read_helix_yaml(yaml_file):
    """
    Read helix definitions from a YAML file.

    Returns
    -------
    dict
        {chain_id: [(start, end), ...], ...}
    """
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

    helix_ranges = {}
    for chain_id, helices in data.items():
        chain_id = str(chain_id)
        ranges = []
        for helix_name, info in helices.items():
            ranges.append((int(info["start"]), int(info["end"]), str(helix_name)))
        helix_ranges[chain_id] = ranges

    return helix_ranges


def compute_new_segids(universe, helix_ranges):
    """
    Compute new segid values for all atoms in the universe.

    For each chain in helix_ranges:
      - Residues within a helix range get segid "H"
      - All other residues in that chain get segid "L"
    Chains not in helix_ranges keep their original segid.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The loaded structure.
    helix_ranges : dict
        {chain_id: [(start, end), ...], ...}

    Returns
    -------
    numpy.ndarray
        Array of segid strings, one per atom.
    """
    new_segids = universe.atoms.segids.copy()

    for chain_id, ranges in helix_ranges.items():
        chain_sel = universe.select_atoms(f"chainID {chain_id}")
        if len(chain_sel) == 0:
            print(f"Warning: No atoms found for chainID {chain_id}, skipping.")
            continue

        # Default: all atoms in this chain get "L"
        new_segids[chain_sel.indices] = "L"

        # Atoms in helix ranges get the helix name as segid
        for start, end, helix_name in ranges:
            helix_sel = universe.select_atoms(
                f"chainID {chain_id} and resid {start}:{end}"
            )
            if len(helix_sel) > 0:
                new_segids[helix_sel.indices] = helix_name

    return new_segids


def write_pdb_with_segids(input_pdb, output_pdb, new_segids):
    """
    Write a PDB file where the segid field (columns 73-76) is replaced
    with the values from new_segids array.

    Parameters
    ----------
    input_pdb : str
        Path to the input PDB file.
    output_pdb : str
        Path to the output PDB file.
    new_segids : array-like
        Array of segid strings, one per ATOM/HETATM record.
    """
    atom_idx = 0
    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                if atom_idx < len(new_segids):
                    segid = new_segids[atom_idx]
                    # PDB format: columns 73-76 are the segment ID (1-indexed)
                    line = line[:72] + f"{segid:<4s}" + line[76:]
                    atom_idx += 1
                fout.write(line)
            else:
                fout.write(line)

    print(f"Wrote {atom_idx} atoms to {output_pdb}")


def main():
    args = parse_arguments()

    # Load PDB
    universe = mda.Universe(args.input_pdb)
    print(f"Loaded {args.input_pdb}: {universe.atoms.n_atoms} atoms, "
          f"{universe.residues.n_residues} residues")

    # Read helix definitions
    helix_ranges = read_helix_yaml(args.helix_yaml)
    print(f"Read helix definitions for {len(helix_ranges)} chains from {args.helix_yaml}")

    # Compute new segids
    new_segids = compute_new_segids(universe, helix_ranges)

    # Write output PDB with updated segids
    write_pdb_with_segids(args.input_pdb, args.output_pdb, new_segids)


if __name__ == "__main__":
    main()
