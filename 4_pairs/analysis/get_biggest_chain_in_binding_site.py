"""
Identify the biggest chain in each binding site and write to CSV.

Reads:
  - A YAML file with chain sizes (chain_id -> {n_residues: int}),
    produced by write_yaml_chainid_and_size.py.
  - A basenames_binding CSV with columns including 'unique_basename' and 'tag',
    where 'tag' is an underscore-separated list of chain IDs (e.g. "6_7", "c_s_z").

Writes a CSV with the columns:
  unique_basename, tag, biggest_chain, biggest_chain_n_residues

Usage:
  python3 get_biggest_chain_in_binding_site.py \
    -yaml chain_sizes.yaml \
    -csv  basenames_binding.csv \
    -ocsv basenames_biggest_chain.csv
"""

import argparse
import csv
import yaml


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find the biggest chain in each binding site."
    )
    parser.add_argument(
        "-yaml", required=True,
        help="YAML file mapping chain IDs to {n_residues: int}."
    )
    parser.add_argument(
        "-csv", required=True,
        help="Input CSV with 'unique_basename' and 'tag' columns."
    )
    parser.add_argument(
        "-ocsv", required=True,
        help="Output CSV file."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Load chain sizes
    with open(args.yaml, "r") as f:
        chain_sizes = yaml.safe_load(f)
    # Normalise keys to str
    chain_sizes = {str(k): v["n_residues"] for k, v in chain_sizes.items()}

    # Read input CSV
    with open(args.csv, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Process each binding site
    results = []
    for row in rows:
        basename = row["unique_basename"]
        tag = row["tag"]
        chains = tag.split("_")

        best_chain = None
        best_size = -1
        for c in chains:
            size = chain_sizes.get(c, 0)
            if size > best_size:
                best_size = size
                best_chain = c

        results.append({
            "unique_basename": basename,
            "tag": tag,
            "biggest_chain": best_chain if best_chain else "NA",
            "biggest_chain_n_residues": best_size if best_size > 0 else 0,
        })

    # Write output
    fieldnames = ["unique_basename", "tag", "biggest_chain", "biggest_chain_n_residues"]
    with open(args.ocsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {args.ocsv}")
    for r in results:
        print(f"  {r['unique_basename']}: chain {r['biggest_chain']} ({r['biggest_chain_n_residues']} residues)")


if __name__ == "__main__":
    main()
