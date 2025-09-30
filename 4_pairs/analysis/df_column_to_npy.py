"""
Extract columns from a DataFrame and save them as .npy files.

Arguments:
    -df <dataframe>: The input DataFrame from which to extract columns.
    -columns <list>: List of columns to be extracted.
    -output <str>: The output directory where the .npy files will be saved. By default the column names will be used as the file names.

Example usage:
    python df_column_to_npy.py -df input.csv -columns "col1" "col2" -output /path/to/output
"""

import pandas as pd
import argparse
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Extract columns from a DataFrame and save them as .npy files.")
    parser.add_argument("-df", "--dataframe", required=True, help="The input DataFrame from which to extract columns.")
    parser.add_argument("-columns", "--columns", nargs="+", required=True, help="List of columns to be extracted.")
    parser.add_argument("-output", "--output", default=".", help="The output directory where the .npy files will be saved.")
    return parser.parse_args()

def main():
    args = parse_args()
    df = pd.read_csv(args.dataframe)
    for col in args.columns:
        if col in df.columns:
            np.save(os.path.join(args.output, f"{col}.npy"), df[col].values)
        else:
            import sys
            sys.stderr.write(f"Column '{col}' not found in DataFrame.\n")
            print(f"Column '{col}' not found in DataFrame.")

if __name__ == "__main__":
    main()
