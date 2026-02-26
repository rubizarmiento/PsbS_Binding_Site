"""
Add one or more columns with fixed values to one or more CSV files.

Arguments:
    -c:           CSV files (one or more)
    -o:           Output CSV file(s). If a single file is given with multiple inputs,
                  all inputs are concatenated into that file. If one output per input
                  is given, each input is saved separately.
    -colname:     Column name(s) to add (array)
    -value:       Value(s) for each column (array, must match length of -colname)
    -sep:         Separator for the output CSV file (default: ',')
"""

import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add columns with fixed values to one or more CSV files."
    )
    parser.add_argument("-c", "--csvs", nargs="+", required=True,
                        help="Input CSV file(s).")
    parser.add_argument("-o", "--output", nargs="+", required=True,
                        help="Output CSV file(s). Provide one file to concatenate all "
                             "inputs, or one file per input to save separately.")
    parser.add_argument("-colname", "--colname", nargs="+", required=True,
                        help="Column name(s) to add.")
    parser.add_argument("-value", "--value", nargs="+", required=True,
                        help="Value(s) for each column (must match number of -colname).")
    parser.add_argument("-sep", "--separator", default=",",
                        help="Separator for the output CSV file (default: ',').")
    return parser.parse_args()


def add_columns(df, colnames, values):
    for colname, value in zip(colnames, values):
        df[colname] = value
    return df


def main():
    args = parse_args()

    if len(args.colname) != len(args.value):
        raise ValueError(
            f"Number of column names ({len(args.colname)}) must match "
            f"number of values ({len(args.value)})."
        )

    # Single output file: concatenate all inputs then save
    if len(args.output) == 1:
        dataframes = []
        for csv_file in args.csvs:
            df = pd.read_csv(csv_file, keep_default_na=False, na_values=[""])
            df = add_columns(df, args.colname, args.value)
            dataframes.append(df)
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(args.output[0], index=False, sep=args.separator)
        print(f"Output saved to {args.output[0]}")

    # One output per input
    elif len(args.output) == len(args.csvs):
        for csv_file, out_file in zip(args.csvs, args.output):
            df = pd.read_csv(csv_file, keep_default_na=False, na_values=[""])
            df = add_columns(df, args.colname, args.value)
            df.to_csv(out_file, index=False, sep=args.separator)
            print(f"Output saved to {out_file}")

    else:
        raise ValueError(
            "Provide either a single output file (to concatenate all inputs) "
            "or one output file per input CSV."
        )


if __name__ == "__main__":
    main()
