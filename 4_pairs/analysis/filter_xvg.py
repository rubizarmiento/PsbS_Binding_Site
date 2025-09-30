"""
Reads an xvg file and splits it by value.

Arguments:
    -i: Input xvg file
    -o: Output directory
    -prefix: Output file name prefix
    -c: Column index to split by, default is 1. Column index starts at 0.

Example:
python3 -i input.xvg -o output_dir -prefix split -c 1
"""
import pandas as pd
import io
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Process xvg files")
    parser.add_argument("-i", "--input", required=True, help="Input xvg file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-prefix", "--prefix", default="output", help="Output file name prefix")
    parser.add_argument("-c", "--column", type=int, default=1, help="Column index to split by, default is 1. Colindex start with 0.")
    return parser.parse_args()


def read_xvg(file,header=None):
    #Read the xvg file
    # Only comment characters used in xvg files
    comments = ["@", "#"]
    #Read the xvg file with readlines, generate a text variable, skip the lines that start with the special characters
    with open(file, "r") as f:
        lines = f.readlines()
    text = ""
    for i in lines:
        if not i.startswith(tuple(comments)):
            text += i
    df = pd.read_csv(io.StringIO(text), sep="\s+", header=header, skiprows=0, skip_blank_lines=True, skipinitialspace=True, engine="python")        
    df = pd.read_csv(io.StringIO(text), sep="\s+", comment="@", header=header, skiprows=0, skip_blank_lines=True, skipinitialspace=True, engine="python")        
    return df

def split_by_colvalue(df, colname):
    #Split the dataframe by the values in the specified column
    groups = df.groupby(colname)
    return groups

def split_colid_by_value(df, colid):
    #Split the dataframe by the values in the specified column
    groups = df.groupby(df.columns[colid])
    return groups

def main():
    args = parse_args()
    df = read_xvg(args.input)
    groups = split_colid_by_value(df, args.column)
    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)
    # Save each group to a separate file
    for name, group in groups:
        group.to_csv(os.path.join(args.output, f"{args.prefix}_{name}.xvg"), index=False, sep="\t")
        output_path = os.path.join(args.output, f"{args.prefix}_{name}.xvg")
        group.to_csv(output_path, index=False, sep="\t")

if __name__ == "__main__":
    main()
