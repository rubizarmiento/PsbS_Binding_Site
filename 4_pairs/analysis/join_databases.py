"""
Reads input files and concatenates them into a single dataframe.

"""

import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Join multiple CSV databases into a single database.")
    parser.add_argument("-i", "--input_files", nargs='+', required=True, help="List of input CSV files to join.")
    parser.add_argument("-o", "--output_file", required=True, help="Output CSV file path.")
    return parser.parse_args()

def join_databases(input_files, output_file):
    dataframes = []
    for file in input_files:
        df = pd.read_csv(file)
        dataframes.append(df)
    combined_df = pd.concat(dataframes, ignore_index=True)
    combined_df.to_csv(output_file, index=False)

def main():
    args = parse_arguments()
    join_databases(args.input_files, args.output_file)
    print(f"Combined database saved to {args.output_file}")
if __name__ == "__main__":
    main()
