"""
Read an array of csvs and concatenate them into a single csv file.

Arguments:
    -c: CSV files array
    -o: Output CSV file name
    -sep: Separator for the output CSV file (default: ',')
"""

import pandas as pd
import argparse

def parser_args():
    parser = argparse.ArgumentParser(description="Concatenate multiple CSV files into a single CSV file.")
    parser.add_argument("-c", "--csvs", nargs="+", required=True, help="Array of CSV files to concatenate.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file name.")
    parser.add_argument("-sep", "--separator", default=",", help="Separator for the output CSV file (default: ',').")
    return parser.parse_args()

def main():
    args = parser_args()
    csv_files = args.csvs
    output_file = args.output
    separator = args.separator

    # Read and concatenate all CSV files
    dataframes = []
    for csv_file in csv_files:
        # Read CSV with keep_default_na=False to treat "NA" as a string literal, not NaN
        df = pd.read_csv(csv_file, keep_default_na=False, na_values=[''])
        dataframes.append(df)
    
    combined_df = pd.concat(dataframes, ignore_index=True)

    # Save the combined dataframe to a new CSV file
    combined_df.to_csv(output_file, index=False, sep=separator)

    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()  
