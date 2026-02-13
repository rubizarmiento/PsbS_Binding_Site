"""
Sums a dataframe by groups and saves the result to a csv file.

Arguments:
-i: CSV input file path
-o: CSV output file path
-category: Row values to group by
-value: Column to sum

"""

import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Sum dataframe by groups and save to CSV.")
    parser.add_argument('-i', '--input', required=True, help='Input CSV file path')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file path')
    parser.add_argument('-c', '--category', required=True, help='Column name to group by')
    parser.add_argument('-v', '--value', required=True, help='Column name to sum')
    return parser.parse_args()

def group_and_sum_df(df, category, value):
    grouped_df = df.groupby(category)[value].sum().reset_index()
    return grouped_df

def main():
    args = parse_arguments()
    print("Input file:", args.input)
    # Read the input CSV file
    df = pd.read_csv(args.input)
    
    # Group and sum the dataframe
    result_df = group_and_sum_df(df, args.category, args.value)
    
    # Save the result to the output CSV file
    result_df.to_csv(args.output, index=False)
    print("Output file saved to:", args.output)

if __name__ == "__main__":
    main()

