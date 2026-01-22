"""
Reads the lifetimes*csv file for PsbS and sums the values of chain A and B. Then writes the csv file of the symmetric values.

Arguments:
    -csv:  PsbS lifetimes with residues from 1 to 424
    -o: Symmetrized PsbS lifetimes with residues from 1 to 424
"""

import pandas as pd
import argparse

def parser_args():
    parser = argparse.ArgumentParser(description='Symmetrize PsbS lifetimes by summing chain A and B values')
    parser.add_argument('-csv', type=str, help='Input PsbS lifetimes CSV file')
    parser.add_argument('-o', type=str, help='Output symmetrized lifetimes CSV file')
    args = parser.parse_args()
    return args

def transform_df(df):
    half = len(df) // 2
    df1 = df.iloc[0:half, :].copy()
    df2 = df.iloc[half:, :].copy()
    
    # Sum corresponding rows from both halves (all numeric columns)
    sum_ns = (df1.iloc[:,1:].values + df2.iloc[:,1:].values).sum(axis=1)
    
    # Keep original DataFrame structure, just replace with sum_ns column
    result = df.iloc[:,0].copy().to_frame()
    # Assign symmetric values: first half gets sum_ns, second half gets same sum_ns
    result['sum_ns'] = 0.0  # Initialize
    result.iloc[0:half, result.columns.get_loc('sum_ns')] = sum_ns
    result.iloc[half:, result.columns.get_loc('sum_ns')] = sum_ns
    
    return result

def main():
    args = parser_args()
    df = pd.read_csv(args.csv)
    df = transform_df(df)
    df.to_csv(args.o,index=False)
    #print("Input file:", args.csv)
    print("Output file saved to:", args.o)

if __name__ == "__main__":
    main()




