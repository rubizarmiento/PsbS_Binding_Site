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

def normalize_columns(df):
    """Normalize column names from raw per-binding-mode format to summary format.
    
    Raw per-binding-mode CSVs have: resid, n_events, occupancy_pct, median_ns, sum_ns, p90_ns, chainID_i, chainID_j, resname_i
    Summary CSVs have: resid, resname, chain, sum_ns
    """
    rename_map = {}
    if 'chainID_i' in df.columns and 'chain' not in df.columns:
        rename_map['chainID_i'] = 'chain'
    if 'resname_i' in df.columns and 'resname' not in df.columns:
        rename_map['resname_i'] = 'resname'
    if rename_map:
        df = df.rename(columns=rename_map)
    return df

def transform_df(df):
    df = normalize_columns(df)
    half = len(df) // 2
    df1 = df.iloc[0:half, :].copy().reset_index(drop=True)
    df2 = df.iloc[half:, :].copy().reset_index(drop=True)
    
    # Sum only the numeric 'sum_ns' column from both halves
    sum_ns = df1['sum_ns'].values + df2['sum_ns'].values
    
    # Build result: keep resid/resname/chain from both halves, assign symmetric sum_ns
    result1 = df1[['resid', 'resname', 'chain']].copy()
    result1['sum_ns'] = sum_ns
    
    result2 = df2[['resid', 'resname', 'chain']].copy()
    result2['sum_ns'] = sum_ns
    
    result = pd.concat([result1, result2], ignore_index=True)
    
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




