"""
Read the lifetimes csv files and perform the following analysis

-max: Get the largest ocuppancy value, usefullt for heatplots
-lifetimes_dir: Directory where the lifetimes csv files are stored


"""


import MDAnalysis
import os
import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analyze lifetimes statistics from CSV files.')
    parser.add_argument('-lifetimes_dir', type=str, required=True, help='Directory where the lifetimes csv files are stored')
    parser.add_argument('-output', type=str, required=True, help='Dat file with statistics')
    args = parser.parse_args()
    return args

def get_max_occupancy(files):
    occupancies = []
    counter=0
    for f in files:
        df = pd.read_csv(f)
        # Check if 'median_ns' column exists
        if 'median_ns' in df.columns:
            occupancies.append(df['median_ns'].max())
        else:
            print(f"Warning: 'median_ns' column not found in {f}")
            counter+=1

    print(f"Files without 'median_ns' column: {counter} total files: {len(files)}")
    return max(occupancies)
    
def get_listdir(path):
    return os.listdir(path)

def arrays_to_dict(keys_arr,values_arr):
    return dict(zip(keys_arr, values_arr))
    

def write_dat(dict,output):
    with open(output, 'w') as f:
        f.write('statistic\tvalue\n')
        for key, value in dict.items():
            f.write(f'{key}\t{value}\n')
    print(f"\nWrote statistics to: {output}")

def main():
    args = parse_arguments()
    # Get list of csv files
    files = [os.path.join(args.lifetimes_dir, f) for f in get_listdir(args.lifetimes_dir) if f.endswith('.csv')]

    max_occupancy = get_max_occupancy(files)

    stats_dict = {
        'max_occupancy': max_occupancy
    }
    write_dat(stats_dict, args.output)

if __name__ == "__main__":
    main()



    
    
    
