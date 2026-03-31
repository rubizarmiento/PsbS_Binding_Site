"""
Reads a csv file and copies the representative pdbs and tpr files for the type of binding event

With the columns:
tag original     new          old_chains count tag_number binding_time_ns
n_s sim_4_A4_N_S sim_4_A4_n_s N_S        6     1          4951

The gets the unique tags and writes a csv file with the columns:
{tag_number}_{tag} {tag} n_events median_ms

"""

import os
import sys
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
import numpy as np


csv_path=sys.argv[1]
out_csv=sys.argv[2] if len(sys.argv)>2 else "unique_basenames.csv"

print(f"Processing CSV file: {csv_path}")
# Read the CSV
df = pd.read_csv(csv_path,header=0,sep=' ')

# Get unique tags
unique_tags = df['tag'].unique()
unique_tag_numbers = df.drop_duplicates(subset=['tag'], keep='first')['tag_number'].values
unique_basenames = [f"{num}_{tag}" for num, tag in zip(unique_tag_numbers, unique_tags)]

# Get tag_rename if it exists
tag_rename_list = []
if 'tag_rename' in df.columns:
    for tag in unique_tags:
        tag_rows = df[df['tag'] == tag]
        # Get first value of tag_rename for this tag (should be same for all rows)
        tag_rename = tag_rows['tag_rename'].iloc[0]
        tag_rename_list.append(tag_rename)

# Calculate n_events (count of each tag) and median_ms (median binding time for each tag)
n_events_list = []
median_ms_list = []
avg_ms_list = []
std_ms_list = []
sum_ms_list = []

for tag in unique_tags:
    # Count occurrences of this tag
    tag_rows = df[df['tag'] == tag]
    n_events = len(tag_rows)
    n_events_list.append(n_events)
    
    # Calculate median of binding_time_ns and convert to microseconds
    binding_times = tag_rows['binding_time_ns'].values / 1000  # Convert ns to microseconds
    median_ms = np.median(binding_times)
    avg_ms = np.mean(binding_times)
    std_ms = np.std(binding_times)
    sum_ms = np.sum(binding_times)
    median_ms_list.append(median_ms)
    avg_ms_list.append(avg_ms)
    std_ms_list.append(std_ms)
    sum_ms_list.append(sum_ms)

# Create a new DataFrame
if tag_rename_list:
    df_unique = pd.DataFrame({
        'unique_basename': unique_basenames,
        'trajectory': unique_basenames,
        'tag': unique_tags,
        'tag_rename': tag_rename_list,
        'n_events': n_events_list,
        'median_ms': median_ms_list,
        'avg_ms': avg_ms_list,
        'std_ms': std_ms_list,
        'sum_ms': sum_ms_list   
    })
else:
    df_unique = pd.DataFrame({
        'unique_basename': unique_basenames,
        'trajectory': unique_basenames,
        'tag': unique_tags,
        'n_events': n_events_list,
        'median_ms': median_ms_list,
        'avg_ms': avg_ms_list,
        'std_ms': std_ms_list,
        'sum_ms': sum_ms_list   
    })

# Round the time columns to 2 decimal places
df_unique['median_ms'] = df_unique['median_ms'].round(2)
df_unique['avg_ms'] = df_unique['avg_ms'].round(2)
df_unique['std_ms'] = df_unique['std_ms'].round(2)
df_unique['sum_ms'] = df_unique['sum_ms'].round(2)

# Write to CSV
df_unique.to_csv(out_csv, index=False, sep=',')
print(f"Wrote unique basenames to {out_csv}")
print(f"Columns: {list(df_unique.columns)}")
print(f"\nFirst few rows:")
print(df_unique.head())
