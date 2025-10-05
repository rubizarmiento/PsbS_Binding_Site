"""
Reads a csv file and copies the representative pdbs and tpr files for the type of binding event

With the columns:
tag original     new          old_chains count tag_number
n_s sim_4_A4_N_S sim_4_A4_n_s N_S        6     1

The gets the unique tags and writes a csv file with the columns:
{tag_number}_{tag} {tag}  

"""

import os
import sys
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd


csv_path=sys.argv[1]
out_csv=sys.argv[2] if len(sys.argv)>2 else "unique_basenames.csv"

print(f"Processing CSV file: {csv_path}")
# Read the CSV
df = pd.read_csv(csv_path,header=0,sep=' ')

# Get unique tags
unique_tags = df['tag'].unique()
unique_tag_numbers = df.drop_duplicates(subset=['tag'], keep='first')['tag_number'].values
unique_basenames = [f"{num}_{tag}" for num, tag in zip(unique_tag_numbers, unique_tags)]

# Create a new DataFrame
df_unique = pd.DataFrame({
    'unique_basename': unique_basenames,
    'tag': unique_tags
})

# Write to CSV
df_unique.to_csv(out_csv, index=False, sep=' ')
print(f"Wrote unique basenames to {out_csv}")
