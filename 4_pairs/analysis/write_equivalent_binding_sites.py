"""
Workflow:
    -Read /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/equivalent_chains.csv into a dictionary
    -Read the .xtc files from /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
    -The chainIDs are found after "A1","A2", "A3","A4". Make a dictionary with the basename changing the chains to the keys in the dictionary
"""\

import os

os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd
from pandas.core import strings



# GEt equivalent chains
csv_file = '/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/equivalent_chains.csv'
df = pd.read_csv(csv_file, sep=' ', header=None, dtype=str)
key_column = df.columns[0]  # Assuming the first column contains the keys
data_dict = df.set_index(key_column).T.to_dict('list')
# Drop NaN values from the lists in the dictionary
for key in data_dict:
    data_dict[key] = [x for x in data_dict[key] if pd.notna(x)]
# Get basenames
dir="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj"
files = [f for f in os.listdir(dir) if f.endswith('.xtc')]
files = [os.path.join(dir, f) for f in files]
basenames = [os.path.basename(f) for f in files]
basenames = [s.replace(".xtc", "") for s in basenames]

# Split
parts = [basename.split('_') for basename in basenames]

chains_arr_of_arr = []
for part in parts:
    for i in ["A1", "A2", "A3", "A4"]:
        if i in part:
            index_i = part.index(i)
            chains_part = part[index_i + 1:]  # Get all parts after A1, A2, A3, or A4
            chains_arr_of_arr.append(chains_part)
            break  # Exit the inner loop once a match is found

# Invert the dictionary: map equivalents back to original keys
inverted_dict = {}
for key, equivs in data_dict.items():
    for equiv in equivs:
        inverted_dict[equiv] = key

# Replace chains with original keys from the inverted dictionary
new_chains_arr_of_arr = []
for chains in chains_arr_of_arr:
    new_chains = []
    for chain in chains:
        if chain in inverted_dict:
            new_chains.append(inverted_dict[chain])  # Replace with the original key
        else:
            new_chains.append(chain)  # Keep the original chain if not found in the inverted dictionary
    new_chains_arr_of_arr.append(new_chains)

# Modify the basenames to reflect the changes
new_basenames = []
for i, part in enumerate(parts):
    for j in ["A1", "A2", "A3", "A4"]:
        if j in part:
            index_j = part.index(j)
            new_part = part[:index_j + 1] + new_chains_arr_of_arr[i]  # Replace chains with new chains
            new_basename = '_'.join(new_part)
            new_basenames.append(new_basename)
            break  # Exit the inner loop once a match is found

#Basenames after A1, A2, A3 or A4 for new_basenames
new_parts = [basename.split('_') for basename in new_basenames]
new_basenames_no_A = []
for part in new_parts:
    for i in ["A1", "A2", "A3", "A4"]:
        if i in part:
            index_i = part.index(i)
            chains_part = part[index_i + 1:]  # Get all parts after A1, A2, A3, or A4
            new_basenames_no_A.append('_'.join(chains_part))
            break  # Exit the inner loop once a match is found

# Create a DataFrame with the old and new basenames
chains_arr_of_arr = ['_'.join(chains) for chains in chains_arr_of_arr]

df = pd.DataFrame({
    'tag': new_basenames_no_A,
    'original': basenames,
    'new': new_basenames,
    'old_chains': chains_arr_of_arr
})

# Add column with the number of times a tag appears
df['count'] = df.groupby('tag')['tag'].transform('count')

# Get unique tags sorted by count descending
tag_count = df.groupby('tag')['count'].first().sort_values(ascending=False)
unique_tags = tag_count.index.tolist()

# Add tag number based on sorted order (most frequent first)
tag_numbers = {tag: i + 1 for i, tag in enumerate(unique_tags)}
df['tag_number'] = df['tag'].map(tag_numbers)

# Sort by tag_number
df = df.sort_values(by=['tag_number'])


# Write the DataFrame to a CSV file
output_file = '/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv'
df.to_csv(output_file, sep=' ', index=False, header=True)

print(f"Basenames with equivalent chains written to {output_file}")