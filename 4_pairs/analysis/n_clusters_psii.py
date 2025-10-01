"""
This script summarizes the number of clusters per cut-off

"""



import os
import MDAnalysis as mda
from sys import argv


base_dir = argv[1]


#base_dir = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering_grouped"
dirs = ["clust_c05", "clust_c075"]



arr_of_arr = []
for directory in dirs:
    files_in_dir = os.listdir(f"{base_dir}/{directory}")
    arr = []
    for file in files_in_dir:
        if file.endswith('.pdb') and not file.endswith('_align.log'):
            path = f"{base_dir}/{directory}/{file}"
            if os.path.isfile(path) and os.path.getsize(path) > 0:  # Check file is not empty
                try:
                    u = mda.Universe(path)
                    n_frames = u.trajectory.n_frames
                    arr.append(n_frames)
                except Exception as e:
                    print(f"Warning: Could not read {path}: {e}")
                    continue
    arr_of_arr.append(arr)

#Save arr_of_arr in file
with open(f"{base_dir}/n_clusters_psii.txt","w") as f:
    for i in range(len(arr_of_arr)):
        arr = arr_of_arr[i]
        dir = dirs[i]
        f.write(f"{dir}: {chr(9).join([str(x) for x in arr])}\n")

print(f"Wrote {base_dir}/n_clusters_psii.txt")



