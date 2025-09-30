"""
Read the cluster.log files inside the /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering/{case}/clust_c075/cluster.log

The directories inside the /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering/ contain different cases.

Searches for the line cl. | #st  rmsd | middle rmsd | cluster members
  1 | 1951  0.467 | 939000 .377 |      0   1000   2000   3000   4000   5000   6000

Where the fourth column is the center of the cluster at a given timepoint.

Save the frame using the command:
dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all
gmx trjconv -f {dir}/{case}/trajectory.xtc -s {dir}/{case}/topology.tpr -o {dir}/centers/{case}.pdb -b {timepoint}
"""