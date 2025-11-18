"""
Assigns the resid and chain of a atomistic PDB to a coarse-grained PDB

Arguments:
    - ref: Input atomistic PDB file
    - i: Input coarse-grained PDB file
    - o: Output coarse-grained PDB file

Example:
    python assign_resid_chain_from_pdb.py -ref atomistic.pdb -i coarse-grained.pdb -o output.pdb

Workflow:
1. Read atomistic PDB to get unique resid, resname, chainID combinations
2. Read coarse-grained PDB and assign resid and chainID based on atomistic info
    
"""
import MDAnalysis as mda
import numpy as np
import os
import sys
import argparse
import warnings
#Ignore warnings
warnings.filterwarnings("ignore", category=UserWarning)
from collections import OrderedDict

def read_arguments():
    parser = argparse.ArgumentParser(description='Assigns the resid and chain of a atomistic PDB to a coarse-grained PDB')
    parser.add_argument("-ref", "--reference", type=str, help="Input atomistic PDB file")
    parser.add_argument("-i", "--input", type=str, help="Input coarse-grained PDB file")
    parser.add_argument("-o", "--output", type=str, help="Output coarse-grained PDB file")
    args = parser.parse_args()
    return args

def get_unique_resname_resid_chain(u):
    resid = u.residues.residues.resids
    resname = u.residues.residues.resnames
    chainID = u.residues.residues.chainIDs
    #Join each element in the array
    resname_resid_chain = []
    for i in range(len(resid)):
        #Make an array of chainID_resid_resname
        resname_resid_chain.append(str(resid[i]) + '_' + resname[i]  + '_' + chainID[i])

    #Make flat array
    resname_resid_chain = [item for sublist in resname_resid_chain for item in sublist]
    #Get unique elements in numpy array
    unique_resname_resid_chain = list(OrderedDict.fromkeys(resname_resid_chain))

    #Get resnames from unique_resname_resid_chain
    resnames = [x.split('_')[1] for x in unique_resname_resid_chain]
    #Get unique elements in numpy array
    unique_resnames = np.unique(resnames)
    n_resname = []
    counter = 0
    first_value = resnames[0]
    for i in range(len(resnames)):
        current_value = resnames[i]
        if current_value == first_value:
            counter += 1
        else:
            n_resname.append(counter)
            counter = 1
            first_value = current_value
    #Append the last counter
    n_resname.append(counter)

    return unique_resname_resid_chain

#Reads a atomistic PDB and assigns the resid and chains to a coarse-grained PDB
def assign_resid_chain(pdb_aa, pdb_cg):
    """
    Assigns the resid and chain of a atomistic PDB to a coarse-grained PDB
    """
    u_aa = mda.Universe(pdb_aa)
    u_cg = mda.Universe(pdb_cg)
    unique_resname_resid_chain = []
    unique_resname_resid_chain = get_unique_resname_resid_chain(u_aa)

    #Add dummy chainID
    u_cg.residues.atoms.chainIDs = ['A'] * u_cg.atoms.n_atoms

    n_atoms = []
    counter = 0
    first_value = unique_resname_resid_chain[0]
    for i in range(u_cg.atoms.n_atoms):
        current_value = str(u_cg.atoms.resids[i]) + '_' + str(u_cg.atoms.chainIDs[i])
        if current_value == first_value:
            counter += 1
        else:
            n_atoms.append(counter)
            counter = 1
            first_value = current_value

    n_atoms.append(counter)


    #Duplicate each unique chainID_resid by n_atoms[i] times
    chainID_cg = []
    resid_cg = []
    for i in range(len(n_atoms)):
        chainID_cg.extend([unique_resname_resid_chain[i].split('_')[2]] * n_atoms[i])
        resid_cg.extend([unique_resname_resid_chain[i].split('_')[1]] * n_atoms[i])

    chainID_cg = [x.split('_')[0] for x in chainID_cg]
    resid_cg = [x.split('_')[1] for x in resid_cg]
    #Make resid a flat array
    resid_cg = [item for sublist in resid_cg for item in sublist]
    

    print('Number of residues: ' + str(len(resid_cg)))
    print('Number of atoms: ' + str(len(chainID_cg)))
    #u_cg.residues.atoms.resids = resid_cg
    u_cg.residues.atoms.chainIDs = chainID_cg
    u_cg.residues.residues.resids = resid_cg


    return u_cg

def add_resid_chainID(u, resid_resnames_chainID):
    """
    Add resid and chainID to a PDB file
    """

    #Get atoms per residue
    resnames = u.atoms.resnames
    resids = u.atoms.resids
    resids_resname = [ str(resids[i]) + '_' + str(resnames[i]) for i in range(len(resids))]
    counts = []
    start = 0
    for i in range(1, len(resids_resname)):
        if resids_resname[i] != resids_resname[i-1]:
            counts.append(i - start)
            start = i
    counts.append(len(resids_resname) - start)  # for the last group

    arr_resid_resnames_chainID = []
    for i in range(len(counts)):
        arr = [resid_resnames_chainID[i]] * counts[i]
        arr_resid_resnames_chainID.append(arr)
    arr_resid_resnames_chainID = [item for sublist in arr_resid_resnames_chainID for item in sublist]

    if len(arr_resid_resnames_chainID) != u.atoms.n_atoms:
        print(f"Error: Length of arr_resid_resnames_chainID ({len(arr_resid_resnames_chainID)}) does not match number of atoms in PDB ({u.atoms.n_atoms})")
        exit()
        

    #Make flat array
    chainIDs = [x.split('_')[2] for x in arr_resid_resnames_chainID]
    resid = [x.split('_')[0] for x in resid_resnames_chainID]
    resname = [x.split('_')[1] for x in arr_resid_resnames_chainID]
    u.residues.resids = resid
    u.atoms.chainIDs = chainIDs

    return u


def main():
    args = read_arguments()
    pdb_aa = args.reference
    pdb_cg = args.input
    output = args.output

    #Check if files exist
    if not os.path.isfile(pdb_aa):
        print('Error: File ' + pdb_aa + ' does not exist')
        sys.exit(1)
    if not os.path.isfile(pdb_cg):
        print('Error: File ' + pdb_cg + ' does not exist')
        sys.exit(1)

    u_aa= mda.Universe(pdb_aa)
    u_cg = mda.Universe(pdb_cg)
    unique_resname_resid_chain = get_unique_resname_resid_chain(u_aa)

    u_new = add_resid_chainID(u_cg, unique_resname_resid_chain)
    #Write the new PDB file
    u_new.atoms.write(output)
    print('Output file: ' + output)
if __name__ == "__main__":
    main()