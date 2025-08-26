

"""
Author: Rubi Zarmiento-Garcia
Date: 25/08/2025

Run contact analysis and saves the results in a boolean matrix.

Arguments
- gro: Path to the GRO file.
- xtc: Path to the XTC file.
- sel1: Selection of atoms for the first group.
- sel2: Selection of atoms for the second group.
- o: Output csv file path.
- cutoff: Distance cutoff for contacts.
- group_by1: Attribute to group the first selection by.
- group_by2: Attribute to group the second selection by.

Returns:
    Saves contact matrix in a csv file

Example 
python3 contact_analysis.py -f topol.gro -traj traj.xtc -sel1 "name CA" -sel2 "name CB" -o contact_matrix.csv
"""

from lib_contacts import *

def parser_args():
    parser = argparse.ArgumentParser(description="Run contact analysis.")
    parser.add_argument("-f", type=str, required=True, help="Path to the GRO file.")
    parser.add_argument("-traj", type=str, required=True, help="Path to the XTC file.")
    parser.add_argument("-sel1", type=str, required=True, help="Selection of atoms for the first group.")
    parser.add_argument("-sel2", type=str, required=True, help="Selection of atoms for the second group.")
    parser.add_argument("-o", type=str, required=True, help="Output csv file path.")
    parser.add_argument("-cutoff", type=float, default=0.35, help="Distance cutoff for contacts.")
    parser.add_argument("-group_by1", type=str, default="resids", help="Attribute to group the first selection by.")
    parser.add_argument("-group_by2", type=str, default="resids", help="Attribute to group the second selection by.")
    return parser.parse_args()

def main():
    args = parser_args()
    gro = args.f
    xtc = args.traj
    sel1 = args.sel1
    sel2 = args.sel2
    cutoff = args.cutoff
    group_by1 = args.group_by1
    group_by2 = args.group_by2
    output_file = args.o

    u = mda.Universe(gro, xtc)
    sel1 = u.select_atoms(sel1)
    sel2 = u.select_atoms(sel2)
    n_frames = len(u.trajectory)
    print("Number of residues in selection 1: ", sel1.n_residues)
    print("Number of residues in selection 2: ", sel2.n_residues)
    print("Number of frames: ", n_frames)

    contact_matrix_obj = ContactMatrix(u, sel1, sel2, cutoff, group_by1=group_by1, group_by2=group_by2)
    contact_matrix_list = contact_matrix_obj.calculate_contact_matrix_per_observation(n_frames=n_frames-1)

    # Print statistics
    for i, matrix in enumerate(contact_matrix_list):
        print(f"Contact matrix for frame {i}:")
        print(matrix.describe())

    #Save avg_matrix and std_matrix to csv
    contact_matrix_list[0].to_csv(output_file)
    print("Contact matrix saved to: ", output_file)

#Run main
if __name__ == "__main__":
    main()