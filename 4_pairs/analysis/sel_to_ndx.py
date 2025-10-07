#Read slection from MDAnalysis and creates and index file for GROMACS.

'''
Arguments:
-f: The strcuture file supported by MDAnalysis, eg. gro, pdb, etc.
-sel: The selections separated by ''
-name: name of the selection, eg. 'Protein'. By defaul the name is sel1, sel2, etc.
-o : The name of the output file, by default index.ndx

Example:

python3 sel_to_ndx.py -f named.gro -sel 'index 0-1000' 'index 1001-3500 and name O31' -name 'Protein1' 'Protein2' -o index.ndx
'''

#Import libraries
import MDAnalysis as mda
import argparse
import os
import warnings

#Ignore warnings
warnings.filterwarnings("ignore")
#Read arguments
parser = argparse.ArgumentParser(description='Reads slection from MDAnalysis and creates and index file for GROMACS. nExample: python3 sel_to_ndx.py -f named.gro -sel "index 0-1000" "index 1001-3500 and name O31" -name "Protein1" "Protein2" -o index.ndx')
parser.add_argument('-f', type=str, help='The strcuture file supported by MDAnalysis, eg. gro, pdb, etc.')
parser.add_argument('-sel', type=str, nargs='+', help='The selections separated by "".')
parser.add_argument('-name', type=str, nargs='+', help='name of the selection, eg. "Protein". By defaul the name is sel1, sel2, etc.')
parser.add_argument('-o', type=str, help='The name of the output file, by default index.ndx', default='index.ndx')


#Fuction to read the arguments
def read_args():
    args = parser.parse_args()
    return args

#Function to check the arguments
def check_args(args):
    #Check if the file exists
    if not os.path.isfile(args.f):
        raise FileNotFoundError(f"The file {args.f} does not exist")
        exit()
    

    #Check ig the names are given, if not, create names for the selections sel1, sel2, etc.
    if args.name is None:
        print("The names of the selections were not given with the -name argument")
        exit()

    #Check if the selections and names are lists, if not, convert them to lists separated by "" or ' '
    if not isinstance(args.sel, list):
        args.sel = [s for sel in args.sel for s in sel.split('""' if '""' in sel else "' '")]
    if not isinstance(args.name, list):
        args.name = [n for name in args.name for n in name.split('""' if '""' in name else "' '")]

    #Check if the number of selections is equal to the number of names
    if len(args.sel) != len(args.name):
        print("The number of selections is not equal to the number of names")
        print("The number of selections is {}".format(len(args.sel)))
        print("The number of names is {}".format(len(args.name)))
        exit()


    return args

def print_info(step=""):
    print("\n----" + step + "----")
    return None


#Function to create the mdanalysis universe, and save the selections as index files
def create_ndx(args):
    print_info("READING UNIVERSE")
    #Create the universe
    universe = mda.Universe(args.f)

    print_info("CREATING INDEX FILE")
    #Create the selections, check if the selection is empty
    with mda.selections.gromacs.SelectionWriter(args.o, mode='w') as ndx:
        for selection, name in zip(args.sel, args.name):
            sel = universe.select_atoms(selection)
            atoms = sel.atoms
            if len(atoms) == 0:
                print("The selection {} is empty".format(selection))
                exit()
            ndx.write(sel, name=name)
            #Print the number of atoms in the selection
            print("The selection {} has {} atoms".format(selection, len(atoms)))
    print("The index file " + args.o + " was created")

#Main function
def main():
    #Read the arguments
    args = read_args()
    #Check the arguments
    args = check_args(args)
    #Create the universe
    create_ndx(args)

#Run the main function
if __name__ == "__main__":
    main()


