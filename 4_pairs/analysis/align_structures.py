"""
Author: Rubi Zarmiento Garcia 
Aligns structures or trajectories using MDAnalysis.
This script aligns a mobile structure/trajectory to a reference structure or trajectory

Arguments:
    - mobile: Path to the mobile structure file (PDB format) or trajectory (XTC/TRR format).
    - mobiletop: Topology file for mobile trajectory (required if mobile is a trajectory).
    - ref: Path to the reference structure file (PDB/GRO format) or trajectory (XTC/TRR format).
    - reftop: Topology file for reference trajectory (required if ref is a trajectory).
    - sel: Selection string to specify which atoms to use for alignment (default is "name BB").
    - o: Output file path for the aligned structure/trajectory (default is "aligned.pdb").

Example usage:
    python align_structures.py -mobile /path/to/mobile.pdb -ref /path/to/reference.pdb -sel "name BB" -o aligned.pdb
    python align_structures.py -mobile /path/to/mobile.xtc -mobiletop /path/to/mobile_top.pdb -ref /path/to/reference.pdb -sel "name BB" -o aligned.xtc
    python align_structures.py -mobile /path/to/mobile.pdb -ref /path/to/ref.xtc -reftop /path/to/ref_top.pdb -sel "name BB" -o aligned.pdb
"""


import MDAnalysis as mda
from MDAnalysis.analysis import align
#Ignore warnings 
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')
def parser():
    import argparse
    parser = argparse.ArgumentParser(description="Aligns a mobile structure/trajectory to a reference structure/trajectory using MDAnalysis.")
    parser.add_argument("-mobile", required=True, help="Path to the mobile structure file (PDB format) or trajectory (XTC/TRR format).")
    parser.add_argument("-mobiletop", help="Topology file for mobile trajectory (required if mobile is a trajectory).")
    parser.add_argument("-ref", required=True, help="Path to the reference structure file (PDB/GRO format) or trajectory (XTC/TRR format).")
    parser.add_argument("-reftop", help="Topology file for reference trajectory (required if ref is a trajectory).")
    parser.add_argument("-sel", default="name BB", help="Selection string for alignment (default: 'name BB').")
    parser.add_argument("-o", default="aligned.pdb", help="Output file path for the aligned structure/trajectory (default: 'aligned.pdb').")

    return parser.parse_args()


def main():
    args = parser()
    
    # Load the mobile structure or trajectory
    if args.mobile.endswith(('.xtc', '.trr', '.dcd', '.nc')):
        # Mobile is a trajectory
        if not args.mobiletop:
            raise ValueError("Topology file (-mobiletop) is required when mobile is a trajectory")
        mobile = mda.Universe(args.mobiletop, args.mobile)
        is_mobile_trajectory = True
        print(f"Mobile is trajectory with {len(mobile.trajectory)} frames")
    else:
        # Mobile is a single structure file
        mobile = mda.Universe(args.mobile)
        is_mobile_trajectory = False
    
    # Load the reference structure or trajectory
    if args.ref.endswith(('.xtc', '.trr', '.dcd', '.nc')):
        # Reference is a trajectory
        if not args.reftop:
            raise ValueError("Topology file (-reftop) is required when reference is a trajectory")
        ref = mda.Universe(args.reftop, args.ref)
        # Use the first frame of the trajectory as reference
        ref.trajectory[0]  # Go to first frame
        print(f"Using first frame of trajectory {args.ref} as reference")
    else:
        # Reference is a single structure file
        ref = mda.Universe(args.ref)

    # Check if the selections have the same number of atoms
    mobile_sel = mobile.select_atoms(args.sel)
    ref_sel = ref.select_atoms(args.sel)
    
    print(f"Mobile selection: {len(mobile_sel)} atoms")
    print(f"Reference selection: {len(ref_sel)} atoms")
    
    if len(mobile_sel) != len(ref_sel):
        print(f"Mobile residues: {mobile_sel.n_residues}")
        print(f"Reference residues: {ref_sel.n_residues}")
        print(f"Mobile chains: {set(mobile_sel.chainIDs)}")
        print(f"Reference chains: {set(ref_sel.chainIDs)}")
        
        # Try to find common residues with same number of atoms
        mobile_resids = set(mobile_sel.resids)
        ref_resids = set(ref_sel.resids)
        common_resids = mobile_resids.intersection(ref_resids)
        
        if len(common_resids) > 0:
            print(f"Found {len(common_resids)} common residues")
            
            # Check each residue to find those with matching atom counts
            matching_resids = []
            for resid in sorted(common_resids):
                mobile_res_atoms = mobile.select_atoms(f"{args.sel} and resid {resid}")
                ref_res_atoms = ref.select_atoms(f"{args.sel} and resid {resid}")
                if len(mobile_res_atoms) == len(ref_res_atoms):
                    matching_resids.append(resid)
            
            if len(matching_resids) > 10:  # Need at least 10 residues for alignment
                print(f"Found {len(matching_resids)} residues with matching atom counts")
                # Create selection with only matching residues
                resid_str = " ".join(map(str, matching_resids))
                new_sel = f"{args.sel} and resid {resid_str}"
                print(f"Trying selection with matching residues")
                
                mobile_sel = mobile.select_atoms(new_sel)
                ref_sel = ref.select_atoms(new_sel)
                print(f"Final mobile selection: {len(mobile_sel)} atoms")
                print(f"Final reference selection: {len(ref_sel)} atoms")
                
                if len(mobile_sel) == len(ref_sel):
                    args.sel = new_sel  # Update the selection
                else:
                    raise ValueError(f"Still mismatched: {len(mobile_sel)} vs {len(ref_sel)}")
            else:
                raise ValueError(f"Not enough matching residues found: {len(matching_resids)}")
        else:
            raise ValueError(f"Selection '{args.sel}' results in different number of atoms: {len(mobile_sel)} (mobile) vs {len(ref_sel)} (ref). No common residues found.")
    
    # Perform the alignment
    if is_mobile_trajectory:
        # Mobile is a trajectory - align each frame
        print(f"Aligning {len(mobile.trajectory)} frames...")
        
        # Determine output format
        if args.o.endswith(('.xtc', '.trr', '.dcd', '.nc')):
            # Output is a trajectory
            with mda.Writer(args.o, mobile.atoms.n_atoms) as w:
                for ts in mobile.trajectory:
                    # Align current frame
                    mobile0 = mobile.select_atoms(args.sel).positions - mobile.select_atoms(args.sel).center_of_mass()
                    ref0 = ref.select_atoms(args.sel).positions - ref.select_atoms(args.sel).center_of_mass()
                    R, rmsd = align.rotation_matrix(mobile0, ref0)
                    
                    mobile.atoms.translate(-mobile.select_atoms(args.sel).center_of_mass())
                    mobile.atoms.rotate(R)
                    mobile.atoms.translate(ref.select_atoms(args.sel).center_of_mass())
                    
                    # Write aligned frame
                    mobile.trajectory.ts.time = ts.time
                    w.write(mobile.atoms)
                    
                print(f"RMSD for last frame: {rmsd:.3f} Angstroms")
        else:
            # Output is a single structure - align first frame only
            mobile.trajectory[0]  # Go to first frame
            mobile0 = mobile.select_atoms(args.sel).positions - mobile.select_atoms(args.sel).center_of_mass()
            ref0 = ref.select_atoms(args.sel).positions - ref.select_atoms(args.sel).center_of_mass()
            R, rmsd = align.rotation_matrix(mobile0, ref0)
            print(f"RMSD after alignment: {rmsd:.3f} Angstroms")
            mobile.atoms.translate(-mobile.select_atoms(args.sel).center_of_mass())
            mobile.atoms.rotate(R)
            mobile.atoms.translate(ref.select_atoms(args.sel).center_of_mass())
    else:
        # Mobile is a single structure
        mobile0 = mobile.select_atoms(args.sel).positions - mobile.select_atoms(args.sel).center_of_mass()
        ref0 = ref.select_atoms(args.sel).positions - ref.select_atoms(args.sel).center_of_mass()
        R, rmsd = align.rotation_matrix(mobile0, ref0)
        print(f"RMSD after alignment: {rmsd:.3f} Angstroms")
        mobile.atoms.translate(-mobile.select_atoms(args.sel).center_of_mass())
        mobile.atoms.rotate(R)
        mobile.atoms.translate(ref.select_atoms(args.sel).center_of_mass())

    # Write the aligned structure/trajectory to the output file
    if not (is_mobile_trajectory and args.o.endswith(('.xtc', '.trr', '.dcd', '.nc'))):
        # Write single structure (either mobile was single, or trajectory output was requested as single)
        mobile.atoms.write(args.o)
        print(f"Wrote aligned structure to {args.o}")
    else:
        print(f"Wrote aligned trajectory to {args.o}")

if __name__ == "__main__":
    main()

