

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
- dt: Time step between frames e.g. 1 ns.
- min_event_ns: Minimum event duration in nanoseconds. Events shorter than this will be discarded. Default 0 ns.
- n_frames: Maximum number of frames to process (for testing). If not specified, processes all frames.
Returns:
    Saves contact matrix in a csv file

Example 
python3 contact_analysis.py -f topol.gro -traj traj.xtc -sel1 "name CA" -sel2 "name CB" -o contact_matrix.csv
python3 contact_analysis.py -f topol.gro -traj traj.xtc -sel1 "name CA" -sel2 "name CB" -o contact_matrix.csv -n_frames 100
"""

from lib_contacts import *

def parser_args():
    parser = argparse.ArgumentParser(description="Run contact analysis.")
    parser.add_argument("-f", type=str, required=True, help="Path to the GRO file.")
    parser.add_argument("-traj", type=str, required=True, help="Path to the XTC file.")
    parser.add_argument("-sel1", type=str, required=True, help="Selection of atoms for the first group.")
    parser.add_argument("-sel2", type=str, required=True, help="Selection of atoms for the second group.")
    parser.add_argument("-o", type=str, required=True, help="Output csv files path.")
    parser.add_argument("-prefix", type=str, required=True, help="Preffix for output files.")
    parser.add_argument("-cutoff", type=float, default=8, help="Distance cutoff for contacts.")
    parser.add_argument("-group_by1", type=str, default="resids", help="Attribute to group the first selection by.")
    parser.add_argument("-group_by2", type=str, default="resids", help="Attribute to group the second selection by.")
    parser.add_argument("-dt", type=float, default=1, help="Time step between frames in ns.")
    parser.add_argument("-min_event_ns", type=float, default=0, help="Minimum event duration in ns.")
    parser.add_argument("-n_frames", type=int, default=None, help="Maximum number of frames to process (for testing). If not specified, processes all frames.")
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
    odir = args.o
    dt = args.dt
    min_event_ns = args.min_event_ns
    max_frames = args.n_frames

    u = mda.Universe(gro, xtc)
    
    # Limit trajectory to first n_frames if specified
    if max_frames is not None and max_frames < len(u.trajectory):
        # Transfer only the first max_frames to memory
        u.transfer_to_memory(start=0, stop=max_frames-1, step=1)
        print(f"Limited trajectory to first {max_frames} frames")
    
    sel1 = u.select_atoms(sel1)
    sel2 = u.select_atoms(sel2)
    
    # Determine number of frames to process
    total_frames = len(u.trajectory)
    n_frames = total_frames
    print(f"Processing {n_frames} frames")

    contact_matrix_obj = ContactMatrix(u, sel1, sel2, cutoff, group_by1=group_by1, group_by2=group_by2)
    contact_matrix_list = contact_matrix_obj.calculate_contact_pairs_matrix_per_observation(n_frames=n_frames-1)

    events_df, residue_summary_df, protein_summary = compute_lifetimes_from_contacts(
    contact_matrix_list[0], dt,min_event_ns)
    
    #Save dataframes 
    #Sort events_df by lifetime
    #events_df = events_df.sort_values(by='lifetime_ns', ascending=False)
    events_df.to_csv(f"{odir}/{args.prefix}_events_df.csv", index=False,float_format='%.2f')
    #Sort summary by resid - convert to numeric for proper sorting
    if 'resid' in residue_summary_df.columns:
        residue_summary_df['resid_numeric'] = pd.to_numeric(residue_summary_df['resid'], errors='coerce')
        residue_summary_df = residue_summary_df.sort_values(by='resid_numeric').drop('resid_numeric', axis=1).reset_index(drop=True)

    residue_summary_df.to_csv(f"{odir}/{args.prefix}_residue_summary_df.csv", index=False,float_format='%.2f')
    print(protein_summary)
    #protein_summary.to_csv(f"{odir}/{args.prefix}_protein_summary.csv", index=False)

#Run main
if __name__ == "__main__":
    main()