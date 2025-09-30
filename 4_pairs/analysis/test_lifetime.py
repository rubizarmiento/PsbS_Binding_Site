import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
odir="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/chain_4"

test_lifetime=f"/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/csv_files/lifetime/lifetime_chain4_events_df.csv"
#Keys: resid_i,resid_j,start_frame,end_frame,frames,lifetime_ns

gro=f"{odir}/initial_fit.pdb"
#xtc=test_ultrashort.xtc
xtc=f"{odir}/aligned_5000ns.xtc"
chain="4"
sel2="chainID A B"
sel1=f"chainID {chain} and not resname CLA CLB CHL *HG* PLQ PL9 *GG* *SQ* *PG* LUT VIO XAT NEO NEX BCR"
cutoff=8
dt=2
min_event_ns=100


df = pd.read_csv(test_lifetime)
u = mda.Universe(gro,xtc)

f0 = df['start_frame'].values
ff = df['end_frame'].values

def check_distance_at_frame(u, frame_num, sel1, sel2, cutoff):
    """Check if distance is below cutoff at a specific frame"""
    try:
        u.trajectory[frame_num]  # Go directly to frame
        pos1 = u.select_atoms(sel1).positions
        pos2 = u.select_atoms(sel2).positions
        dists = distance_array(pos1, pos2)
        return dists.min() < cutoff
    except (IndexError, KeyError):
        print(f"Warning: Frame {frame_num} not found in trajectory")
        return False

passed_events = []
failed_events = []

print(f"Testing {len(f0)} events...")
print(f"Trajectory has {len(u.trajectory)} frames (0 to {len(u.trajectory)-1})")

for i in range(len(f0)):
    start_frame = f0[i]
    end_frame = ff[i]
    
    # Check if frames are within trajectory bounds
    if start_frame >= len(u.trajectory) or end_frame >= len(u.trajectory):
        print(f"Event {i}: Frames {start_frame}-{end_frame} out of bounds")
        failed_events.append(i)
        continue
    
    # Check distance at start and end frames
    start_passed = check_distance_at_frame(u, start_frame, sel1, sel2, cutoff)
    end_passed = check_distance_at_frame(u, end_frame, sel1, sel2, cutoff)
    
    # Only keep events where both start and end frames pass the distance test
    if start_passed and end_passed:
        passed_events.append(i)
    else:
        failed_events.append(i)
        print(f"Event {i}: start_passed={start_passed}, end_passed={end_passed}")

print(f"Passed events: {len(passed_events)}/{len(f0)}")
print(f"Failed events: {len(failed_events)}")
if len(passed_events) < 10:  # Show details if few events pass
    print("Passed event indices:", passed_events)
