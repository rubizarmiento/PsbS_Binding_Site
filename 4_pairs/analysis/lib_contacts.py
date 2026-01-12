import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from MDAnalysis.analysis.distances import distance_array
import os
os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd

#import seaborn as sns
import sys
import time
import warnings
warnings.filterwarnings('ignore')
import re
from typing import Callable, Tuple, Dict
import argparse


class ContactMatrix:
    """
    Class to calculate the contact matrix between two selections of atoms within a molecular universe.
    Author: Rubi Zarmiento Garcia
    """
    def __init__(self, universe, selection1, selection2, cutoff,group_by1='resnames',group_by2='resnames'):
        """
        Initializes the ContactMatrix class.
        
        Parameters:
        universe (MDAnalysis.Universe): The molecular universe containing the trajectory.
        selection1 (MDAnalysis.AtomGroup): The first selection of atoms.
        selection2 (MDAnalysis.AtomGroup): The second selection of atoms.
        cutoff (float): The cutoff distance for contacts in angstroms.
        """
        self.universe = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.cutoff = cutoff
        # Get the attributes to group by, if it is a list, then sum the attributes
        try:
            self.attributes1 = selection1.atoms.__getattribute__(group_by1)
        except:
            if len(selection1.atoms) == 1:
                self.attributes1 = selection1.atoms.__getattribute__(group_by1[0])   
            else: 
                if isinstance(group_by1, list):
                    list_attributes1 = []
                    list_attributes2 = []
                    for attribute in group_by1:
                        list_1 = selection1.atoms.__getattribute__(attribute)
                        #list_2 = selection2.atoms.__getattribute__(attribute)
                        #convert to string
                        list_1 = list(map(str, list_1))
                        list_attributes1.append(list_1)
                        #list_2 = list(map(str, list_2))
                        #list_attributes2.append(list_2)
                    #Sum the first element of each list, then the second element of each list, etc.
                    self.attributes1 = ['%s%s' % x for x in zip(*list_attributes1)]
                    #self.attributes2 = ['%s%s' % x for x in zip(*list_attributes2)]
        try:
            self.attributes2 = selection2.atoms.__getattribute__(group_by2)
        except:
            if len(selection2.atoms) == 1:
                self.attributes2 = selection2.atoms.__getattribute__(group_by2[0])
            else:
                if isinstance(group_by2, list):
                    list_attributes2 = []
                    for attribute in group_by2:
                        list_2 = selection2.atoms.__getattribute__(attribute)
                        #convert to string
                        list_2 = list(map(str, list_2))
                        list_attributes2.append(list_2)
                    #Sum the first element of each list, then the second element of each list, etc.
                    self.attributes2 = ['%s%s' % x for x in zip(*list_attributes2)]

        self.contacts_matrix = np.zeros((len(self.attributes1), len(self.attributes2)), dtype=int)
        self.distances = np.empty((selection1.n_atoms, selection2.n_atoms), dtype=float)
        self.indices1 = self.attributes1
        self.indices2 = self.attributes2
        self.n_frames = len(universe.trajectory)
        
    def calculate_contacts(self):
        """
        Calculates the contact matrix over the entire trajectory and groups by indices1 and indices2.
        """
        for frame in self.universe.trajectory:
            print(f"{frame.frame}/{len(self.universe.trajectory)}", end="\r", file=sys.stderr)
            distance_array(self.selection1.atoms.positions, self.selection2.atoms.positions, result=self.distances)
            contacts = (self.distances < self.cutoff) & (self.distances > 0)
            # Initialize dataframe with unique resnames instead of resnames
            df = pd.DataFrame(contacts, index=self.indices1, columns=self.indices2)
            # Group by resid and if one element is true, then the whole resid is true
            df = df.groupby(df.columns, axis=1).any()
            df = df.groupby(df.index, axis=0).any()
            
            self.contacts_matrix += df.values
        #Set column and row names
        self.contacts_matrix = pd.DataFrame(self.contacts_matrix, index=np.unique(self.attributes1), columns=np.unique(self.selection2.attributes2))
        return self.contacts_matrix
    
    def calculate_contact_matrix_per_observation(self, n_frames):
        """
        Save a contact matrix for every n_frames and stop if the largest value of saves_at is reached.
        """
        contact_matrix_list = []
        saves_at = np.arange(n_frames, self.n_frames + 1, n_frames)
        largest_save_point = saves_at[-1]  # Get the largest value from saves_at

        # Initialize contact matrix outside the loop
        contacts_matrix = np.zeros((len(np.unique(self.attributes1)), len(np.unique(self.attributes2))), dtype=int)
        
        for frame in self.universe.trajectory:
            # Stop if the largest save point is reached
            if frame.frame > largest_save_point:
                print(f"Stopping at frame: {frame.frame}")
                break
            
            # Calculate distances and contacts for the current frame
            distance_array(self.selection1.atoms.positions, self.selection2.atoms.positions, result=self.distances)
            contacts = (self.distances < self.cutoff) & (self.distances > 0)
            df = pd.DataFrame(contacts, index=self.indices1, columns=self.indices2)
            df = df.groupby(df.columns, axis=1).any()
            df = df.groupby(df.index, axis=0).any()
            contacts_matrix += df.values
            
            # Check if the current frame is a point at which to save the contact matrix
            if frame.frame in saves_at:
                print(f"Saving matrix at frame: {frame.frame}")
                contact_matrix_df = pd.DataFrame(contacts_matrix, index=np.unique(self.attributes1), columns=np.unique(self.attributes2))
                contact_matrix_list.append(contact_matrix_df)
                
                # Reset contact matrix for the next period
                contacts_matrix = np.zeros((len(np.unique(self.attributes1)), len(np.unique(self.attributes2))), dtype=int)
        print("Number of contact matrices saved: ", len(contact_matrix_list))
        return contact_matrix_list
    
    def calculate_contact_pairs_matrix_per_observation(self, n_frames):
        """
        Save a contact matrix for every n_frames and stop if the largest value of saves_at is reached.
        """
        contact_matrix_list = []
        saves_at = np.arange(n_frames, self.n_frames + 1, n_frames)

        largest_save_point = saves_at[-1]  # Get the largest value from saves_at

        # Initialize contact matrix outside the loop
        contacts_matrix = np.zeros((len(np.unique(self.attributes1)), len(np.unique(self.attributes2))), dtype=int)
        
        #Initialize attribute1_attribute2 x frames matrix e.g. resname1_resname2 x frames
        contact_pairs_matrix = []
        # Merge attributes1 and attributes2 with the format 'attribute1_attribute2 and the shape is (n_attributes1 * n_attributes2
        names = []
        for attr1 in np.unique(self.attributes1):
            for attr2 in np.unique(self.attributes2):
                names.append(f"{attr1}_{attr2}")

        print("Shape of contact pairs matrix: ", len(names), self.n_frames)
        chunk_start_frame = 0  # Track the starting frame of the current chunk
        for frame in self.universe.trajectory:
            # Stop if the largest save point is reached
            if frame.frame > largest_save_point:
                print(f"Stopping at frame: {frame.frame}")
                break
            
            # Calculate distances and contacts for the current frame
            distance_array(self.selection1.atoms.positions, self.selection2.atoms.positions, result=self.distances)
            contacts = (self.distances < self.cutoff) & (self.distances > 0)
            df = pd.DataFrame(contacts, index=self.indices1, columns=self.indices2)
            df = df.groupby(df.columns, axis=1).any()
            df = df.groupby(df.index, axis=0).any()
            
            # Flatten the dataframe and convert to int (0/1 instead of False/True)
            contact_pairs = df.values.flatten().astype(np.uint8)
            # Update the contact pairs matrix for the current frame
            contact_pairs_matrix.append(contact_pairs)

            # Check if the current frame is a point at which to save the contact matrix
            if frame.frame in saves_at:
                print(f"Saving matrix at frame: {frame.frame}")
                # Convert to numpy array first for better performance, then to DataFrame
                contact_matrix_array = np.array(contact_pairs_matrix, dtype=np.uint8)
                # Create frame indices starting from chunk_start_frame, add 1 to skip frame 0 (structure file)
                frame_indices = np.arange(chunk_start_frame + 1, chunk_start_frame + 1 + len(contact_pairs_matrix))
                contact_matrix_df = pd.DataFrame(contact_matrix_array, index=frame_indices, columns=names)
                contact_matrix_list.append(contact_matrix_df)
                # Reset contact matrix for the next period
                contact_pairs_matrix = []
                chunk_start_frame = frame.frame + 1  # Next chunk starts after this frame
        print("Number of contact matrices saved: ", len(contact_matrix_list))
        
        return contact_matrix_list
    
    def print_contact_matrix_list(self, contact_matrix_list):
        """
        Prints the contact matrix list.
        
        Parameters:
        contact_matrix_list (list): The list of contact matrices.
        """
        for i, contact_matrix in enumerate(contact_matrix_list):
            print(f"Contact matrix {i+1}:")
            print(contact_matrix)
            print()

    def avg_std_matrix_list(self, contact_matrix_list):
        """
        Calculates the average and standard deviation of the contact matrix list.
        
        Parameters:
        contact_matrix_list (list): The list of contact matrices.
        """
        avg_matrix = np.mean(contact_matrix_list, axis=0)
        std_matrix = np.std(contact_matrix_list, axis=0)
        #Add column and row names
        avg_matrix = pd.DataFrame(avg_matrix, index=np.unique(self.attributes1), columns=np.unique(self.attributes2))
        std_matrix = pd.DataFrame(std_matrix, index=np.unique(self.attributes1), columns=np.unique(self.attributes2))
        return avg_matrix, std_matrix
    
    def time_matrix_list(self, contact_matrix_list,universe):
        """
        Calculates the average and standard deviation of the contact matrix list.
        
        Parameters:
        contact_matrix_list (list): The list of contact matrices.
        universe (MDAnalysis.Universe): The molecular universe containing the trajectory.
        """
        total_time = universe.trajectory.totaltime
        total_frames = universe.trajectory.n_frames
        time_per_frame = total_time/total_frames
        time_matrix_list = []
        for matrix in contact_matrix_list:
            time_matrix = matrix*time_per_frame
            time_matrix_list.append(time_matrix)
        return time_matrix_list
            
    def distance_vs_time(self):
        """
        Calculates distance vs time between two selections.
        """
        #Distance matrix shape: n_frames x n_atoms1 
        distance_matrix = np.empty((self.n_frames, len(np.unique(self.attributes1))), dtype=float)
        print("Shape of distance matrix: ", distance_matrix.shape)
        #distance_array(self.selection1.atoms.positions, self.selection2.atoms.positions, result=self.distances)
        for frame in self.universe.trajectory:
            print(f"{frame.frame}/{len(self.universe.trajectory)}", end="\r", file=sys.stderr)
            distance_array(self.selection1.atoms.positions, self.selection2.atoms.positions, result=self.distances)
            # Initialize dataframe with unique resnames instead of resnames
            df = pd.DataFrame(self.distances, index=self.indices1, columns=self.indices2)
            # Group by resid and only conserve the minimum distance of selection2
            df = df.groupby(df.index, axis=0).min()
            df = df.groupby(df.columns, axis=1).min()
            #df = df.groupby(df.columns, axis=1).any()
            #df = df.groupby(df.index, axis=0).any()
            #Save in distance matrix
            distance_matrix[frame.frame,:] = df.values.flatten()
        #Convert to dataframe
        distance_matrix = pd.DataFrame(distance_matrix, index=np.arange(1,self.n_frames+1), columns=np.unique(self.attributes1))
        

        return distance_matrix
        
#Function to plot the contact matrix
def plot_contact_matrix(contact_matrix,n_frames,figsize=(5,5),cmap="coolwarm"):
    """
    Plots the contact matrix.

    Parameters:
    contact_matrix (pd.DataFrame): The contact matrix.
    n_frames (int): The number of frames in the trajectory.
    figsize (tuple): The figure size.
    cmap (str): The colormap.
    """
    plt.figure(figsize=figsize)
    #Set the axis to start from lower to higher

    sns.heatmap(contact_matrix, cmap=cmap, vmin=0, vmax=n_frames)
    plt.axes().invert_yaxis()
    plt.show()

def plot_one_subplot_per_column_dataframe(df,ylabel="ylabel",xlabel="xlabel",ymax=None,pallette='tab10',hline=None):
    """
    Plots one lineplot per column of a dataframe in a single figure with subplots.

    Parameters:
    df (pd.DataFrame): The dataframe.
    ylabel (str): The ylabel.
    xlabel (str): The xlabel.
    ymax (float): The maximum value of the y axis.
    """
    n_subplots = df.shape[1]
    subplot_per_row = 3
    subplot_per_column = int(np.ceil(n_subplots/subplot_per_row))
    fig, axes = plt.subplots(subplot_per_column, subplot_per_row, figsize=(15, 5*subplot_per_column))
    for i, column in enumerate(df.columns):
        sns.lineplot(x=df.index, y=column, data=df, ax=axes[i//subplot_per_row, i%subplot_per_row], palette=pallette)
        #Set labels
        axes[i//subplot_per_row, i%subplot_per_row].set_title(column)
        axes[i//subplot_per_row, i%subplot_per_row].set_ylabel(ylabel)
        axes[i//subplot_per_row, i%subplot_per_row].set_xlabel(xlabel)
        #set y limit
        if ymax is not None:
            axes[i//subplot_per_row, i%subplot_per_row].set_ylim(0,ymax)
        if hline is not None:
            axes[i//subplot_per_row, i%subplot_per_row].axhline(y=hline, color="black", linestyle='--')

    plt.tight_layout()
    #Turn off any unused subplots
    for i in range(n_subplots,subplot_per_row*subplot_per_column):
        axes[i//subplot_per_row, i%subplot_per_row].axis('off')

    
    plt.show()


def extract_lifetimes(df):
    """
    Extract lifetimes from a DataFrame where each column contains binary time series data.
    
    Parameters:
    df (pd.DataFrame): DataFrame with shape (labels x time) containing binary data (0s and 1s)
    
    Returns:
    list: List of arrays, where each array contains the lifetimes (consecutive 1s) for each column
    
    Example:
    If a column has data [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1], 
    it returns [3, 2] representing two lifetimes of length 3 and 2.
    """
    lifetimes_per_column = []
    
    for column in df.columns:
        column_data = df[column].values
        lifetimes = []
        current_lifetime = 0
        
        for value in column_data:
            if value == 1:
                current_lifetime += 1
            else:
                if current_lifetime > 0:
                    lifetimes.append(current_lifetime)
                    current_lifetime = 0
        
        # Don't forget the last lifetime if the sequence ends with 1s
        if current_lifetime > 0:
            lifetimes.append(current_lifetime)
        
        lifetimes_per_column.append(np.array(lifetimes))
    
    return lifetimes_per_column

# Alternative more compact version using groupby
def extract_lifetimes_vectorized(df):
    """
    Vectorized version using pandas groupby for better performance.
    """
    lifetimes_per_column = []
    
    for column in df.columns:
        series = df[column]
        # Create groups where consecutive identical values are grouped together
        groups = (series != series.shift()).cumsum()
        
        # Filter for groups where the value is 1 and get their sizes
        lifetimes = series.groupby(groups).size()[series.groupby(groups).first() == 1].values
        lifetimes_per_column.append(lifetimes)
    
    return lifetimes_per_column

def assign_bfactor_to_universe(universe, resids, bfactor, decimal_places=0):
    """
    Assigns B-factors to the atoms in a universe based on the provided residue IDs and B-factor values.

    Parameters:
    universe (MDAnalysis.Universe): The molecular universe containing the trajectory.
    resids (np.ndarray): Array of residue IDs.
    bfactor (np.ndarray): Array of B-factor values corresponding to the residues.

    Returns:
    MDAnalysis.Universe: The updated universe with assigned B-factors.
    """
    for resid, bf in zip(resids, bfactor):
        # Select atoms by residue ID
        atoms = universe.select_atoms(f'resid {resid}')

        # Raise warning and skip if len(atoms) == 0
        if len(atoms) == 0:
            print(f"Warning: No atoms found for resid {resid}. Skipping B-factor assignment for this resid.")
            continue
        # Assign B-factor
        if decimal_places == 0:
            # Round to nearest integer and convert to int
            atoms.bfactors = int(round(bf))
        else:
            atoms.bfactors = np.round(bf, decimals=decimal_places)
    return universe

def save_universe_to_pdb(universe, filename):
    """
    Saves the universe to a PDB file with B-factors preserved.

    Parameters:
    universe (MDAnalysis.Universe or AtomGroup): The molecular universe or atom group to save.
    filename (str): The name of the output PDB file.
    """
    universe.write(filename)
    print(f"Universe saved to {filename}")


def compute_lifetimes_from_contacts(
    contacts_df: pd.DataFrame,
    dt: float,
    min_event_ns: float = 0.0,
    exclude_single_letters: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    """
    Compute lifetimes from a wide 'frames × resid_pairs' contact matrix.

    Parameters
    ----------
    contacts_df : pd.DataFrame
        Index: frames (int). Columns: residue-pair labels. Values: bool (or 0/1).
        Example columns: "101_315", "A_B", "LYS_ARG", etc.
    dt : float
        Time between saved frames in nanoseconds.
    min_event_ns : float, default 0.0
        Minimum lifetime (ns) to keep an event; set e.g. 5.0 to drop micro-events.
    exclude_single_letters : bool, default True
        If True, exclude single-letter residue IDs (likely chain IDs) from residue summary.

    Returns
    -------
    events_df : DataFrame
        Per-pair events with ['resid_i','resid_j','start_frame','end_frame','frames','lifetime_ns'].
    residue_summary_df : DataFrame
        Per-residue stats: ['resid','n_events','occupancy_pct','median_ns','p90_ns'].
    protein_summary : dict
        {'total_bound_ns': float, 'bound_fraction_pct': float}
    """
    # Ensure index are integer frames and values are boolean
    df = contacts_df.copy()
    # Convert any numeric to boolean; treat >0 as True
    if not np.issubdtype(df.dtypes.values[0], np.bool_):
        df = df.astype(bool)

    frame_min, frame_max = int(df.index.min()), int(df.index.max())
    n_frames = frame_max - frame_min + 1
    if n_frames <= 0:
        return (
            pd.DataFrame(columns=['resid_i','resid_j','start_frame','end_frame','frames','lifetime_ns']),
            pd.DataFrame(columns=['resid','n_events','occupancy_pct','median_ns','p90_ns']),
            dict(total_bound_ns=0.0, bound_fraction_pct=0.0)
        )

    # Parse residue pairs from column names (handles both string and numeric residues)
    def _parse_pair(col: str) -> tuple:
        # Split by underscore or common separators
        parts = str(col).replace('-', '_').split('_')
        if len(parts) >= 2:
            return parts[0], parts[1]
        else:
            raise ValueError(f"Cannot parse residue pair from column name: {col}")

    # Helper function to check if a residue ID should be excluded
    def _should_exclude_resid(resid: str) -> bool:
        if exclude_single_letters:
            # Exclude single letters (likely chain IDs)
            return len(str(resid)) == 1 and str(resid).isalpha()
        return False

    # Helper: ON-segment extraction
    def _segments(on_bool: np.ndarray) -> np.ndarray:
        """Return array of (start_idx, end_idx) inclusive (indices relative to df.index order)."""
        x = on_bool.astype(int)
        diff = np.diff(x, prepend=0, append=0)
        starts = np.where(diff == 1)[0]
        ends   = np.where(diff == -1)[0] - 1
        return np.stack([starts, ends], axis=1) if starts.size else np.empty((0,2), dtype=int)

    # Protein-level ON timeline (any pair ON)
    protein_on = df.any(axis=1).to_numpy()
    total_on_frames = int(protein_on.sum())
    total_bound_ns = total_on_frames * dt 
    bound_fraction_pct = 100.0 * (total_on_frames / n_frames)
    protein_summary = dict(total_bound_ns=float(total_bound_ns),
                           bound_fraction_pct=float(bound_fraction_pct))

    # Per-pair events
    events_rows = []
    resid_on: Dict[str, np.ndarray] = {}

    for col in df.columns:
        ri, rj = _parse_pair(col)
        on_series = df[col].to_numpy()

        segs = _segments(on_series)
        if segs.size:
            # Convert segments to lifetimes
            lengths = segs[:,1] - segs[:,0] + 1
            lifetimes_ns = lengths * dt 
            # Optional: drop short events
            if min_event_ns > 0.0:
                keep = lifetimes_ns >= float(min_event_ns)
                segs = segs[keep]
                lifetimes_ns = lifetimes_ns[keep]

            for (s_idx, e_idx), L_ns in zip(segs, lifetimes_ns):
                start_frame = int(df.index[s_idx])
                end_frame   = int(df.index[e_idx])
                events_rows.append(dict(
                    resid_i=ri, resid_j=rj,
                    start_frame=start_frame,
                    end_frame=end_frame,
                    frames=int(e_idx - s_idx + 1),
                    lifetime_ns=float(L_ns)
                ))

        # Accumulate per-residue ON (any partner) - but exclude single letters
        if not _should_exclude_resid(ri):
            if ri not in resid_on:
                resid_on[ri] = np.zeros(n_frames, dtype=bool)
            resid_on[ri] |= on_series

    events_df = (pd.DataFrame(events_rows) if events_rows
                 else pd.DataFrame(columns=['resid_i','resid_j','start_frame','end_frame','frames','lifetime_ns']))

    # Per-residue stats (excluding single letters)
    res_rows = []
    for resid, on_vec in resid_on.items(): 
        occ_pct = 100.0 * (on_vec.sum() / n_frames)
        segs = _segments(on_vec)
        if segs.size:
            lengths = segs[:,1] - segs[:,0] + 1
            lifetimes_ns = lengths * dt
            if min_event_ns > 0.0:
                lifetimes_ns = lifetimes_ns[lifetimes_ns >= float(min_event_ns)]
        else:
            lifetimes_ns = np.array([], dtype=float)

        if lifetimes_ns.size == 0:
            med = p90 = sum_ns = 0.0
            n_ev = 0
        else:
            med = float(np.median(lifetimes_ns))
            p90 = float(np.percentile(lifetimes_ns, 90))
            sum_ns = float(np.sum(lifetimes_ns))
            n_ev = int(lifetimes_ns.size)

        res_rows.append(dict(
            resid=resid,
            n_events=n_ev,
            occupancy_pct=float(occ_pct),
            median_ns=med,
            sum_ns=sum_ns,
            p90_ns=p90
        ))

        

    residue_summary_df = pd.DataFrame(res_rows)
    #residue_summary_df.head()
    #.sort_values('resid').reset_index(drop=True)
    return events_df, residue_summary_df, protein_summary