"""
Workflow:
    1. Read {ref_pdb} to get a dictionary with the number of residues per chainID with MDAnalysis.
    2. Preprocess each PDB file to merge chains where residue numbering restarts at 1.
    3. For all the pdbs in {dir} get a dictionary with the number of residues per chain.
    4. The chains in the pdbs have a wrong name. To fix it:
        -First, the chainIDs are indicated in the filename after A1, A2, A3 or A4. 
            For example, in the file 14_sim_6_A4_7_8.pdb, the chainIDs are 7 and 8.
        -Then, using the first dictionary and comparing the number of resiudes, the correct chainID is assign.
"""

import MDAnalysis as mda
import os
import glob
import re
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

ref_pdb = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb"
dir = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering_grouped/clust_c075"

def get_residue_counts_from_pdb(pdb_file):
    """Get a dictionary with the number of residues per chainID from a PDB file."""
    u = mda.Universe(pdb_file)
    residue_counts = {}
    
    for chain in u.segments:
        chain_id = chain.segid
        residue_count = len(chain.residues)
        # Keep the maximum residue count for each chain ID
        if chain_id not in residue_counts or residue_count > residue_counts[chain_id]:
            residue_counts[chain_id] = residue_count
        print(f"Chain {chain_id}: {residue_count} residues (keeping max)")
    
    return residue_counts

def parse_chain_ids_from_filename(filename):
    """Extract chain IDs from filename after A1, A2, A3 or A4."""
    # Pattern: something_A[1-4]_chain1_chain2_..._chainN.pdb
    basename = os.path.basename(filename)
    name_without_ext = os.path.splitext(basename)[0]
    
    # Find the pattern A1, A2, A3, or A4 followed by chain IDs
    pattern = r'A[1-4]_(.+)'
    match = re.search(pattern, name_without_ext)
    
    if match:
        chains_part = match.group(1)
        # Split by underscore to get individual chain IDs
        chain_ids = chains_part.split('_')
        return chain_ids
    else:
        print(f"Warning: Could not parse chain IDs from filename {filename}")
        return []

def get_residue_counts_from_pdb_chains(pdb_file):
    """Get residue counts for each chain in a PDB file using atoms instead of segments."""
    u = mda.Universe(pdb_file)
    
    # Group atoms by chainID
    chains = {}
    for atom in u.atoms:
        chain_id = atom.chainID
        if chain_id not in chains:
            chains[chain_id] = set()
        chains[chain_id].add(atom.resid)
    
    # Count unique residues per chain
    residue_counts = {chain_id: len(resids) for chain_id, resids in chains.items()}
    
    return residue_counts

def preprocess_pdb_chains(pdb_file):
    """Preprocess PDB file to fix chain assignments where residue numbering restarts at 1."""
    u = mda.Universe(pdb_file)
    
    print(f"Preprocessing {os.path.basename(pdb_file)} for chain merging...")
    
    # Group atoms by residue number first
    residues_by_number = {}
    for atom in u.atoms:
        resid = atom.resid
        if resid not in residues_by_number:
            residues_by_number[resid] = []
        residues_by_number[resid].append(atom)
    
    # Sort residue numbers
    sorted_resids = sorted(residues_by_number.keys())
    
    # Group consecutive residues into logical chains
    logical_chains = []
    current_chain = []
    
    for i, resid in enumerate(sorted_resids):
        if not current_chain:
            # Start new chain
            current_chain.append(resid)
        else:
            # Check if this residue is consecutive with the previous one
            prev_resid = current_chain[-1]
            if resid == prev_resid + 1:
                # Consecutive residue, add to current chain
                current_chain.append(resid)
            else:
                # Gap in residue numbering, start new chain
                if len(current_chain) > 0:
                    logical_chains.append(current_chain)
                current_chain = [resid]
    
    # Add the last chain
    if current_chain:
        logical_chains.append(current_chain)
    
    print(f"Found {len(logical_chains)} logical chains:")
    for i, chain_resids in enumerate(logical_chains):
        first_resid = chain_resids[0]
        last_resid = chain_resids[-1]
        total_atoms = sum(len(residues_by_number[resid]) for resid in chain_resids)
        print(f"  Chain {i}: residues {first_resid}-{last_resid} ({len(chain_resids)} residues, {total_atoms} atoms)")
    
    # Assign consistent chain IDs and segids to each logical chain
    # First, reset all chainIDs to a temporary value
    for atom in u.atoms:
        atom.chainID = 'X'  # Temporary chainID
    
    # Then assign correct chainIDs and segids to logical chains
    chain_id_counter = 0
    for chain_resids in logical_chains:
        # Use a simple letter-based chain ID
        chain_id = chr(ord('A') + chain_id_counter)
        
        # Collect all atoms in this chain
        chain_atoms = []
        for resid in chain_resids:
            chain_atoms.extend(residues_by_number[resid])
        
        # Set chainIDs for atoms
        for atom in chain_atoms:
            atom.chainID = chain_id
        
        # Set segid for the segment (if it exists)
        # Find segments that contain these atoms
        segments_to_update = set()
        for atom in chain_atoms:
            if hasattr(atom, 'segment') and atom.segment is not None:
                segments_to_update.add(atom.segment)
        
        for segment in segments_to_update:
            segment.segid = chain_id
        
        chain_id_counter += 1
    
    # Create preprocessed filename
    dirname = os.path.dirname(pdb_file)
    basename = os.path.basename(pdb_file)
    preprocessed_basename = f"preprocessed_{basename}"
    preprocessed_file = os.path.join(dirname, preprocessed_basename)
    
    # Write the preprocessed structure
    with mda.Writer(preprocessed_file) as writer:
        writer.write(u.atoms)
    
    print(f"Preprocessed PDB saved as: {preprocessed_file}")
    
    # Reload the universe to ensure segments are created correctly
    u = mda.Universe(preprocessed_file)
    
    return preprocessed_file

def map_wrong_chains_to_correct(ref_counts, pdb_file, expected_chain_ids):
    """Map the current chain IDs in PDB to the correct ones based on residue counts."""
    current_counts = get_residue_counts_from_pdb_chains(pdb_file)
    
    print(f"Current chains in {os.path.basename(pdb_file)}: {current_counts}")
    print(f"Expected chain IDs: {expected_chain_ids}")
    
    # Create mapping from current chain ID to correct chain ID
    mapping = {}
    
    # Create a list of all possible mappings with their differences
    possible_mappings = []
    for expected_id in expected_chain_ids:
        if expected_id in ref_counts:
            expected_count = ref_counts[expected_id]
            for current_id, current_count in current_counts.items():
                if current_id not in mapping.values():  # Don't reuse current chains
                    diff = abs(current_count - expected_count)
                    possible_mappings.append((diff, current_id, expected_id, current_count, expected_count))
    
    # Sort by difference (best matches first)
    possible_mappings.sort()
    
    # Assign mappings, ensuring no current chain is used twice
    used_current_chains = set()
    for diff, current_id, expected_id, current_count, expected_count in possible_mappings:
        if current_id not in used_current_chains and expected_id not in mapping.values():
            mapping[current_id] = expected_id
            used_current_chains.add(current_id)
            print(f"Mapping chain {current_id} ({current_count} residues) -> {expected_id} ({expected_count} residues), diff: {diff}")
    
    return mapping

def fix_pdb_chain_ids(pdb_file, mapping):
    """Fix the chain IDs in a PDB file based on the mapping using MDAnalysis."""
    u = mda.Universe(pdb_file)
    
    # Apply the mapping to change chainIDs
    for current_chain, correct_chain in mapping.items():
        # Select atoms with the current chain ID
        atoms_to_change = u.select_atoms(f"chainID {current_chain}")
        if len(atoms_to_change) > 0:
            # Change the chainID for these atoms
            atoms_to_change.chainIDs = [correct_chain] * len(atoms_to_change)
            print(f"Changed {len(atoms_to_change)} atoms from chain {current_chain} to {correct_chain}")
    
    # Create output filename with "fix" prefix
    dirname = os.path.dirname(pdb_file)
    basename = os.path.basename(pdb_file)
    fixed_basename = f"fix_{basename}"
    output_file = os.path.join(dirname, fixed_basename)
    
    # Write the modified structure
    with mda.Writer(output_file) as writer:
        writer.write(u.atoms)
    
    print(f"Fixed PDB saved as: {output_file}")
    print(f"Applied mapping: {mapping}")
    
    # Reload the universe to verify the changes
    u_fixed = mda.Universe(output_file)
    print(f"Final chains in fixed file:")
    for seg in u_fixed.segments:
        print(f"  Chain {seg.segid}: {len(seg.residues)} residues")
    
    return output_file

def main():
    print("Reading reference PDB to get residue counts per chain...")
    ref_counts = get_residue_counts_from_pdb(ref_pdb)
    print(f"Reference chain counts: {ref_counts}")
    
    # Add custom PsbS chains since they are not defined in ref_pdb
    ref_counts['9'] = 424
    #ref_counts['B1'] = 212
    print(f"Added PsbS chains: A1 (212 residues), B1 (212 residues)")
    print(f"Updated reference chain counts: {ref_counts}")
    
    # Get all PDB files in the directory
    pdb_files = glob.glob(os.path.join(dir, "*.pdb"))
    print(f"Found {len(pdb_files)} PDB files in {dir}")
    
    processed_count = 0
    fixed_count = 0
    
    for pdb_file in pdb_files:
        print(f"\nProcessing {os.path.basename(pdb_file)}...")
        processed_count += 1
        
        # Preprocess the PDB file to fix chain assignments
        preprocessed_file = preprocess_pdb_chains(pdb_file)
        
        # Parse expected chain IDs from filename
        expected_chain_ids = parse_chain_ids_from_filename(pdb_file)
        
        if not expected_chain_ids:
            print(f"Skipping {pdb_file} - could not parse chain IDs")
            continue
        
        # Create mapping from current to correct chain IDs using preprocessed file
        mapping = map_wrong_chains_to_correct(ref_counts, preprocessed_file, expected_chain_ids)
        
        if mapping:
            # Fix the PDB file and get the new filename
            fixed_file = fix_pdb_chain_ids(preprocessed_file, mapping)
            print(f"Original: {pdb_file}")
            print(f"Preprocessed: {preprocessed_file}")
            print(f"Fixed: {fixed_file}")
            fixed_count += 1
        else:
            print(f"No mapping found for {pdb_file}")
    
    print(f"\nSummary: Processed {processed_count} files, fixed {fixed_count} files")

if __name__ == "__main__":
    main()