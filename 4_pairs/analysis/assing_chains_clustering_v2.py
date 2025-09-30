"""
Workflow:
    1. Read {ref_pdb} to get a dictionary with the number of residues per chainID with MDAnalysis.
    2. For all the pdbs in {dir} get a dictionary with the number of residues per chain.
    3. The chains in the pdbs have a wrong name. To fix it:
        -Infer the segments from resid 1 to the next resid 1 and assign a temporary segid (A, B, C...).
        -First, the chainIDs are indicated in the filename after A1, A2, A3 or A4.
            For example, in the file 14_sim_6_A4_7_8.pdb, the chainIDs are 7 and 8.
        -Then, using the first dictionary and comparing the number of residues of each segment, the correct chainID is assigned.

"""

import MDAnalysis as mda
import os
import glob
import re
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# Configuration
ref_pdb = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb"
dir_path = "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering_grouped/clust_c075"


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


def infer_logical_chains(pdb_file):
    """Infer logical chains from existing chain IDs in the PDB file."""
    u = mda.Universe(pdb_file)

    print(f"Reading existing chains from {os.path.basename(pdb_file)}...")

    # Group atoms by existing chainID
    chains_by_id = {}
    for atom in u.atoms:
        chain_id = atom.chainID
        if chain_id not in chains_by_id:
            chains_by_id[chain_id] = []
        chains_by_id[chain_id].append(atom)

    # Convert to logical chains format
    logical_chains = []
    chain_id_counter = 0

    for chain_id, atoms in chains_by_id.items():
        if not atoms:
            continue

        # Get residue numbers for this chain
        residues_in_chain = set()
        for atom in atoms:
            residues_in_chain.add(atom.resid)

        sorted_resids = sorted(residues_in_chain)
        logical_chains.append(sorted_resids)

        # Count atoms
        total_atoms = len(atoms)
        chain_label = chr(ord('A') + chain_id_counter)
        print(f"  Chain {chain_label}: chainID '{chain_id}', residues {sorted_resids[0]}-{sorted_resids[-1]} ({len(sorted_resids)} residues, {total_atoms} atoms)")

        chain_id_counter += 1

    return logical_chains, None  # Return None for residues_by_number since we don't need it


def assign_temporary_chain_ids(u, logical_chains, residues_by_number=None):
    """Assign temporary chain IDs based on existing chain structure."""
    # The chains are already properly defined in the PDB file
    # We just need to ensure they have consistent temporary IDs for mapping
    print("Chains already defined in PDB file - using existing structure")

    # No need to modify chainIDs since they're already defined
    # Just ensure segments are properly recognized
    pass


def map_chains_to_correct_ids(ref_counts, current_counts, expected_chain_ids):
    """Map current chain IDs to correct ones based on residue counts."""
    print(f"Current chains: {current_counts}")
    print(f"Expected chain IDs: {expected_chain_ids}")
    print(f"Reference counts: {ref_counts}")

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


def fix_pdb_chain_ids(pdb_file, mapping, output_prefix="fix"):
    """Fix the chain IDs in a PDB file based on the mapping."""
    u = mda.Universe(pdb_file)

    # Apply the mapping to change chainIDs
    for current_chain, correct_chain in mapping.items():
        # Handle empty chainID specially
        if current_chain == '':
            atoms_to_change = u.select_atoms("chainID ''")
        else:
            atoms_to_change = u.select_atoms(f"chainID {current_chain}")

        if len(atoms_to_change) > 0:
            # Change the chainID for these atoms
            atoms_to_change.chainIDs = [correct_chain] * len(atoms_to_change)
            print(f"Changed {len(atoms_to_change)} atoms from chain '{current_chain}' to {correct_chain}")

    # Create output filename
    dirname = os.path.dirname(pdb_file)
    basename = os.path.basename(pdb_file)
    name_without_ext = os.path.splitext(basename)[0]
    output_basename = f"{output_prefix}_{name_without_ext}.pdb"
    output_file = os.path.join(dirname, output_basename)

    # Write the modified structure
    with mda.Writer(output_file) as writer:
        writer.write(u.atoms)

    print(f"Fixed PDB saved as: {output_file}")
    print(f"Applied mapping: {mapping}")

    return output_file


def process_single_pdb(pdb_file, ref_counts):
    """Process a single PDB file: infer chains, map to correct IDs, and save fixed version."""
    print(f"\n{'='*60}")
    print(f"Processing {os.path.basename(pdb_file)}")
    print(f"{'='*60}")

    # Parse expected chain IDs from filename
    expected_chain_ids = parse_chain_ids_from_filename(pdb_file)

    if not expected_chain_ids:
        print(f"Skipping {pdb_file} - could not parse chain IDs")
        return None

    # Infer logical chains from existing chain IDs
    logical_chains, _ = infer_logical_chains(pdb_file)

    # Load universe
    u = mda.Universe(pdb_file)

    # Get current chain counts from existing chain IDs
    current_counts = {}
    for atom in u.atoms:
        chain_id = atom.chainID
        if chain_id not in current_counts:
            current_counts[chain_id] = 0
        current_counts[chain_id] += 1

    print(f"Current chain counts: {current_counts}")

    # Map chains to correct IDs
    mapping = map_chains_to_correct_ids(ref_counts, current_counts, expected_chain_ids)

    # Fix the PDB file
    fixed_file = fix_pdb_chain_ids(pdb_file, mapping)

    return fixed_file

    if mapping:
        # Fix the PDB file and get the new filename
        fixed_file = fix_pdb_chain_ids(pdb_file, mapping)
        print(f"Original: {pdb_file}")
        print(f"Fixed: {fixed_file}")
        return fixed_file
    else:
        print(f"No mapping found for {pdb_file}")
        return None


def main():
    """Main function to process all PDB files in the directory."""
    print("Cleaning up previous output files...")
    
    # Remove files starting with 'fix' or 'preprocessing' in the target directory
    cleanup_patterns = ["fix*", "preprocessing*"]
    for pattern in cleanup_patterns:
        files_to_remove = glob.glob(os.path.join(dir_path, pattern))
        for file_path in files_to_remove:
            try:
                os.remove(file_path)
                print(f"Removed: {os.path.basename(file_path)}")
            except OSError as e:
                print(f"Error removing {os.path.basename(file_path)}: {e}")
    
    print("Reading reference PDB to get residue counts per chain...")
    ref_counts = get_residue_counts_from_pdb(ref_pdb)
    print(f"\nReference chain counts: {ref_counts}")

    # Add custom PsbS chains since they are not defined in ref_pdb
    ref_counts['9'] = 424
    print(f"Added PsbS chains: 9 (424 residues)")
    print(f"Updated reference chain counts: {ref_counts}")

    # Get all PDB files in the directory
    pdb_files = glob.glob(os.path.join(dir_path, "*.pdb"))
    print(f"\nFound {len(pdb_files)} PDB files in {dir_path}")

    processed_count = 0
    fixed_count = 0

    for pdb_file in pdb_files:
        processed_count += 1
        fixed_file = process_single_pdb(pdb_file, ref_counts)
        if fixed_file:
            fixed_count += 1

    print(f"\n{'='*60}")
    print(f"SUMMARY: Processed {processed_count} files, fixed {fixed_count} files")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()