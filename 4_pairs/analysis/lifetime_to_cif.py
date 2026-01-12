"""
Extract B-factors from lifetime analysis and assign to PDB structures.

This module reads PDB files from atomistic simulations (after CG-to-AT conversion)
and assigns B-factors based on residue lifetimes from contact analysis. 

Input Files:
------------
Lifetime CSV files (lifetimes_dir):
  - For PROTEIN/COFACTORS: {lifetimes_dir}/chain_X_{case}_residue_summary_df.csv
    * Where X is the chain identifier (e.g., "4") from case name "chain_4_3"
    * Contains columns: resid, sum_ns (used as B-factor)
  - For PSBS: {lifetimes_dir}/psbs_{case}_residue_summary_df.csv
    * Contains columns: resid, sum_ns (used as B-factor)

Structure PDB files (pdb_dir):
  - PROTEIN: {pdb_dir}/{case}/FINAL/final_aligned.pdb
  - COFACTORS: {pdb_dir}/{case}_cofactors_cg.pdb
  - PSBS: {pdb_dir}/{case}/FINAL/final_aligned.pdb
  * Where {case} is the case directory name (e.g., "chain_4_3")

Output Files:
-------------
mmCIF files (odir):
  - PROTEIN: {odir}/{case}_protein.cif
  - COFACTORS: {odir}/{case}_cofactors.cif
  - PSBS: {odir}/{case}_psbs.cif

The script processes three types of structures:
- PROTEIN: Main protein chains with contact lifetime data per chain
- COFACTORS: Cofactor molecules (pigments, metal ions) with lifetime data per chain
- PSBS: PsbS protein with unified lifetime data across all residues

If B-factor data is missing, all residues default to B-factor=0.
mmCIF format supports B-factors > 999 (PDB format limitation).

Author: Rubi Zarmiento-Garcia

"""

import os
import sys
import warnings
import pandas as pd
import MDAnalysis as mda
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.mmcifio import MMCIFIO

# Suppress MDAnalysis element guessing warnings
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

def validate_input_files(pdb_dir, lifetimes_dir, case, struct_type, chain):
    """
    Validate that required input files exist and raise errors if missing.
    
    Parameters
    ----------
    pdb_dir : str
        Directory containing case subdirectories with PDB files.
    lifetimes_dir : str
        Directory containing residue summary CSV files.
    case : str
        Case identifier (e.g., "chain_4_3").
    struct_type : str
        Structure type: "PROTEIN", "COFACTORS", or "PSBS".
    chain : str
        Chain identifier (e.g., "4").
        
    Raises
    ------
    FileNotFoundError
        If required PDB or CSV files are missing (except for COFACTORS).
    """
    # Check PDB files
    if struct_type in ["PROTEIN", "PSBS"]:
        pdb_file = f"{pdb_dir}/{case}/FINAL/final_aligned.pdb"
        if not os.path.isfile(pdb_file):
            raise FileNotFoundError(f"Required PDB file not found: {pdb_file}")
    elif struct_type == "COFACTORS":
        pdb_file = f"{pdb_dir}/{case}_cofactors_cg.pdb"
        # Don't raise error for cofactors, just skip
        if not os.path.isfile(pdb_file):
            print(f"Cofactors PDB not found (skipping): {pdb_file}")
            return False
    
    # Check CSV files
    if struct_type == "PSBS":
        csv_file = f"{lifetimes_dir}/psbs_{case}_residue_summary_df.csv"
        if not os.path.isfile(csv_file):
            raise FileNotFoundError(f"Required CSV file not found: {csv_file}")
    elif struct_type in ["PROTEIN", "COFACTORS"]:
        csv_file = f"{lifetimes_dir}/chain_{chain}_{case}_residue_summary_df.csv"
        if not os.path.isfile(csv_file):
            if struct_type == "COFACTORS":
                # Don't raise error for cofactors, just skip
                print(f"Cofactors CSV not found (skipping): {csv_file}")
                return False
            else:
                raise FileNotFoundError(f"Required CSV file not found: {csv_file}")
    
    return True

def load_bfactor_mapping(lifetimes_dir, case, chains=None):
    """
    Load residue-to-bfactor mapping from CSV files.
    
    Parameters
    ----------
    lifetimes_dir : str
        Directory containing residue summary CSV files.
    case : str
        Case identifier (e.g., "5_6", "9_c_s_z").
    chains : list of str, optional
        Chain identifiers to load. If None, assumes PSBS case.
        
    Returns
    -------
    dict
        Mapping of (chain, residue ID) tuple to B-factor (sum_ns values) for chains,
        or residue ID to B-factor for PSBS.
        If no data found, returns empty dict (all atoms get B-factor 0).
    """
    resid_to_bfactor = {}
    
    if chains is None:
        # PSBS case: single file, use resid as key
        csv_file = f"{lifetimes_dir}/psbs_{case}_residue_summary_df.csv"
        if not os.path.isfile(csv_file):
            print(f"CSV file not found: {csv_file}. Using B-factor=0 for all residues.")
            return {}
        df = pd.read_csv(csv_file)
        resid_to_bfactor = dict(zip(df['resid'], df['sum_ns']))
        print(f"Loaded {len(resid_to_bfactor)} residue-bfactor mappings from PSBS")
    else:
        # Protein/cofactors case: multiple chain files, use (chain, resid) as key
        for chain in chains:
            csv_file = f"{lifetimes_dir}/chain_{chain}_{case}_residue_summary_df.csv"
            if not os.path.isfile(csv_file):
                print(f"CSV file not found for chain {chain}: {csv_file}. Using B-factor=0.")
                continue
            df_chain = pd.read_csv(csv_file)
            for resid, bfactor in zip(df_chain['resid'], df_chain['sum_ns']):
                # Use (chain, resid) tuple as key to avoid conflicts
                resid_to_bfactor[(chain, resid)] = bfactor
        
        if resid_to_bfactor:
            print(f"Loaded {len(resid_to_bfactor)} residue-bfactor mappings from {len(chains)} chains")
            print(f"  CSV file: {csv_file}")
        else:
            print(f"No B-factor data found for any chains. Using B-factor=0 for all residues.")
    
    return resid_to_bfactor


def assign_bfactors_to_structure(structure, resid_to_bfactor):
    """
    Assign B-factors to atoms in a structure based on residue mapping.
    
    Parameters
    ----------
    structure : Bio.PDB.Structure
        The structure to modify.
    resid_to_bfactor : dict
        Mapping of residue ID (or (chain, residue ID) tuple) to B-factor value.
        
    Returns
    -------
    Bio.PDB.Structure
        The modified structure with assigned B-factors.
    """
    for residue in structure.get_residues():
        residue_id = residue.get_id()[1]
        chain_id = residue.parent.id
        
        # Try (chain, resid) first, then fall back to resid only
        bfactor = resid_to_bfactor.get((chain_id, residue_id), 
                                       resid_to_bfactor.get(residue_id, 0.0))
        
        for atom in residue.get_atoms():
            atom.set_bfactor(bfactor)
            element = atom.element
            if element == "X":
                atom.element = "C"  # Change unknown elements to Carbon
    return structure

def save_structure_to_mmcif(structure, output_file, bonds=None):
    """
    Save a Bio.PDB structure to mmCIF format with optional bond information.
    
    Parameters
    ----------
    structure : Bio.PDB.Structure
        The structure to save.
    output_file : str
        Output file path (should end with .cif).
    bonds : list of tuple, optional
        List of bonds as tuples of atom indices (1-based).
        
    Notes
    -----
    Bond information is written to the _struct_conn category in mmCIF format.
    Bonds are manually appended to the mmCIF file after initial save.
    """
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)
    
    # Add bonds if provided
    if bonds:
        # Read the saved file and append bond information
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Create bond entries in mmCIF format
        bond_lines = ["#\nloop_\n"]
        bond_lines.append("_struct_conn.id\n")
        bond_lines.append("_struct_conn.conn_type_id\n")
        bond_lines.append("_struct_conn.ptnr1_label_atom_id\n")
        bond_lines.append("_struct_conn.ptnr1_label_comp_id\n")
        bond_lines.append("_struct_conn.ptnr1_label_asym_id\n")
        bond_lines.append("_struct_conn.ptnr1_label_seq_id\n")
        bond_lines.append("_struct_conn.ptnr2_label_atom_id\n")
        bond_lines.append("_struct_conn.ptnr2_label_comp_id\n")
        bond_lines.append("_struct_conn.ptnr2_label_asym_id\n")
        bond_lines.append("_struct_conn.ptnr2_label_seq_id\n")
        
        atom_list = list(structure.get_atoms())
        bonds_added = 0
        for idx, (atom1_idx, atom2_idx) in enumerate(bonds):
            try:
                atom1 = atom_list[atom1_idx - 1]
                atom2 = atom_list[atom2_idx - 1]
                
                # Format as single line: id type atom1name res1name chain1 res1id atom2name res2name chain2 res2id
                bond_record = (
                    f"b{idx + 1} covale {atom1.name} {atom1.parent.resname} "
                    f"{atom1.parent.parent.id} {atom1.parent.id[1]} "
                    f"{atom2.name} {atom2.parent.resname} "
                    f"{atom2.parent.parent.id} {atom2.parent.id[1]}\n"
                )
                bond_lines.append(bond_record)
                bonds_added += 1
            except (IndexError, AttributeError):
                # Silently skip bonds that reference atoms outside the Bio.PDB structure
                continue
        
        # Append bonds to file
        with open(output_file, 'a') as f:
            f.writelines(bond_lines)
        
        print(f"Wrote mmCIF to: {output_file} with {bonds_added} bonds")
    else:
        print(f"Wrote mmCIF to: {output_file}")

def save_structure_to_pdb(universe, output_file):
    """
    Save a structure to PDB format using MDAnalysis with bonds capped at 12 per atom.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis Universe to save (includes structure and bonds).
    output_file : str
        Output file path (should end with .pdb).
        
    Notes
    -----
    MDAnalysis writes CONECT records automatically.
    Bonds are post-processed to cap at maximum 12 per atom for VMD compatibility.
    B-factors are preserved from the universe atoms.
    """
    # Write PDB with all bonds first
    universe.atoms.write(output_file, bonds='all')
    
    print(f"Wrote PDB to: {output_file}")

def main():
    """
    Main function to extract lifetimes and assign to PDB structures.
    
    Processes all case directories and outputs mmCIF files with B-factors.
    """
    # Input
    pdb_dir = sys.argv[1]
    lifetimes_dir = sys.argv[2]
    odir = sys.argv[3]
    sel_protein = sys.argv[4]
    sel_cofactors = sys.argv[5]
    sel_psbs = sys.argv[6]
    basenames_csv = sys.argv[7]  # CSV file with exact case names to process

    # Create output directory
    os.makedirs(odir, exist_ok=True)
    
    # Read case names from CSV file (one per line)
    if not os.path.isfile(basenames_csv):
        raise FileNotFoundError(f"Basenames CSV file not found: {basenames_csv}")
    
    with open(basenames_csv, 'r') as f:
        cases = [line.strip() for line in f if line.strip()]
    
    if not cases:
        raise ValueError(f"No case names found in {basenames_csv}")
    
    print(f"Found {len(cases)} cases from {basenames_csv}.")
    
    # Configuration for three structure types
    config = {
        "PROTEIN": {
            "pdb_pattern": "{pdb_dir}/{case}/FINAL/final_aligned.pdb",
            "output": "{odir}/{case}_protein.cif",
            "format": "mmcif",
        },
        "COFACTORS": {
            "pdb_pattern": "{pdb_dir}/{case}_cofactors_cg.pdb",
            "output": "{odir}/{case}_cofactors.cif",
            "format": "mmcif",
        },
        "PSBS": {
            "pdb_pattern": "{pdb_dir}/{case}/FINAL/final_aligned.pdb",
            "output": "{odir}/{case}_psbs.cif",
            "format": "mmcif",
        }
    }
    
    # Process each structure type for each case
    for struct_type in ["PROTEIN", "COFACTORS", "PSBS"]:
        for case in cases:
            # Get chain identifier from case name (only first component after "chain_")
            # e.g., "chain_4_3" → chain = "4"
            chain = case.split("_")[1]
            
            # Validate input files (raises error for PROTEIN/PSBS, skips for COFACTORS)
            try:
                if not validate_input_files(pdb_dir, lifetimes_dir, case, struct_type, chain):
                    continue  # Skip if cofactors files are missing
            except FileNotFoundError as e:
                print(f"ERROR for {struct_type} {case}: {e}")
                raise
            
            pdb_file = config[struct_type]["pdb_pattern"].format(
                pdb_dir=pdb_dir, case=case
            )
            output_file = config[struct_type]["output"].format(
                odir=odir, case=case
            )
            
            print(f"\n---TYPE: {struct_type} {case}---")
            
            chains = [chain]  # Wrap in list for load_bfactor_mapping compatibility
            
            # Load B-factor mapping and create intermediate PDB with selection
            pdb_output_file = output_file.replace('.cif', '.pdb')
            
            if struct_type == "PSBS":
                resid_to_bfactor = load_bfactor_mapping(lifetimes_dir, case, chains=None)
                # Save PDB file with the selected atoms only
                u = mda.Universe(pdb_file)
                psbs_sel = u.select_atoms(sel_psbs)
                save_structure_to_pdb(psbs_sel, pdb_output_file)
            elif struct_type == "COFACTORS":
                # For cofactors, load B-factor data from chain CSV files
                resid_to_bfactor = load_bfactor_mapping(lifetimes_dir, case, chains=chains)
                # Save PDB file with the selected cofactor atoms only
                u = mda.Universe(pdb_file)
                cofactors_sel = u.select_atoms(sel_cofactors)
                save_structure_to_pdb(cofactors_sel, pdb_output_file)
            else:  # PROTEIN
                resid_to_bfactor = load_bfactor_mapping(lifetimes_dir, case, chains=chains)
                u = mda.Universe(pdb_file)
                protein_sel = u.select_atoms(sel_protein)
                save_structure_to_pdb(protein_sel, pdb_output_file)

            # Save in appropriate format
            if config[struct_type]["format"] == "pdb":
                # For PDB: use MDAnalysis (automatically handles bonds and capping at 12)
                print(f"Loading with MDAnalysis: {pdb_output_file}")
                u = mda.Universe(pdb_output_file)
                
                # Assign B-factors to the universe atoms from residue mapping
                for residue in u.residues:
                    residue_id = residue.resnum
                    
                    # Try chainID first (more reliable), then fall back to segid
                    chain_id = (residue.atoms[0].chainID if hasattr(residue.atoms[0], 'chainID') else '') or residue.segid
                    
                    # Try (chain, resid) first, then fall back to resid only
                    bfactor = resid_to_bfactor.get((chain_id, residue_id), 
                                                   resid_to_bfactor.get(residue_id, 0.0))
                    
                    for atom in residue.atoms:
                        atom.tempfactor = bfactor
                
                # Save using MDAnalysis (automatically includes bonds with 12-bond limit)
                save_structure_to_pdb(u, output_file)
            else:
                # For mmCIF: use Bio.PDB
                print(f"Loading with Bio.PDB: {pdb_output_file}")
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("struct", pdb_output_file)
                
                # Assign B-factors to structure
                structure = assign_bfactors_to_structure(structure, resid_to_bfactor)
                
                # Save as mmCIF
                save_structure_to_mmcif(structure, output_file, bonds=None)

if __name__ == "__main__":
    main()
