#!/usr/bin/env python3
"""
Convert a coarse-grained PDB segment to a cartoon-ready pseudo-protein PDB.

- Picks one bead per residue using a provided MDAnalysis selection string
- Builds N/CA/C atoms per residue at the bead position (small x-offset for N/C)
- (Optional) Assigns DSSP-style codes per residue: H,G,I->H; E,B->E; T->T; S->S; else->C
- (Optional) Writes per-residue values into B-factor (defaults to BB atom B-factors)
- Writes peptide CONECTs so VMD sees a bonded backbone
- Output loads in VMD as Protein: 1 and renders with NewCartoon

Usage:
    python cg2cartoon.py \
        --pdb 1_n_s.pdb --sel "name BB and resid 10:30" \
        --dssp HHHHHHHHHHHHHHHHHHHHH \
        --dssp_file single_chain.dssp \
        --dssp_dir /path/to/dssp/files \
        --bfac bfac.txt \
        --out cg_cartoon_10_30.pdb
"""

import argparse
from pathlib import Path
import sys

import numpy as np
import MDAnalysis as mda

SS_MAP = {
    'H': 'H', 'G': 'H', 'I': 'H',
    'E': 'E', 'B': 'E',
    'T': 'T',
    'S': 'S'
}
# default -> C (coil)

def load_dssp_file(path: Path):
    """
    Read DSSP sequence from a file, handling FASTA-like format.
    Skips header lines starting with '>' and concatenates sequence lines.
    Returns the DSSP sequence string.
    """
    dssp_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue  # Skip empty lines and headers
            dssp_lines.append(line)
    return ''.join(dssp_lines)

def load_dssp_from_dir(dssp_dir: Path, chains):
    """
    Load DSSP sequences for multiple chains from a directory.
    Looks for files matching *_chain_{chain}.dssp pattern.
    Returns dict[chain_id -> dssp_sequence]
    """
    dssp_dict = {}
    for chain in chains:
        # Look for files matching the pattern
        pattern = f"*_chain_{chain}.dssp"
        matching_files = list(dssp_dir.glob(pattern))
        
        if matching_files:
            dssp_file = matching_files[0]  # Take the first match
            try:
                dssp_seq = load_dssp_file(dssp_file)
                dssp_dict[chain] = dssp_seq
                print(f"Loaded DSSP for chain {chain} from {dssp_file.name}")
            except Exception as e:
                print(f"Warning: Could not load DSSP for chain {chain}: {e}")
        else:
            print(f"Warning: No DSSP file found for chain {chain} (looked for {pattern})")
    
    return dssp_dict

def load_bfactors(path: Path):
    """
    Read per-residue B-factors from a simple two-column text file:
        resid value
    Returns dict[int->float]
    """
    bmap = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            r = int(parts[0]); v = float(parts[1])
            bmap[r] = v
    return bmap

def pdb_atom_line(serial, name, resn, chain, resid, x, y, z, occ=1.00, bfac=0.00, element=''):
    # PDB ATOM format (fixed width). Element is optional (right-justified).
    return (f"ATOM  {serial:5d} {name:^4s}{resn:>3s} {chain:1s}"
            f"{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}"
            f"{'':>10s}{element:>2s}")

def pdb_conect_line(i, j):
    return f"CONECT{i:5d}{j:5d}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdb', required=True, help='Input CG PDB')
    ap.add_argument('--sel', required=True, help='MDAnalysis selection string (e.g., "name BB and resid 10:30")')
    ap.add_argument('--bbname', default='BB', help='Backbone bead name (default: BB)')
    ap.add_argument('--chain', default='A', help='Chain ID for output (default: A)')
    ap.add_argument('--chains', default=None, help='Comma-separated list of chains to process (optional; processes all if not specified)')
    ap.add_argument('--resn',  default='ALA', help='Residue name for output (default: ALA)')
    ap.add_argument('--dssp',  default=None, help='DSSP string for the residues (optional)')
    ap.add_argument('--dssp_file', default=None, help='File containing DSSP string (applied to all chains)')
    ap.add_argument('--dssp_dir', default=None, help='Directory containing *_chain_{chain}.dssp files (optional)')
    ap.add_argument('--bfac',  default=None, help='Two-column file "resid value" for B-factors (optional; defaults to BB atom B-factors)')
    ap.add_argument('--offset', type=float, default=0.10, help='N/C offset from CA (default 0.10)')
    ap.add_argument('--out',   required=True, help='Output PDB path')
    args = ap.parse_args()

    outpath = Path(args.out)

    # Load CG universe
    u = mda.Universe(args.pdb)

    # Parse chains argument
    target_chains = None
    if args.chains:
        target_chains = [c.strip() for c in args.chains.split(',')]
        print(f"Processing only chains: {target_chains}")

    # Apply the user-provided selection string
    sel = u.select_atoms(args.sel)

    # Apply chain filter if specified
    if target_chains:
        chain_filter = " or ".join([f"chainID {chain}" for chain in target_chains])
        sel = sel.select_atoms(chain_filter)

    # We want exactly one bead per residue. If multiple per residue, pick the first by index.
    # Build (resid -> first atom) mapping ordered by resid
    res_to_atom = {}
    for a in sel.atoms:
        r = int(a.resid)
        if r not in res_to_atom:
            res_to_atom[r] = a

    residues = sorted(res_to_atom.keys())
    if len(residues) < 2:
        print("ERROR: Need at least 2 residues to build a cartoon trace.", file=sys.stderr)
        sys.exit(1)

    # Get unique chains from selected atoms
    unique_chains = set(sel.atoms.chainIDs)
    unique_chains = sorted([c for c in unique_chains if c.strip()])
    print(f"Found chains in selection: {unique_chains}")

    # Optional B-factor map
    bmap = {}
    if args.bfac:
        bmap = load_bfactors(Path(args.bfac))

    # Load combined DSSP sequence if available
    combined_dssp = None
    if args.dssp_dir:
        dssp_dir_path = Path(args.dssp_dir)
        if dssp_dir_path.exists():
            # Look for combined DSSP file matching output pattern
            combined_pattern = f"*{Path(args.out).stem}*.dssp"
            combined_files = list(dssp_dir_path.glob(combined_pattern))
            if combined_files:
                try:
                    combined_dssp = load_dssp_file(combined_files[0])
                    print(f"Loaded combined DSSP: {len(combined_dssp)} residues from {combined_files[0].name}")
                except Exception as e:
                    print(f"Warning: Could not load combined DSSP: {e}")

    # Process all chains and combine into single output
    all_lines = []
    all_dssp = []
    global_atom_serial = 1
    global_residue_offset = 0

    # HELIX record for the combined structure
    all_lines.append(f"HELIX    1 H1  {args.resn:>3s} {unique_chains[0]}{1:4d} {args.resn:>3s} {unique_chains[-1]}{500:4d}  1{'':30s}{500:5d}")

    for chain_id in unique_chains:
        print(f"Processing chain {chain_id}...")

        # Select atoms for this specific chain
        chain_sel = sel.select_atoms(f"chainID {chain_id}")

        # Build (resid -> first atom) mapping for this chain
        res_to_atom = {}
        for a in chain_sel.atoms:
            r = int(a.resid)
            if r not in res_to_atom:
                res_to_atom[r] = a

        chain_residues = sorted(res_to_atom.keys())
        if len(chain_residues) < 2:
            print(f"Warning: Chain {chain_id} has fewer than 2 residues, skipping.")
            continue

        print(f"Chain {chain_id}: {len(chain_residues)} residues")

        # Load DSSP for this chain - try combined file first, then per-chain
        dssp_seq = None
        if combined_dssp:
            # Extract the portion for this chain from the combined sequence
            chain_start_idx = global_residue_offset
            chain_end_idx = chain_start_idx + len(chain_residues)
            if chain_end_idx <= len(combined_dssp):
                dssp_seq = combined_dssp[chain_start_idx:chain_end_idx]
                print(f"Using DSSP portion for chain {chain_id}: residues {chain_start_idx+1}-{chain_end_idx}")
        else:
            # Fallback to per-chain files
            if args.dssp_dir:
                dssp_dir_path = Path(args.dssp_dir)
                if dssp_dir_path.exists():
                    pattern = f"*_chain_{chain_id}.dssp"
                    matching_files = list(dssp_dir_path.glob(pattern))
                    if matching_files:
                        dssp_file = matching_files[0]
                        try:
                            dssp_seq = load_dssp_file(dssp_file)
                            print(f"Loaded per-chain DSSP for chain {chain_id} from {dssp_file.name}")
                        except Exception as e:
                            print(f"Warning: Could not load DSSP for chain {chain_id}: {e}")
                    else:
                        print(f"Warning: No DSSP file found for chain {chain_id}")

        # Process residues for this chain
        for local_res_idx, r in enumerate(chain_residues):
            a = res_to_atom[r]
            x, y, z = a.position
            bf = float(bmap.get(r, a.tempfactor))

            # Get DSSP code for this residue
            ss_code = 'C'  # Default to coil
            if dssp_seq and local_res_idx < len(dssp_seq):
                raw = dssp_seq[local_res_idx].upper()
                ss_code = SS_MAP.get(raw, 'C')
            elif combined_dssp and (global_residue_offset + local_res_idx) < len(combined_dssp):
                # Use combined DSSP with global indexing
                idx = global_residue_offset + local_res_idx
                raw = combined_dssp[idx].upper()
                ss_code = SS_MAP.get(raw, 'C')

            all_dssp.append(ss_code)

            # Global residue number (continuing from previous chains)
            global_resid = global_residue_offset + local_res_idx + 1

            # Create atoms with global numbering
            all_lines.append(pdb_atom_line(global_atom_serial, "N",  args.resn, chain_id, global_resid, x - args.offset, y, z, 1.00, bf, 'N'))
            global_atom_serial += 1
            all_lines.append(pdb_atom_line(global_atom_serial, "CA", args.resn, chain_id, global_resid, x,                  y, z, 1.00, bf, 'C'))
            global_atom_serial += 1
            all_lines.append(pdb_atom_line(global_atom_serial, "C",  args.resn, chain_id, global_resid, x + args.offset, y, z, 1.00, bf, 'C'))
            global_atom_serial += 1

        # Update global residue offset for next chain
        global_residue_offset += len(chain_residues)

    # Add CONECT records for all atoms
    serials_N = []
    serials_CA = []
    serials_C = []

    # Parse the lines to get serial numbers for CONECT records
    for line in all_lines:
        if line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 3:
                atom_name = parts[2]
                serial = int(parts[1])
                if atom_name == 'N':
                    serials_N.append(serial)
                elif atom_name == 'CA':
                    serials_CA.append(serial)
                elif atom_name == 'C':
                    serials_C.append(serial)

    # Intra-residue CONECTs
    for N, CA, C in zip(serials_N, serials_CA, serials_C):
        all_lines.append(pdb_conect_line(N, CA))
        all_lines.append(pdb_conect_line(CA, C))

    # Peptide CONECTs (only within each chain, not between chains)
    chain_start = 0
    for chain_id in unique_chains:
        chain_sel = sel.select_atoms(f"chainID {chain_id}")
        res_to_atom = {}
        for a in chain_sel.atoms:
            r = int(a.resid)
            if r not in res_to_atom:
                res_to_atom[r] = a

        chain_residues = sorted(res_to_atom.keys())
        if len(chain_residues) < 2:
            continue

        chain_length = len(chain_residues)
        chain_end = chain_start + chain_length

        # Peptide bonds within this chain
        for i in range(chain_start, chain_end - 1):
            if i + 1 < len(serials_C) and i < len(serials_N):
                all_lines.append(pdb_conect_line(serials_C[i], serials_N[i + 1]))

        chain_start = chain_end

    all_lines.append("END")

    # Write the combined PDB file
    outpath.write_text("\n".join(all_lines) + "\n")
    total_atoms = global_atom_serial - 1
    print(f"Wrote {outpath} with {global_residue_offset} residues ({total_atoms} atoms).")

    # Write combined DSSP sequence to file
    dssp_outpath = outpath.with_suffix('.dssp')
    with open(dssp_outpath, 'w') as f:
        f.write(''.join(all_dssp))
    print(f"Wrote combined DSSP sequence to {dssp_outpath}")

    # Helpful hint for VMD
    print("\nVMD tips:")
    print("  mol new", outpath, "type pdb")
    print("  mol representation NewCartoon")
    print("  mol selection {all}")
    print("  mol color Beta   # if you provided B-factors")
    print("  mol addrep top")

    print(f"\nCompleted processing {len(unique_chains)} chains into single output file.")

if __name__ == "__main__":
    main()