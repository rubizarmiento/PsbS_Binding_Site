"""
Workflow:
    -Read the HELIX lines in /martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb
    -odir: /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/bendix
    -Write .dat files with helix residue ranges for each chain in format:
     first_residue last_residue first_residue last_residue ...
"""

import os
from collections import defaultdict

def parse_helix_records(pdb_file):
    """Parse HELIX records from PDB file and organize by chain"""
    helices_by_chain = defaultdict(list)
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HELIX'):
                # Parse HELIX record according to PDB format
                helix_serial = line[7:10].strip()
                helix_id = line[11:14].strip()
                start_res_name = line[15:18].strip()
                start_chain = line[19].strip()
                start_res_num = int(line[21:25].strip())
                end_res_name = line[27:30].strip()
                end_chain = line[31].strip()
                end_res_num = int(line[33:37].strip())
                helix_class = int(line[38:40].strip()) if line[38:40].strip() else None
                length = int(line[71:76].strip()) if line[71:76].strip() else None
                
                # Store helix info by chain
                helix_info = {
                    'serial': helix_serial,
                    'id': helix_id,
                    'start_res': start_res_num,
                    'end_res': end_res_num,
                    'start_res_name': start_res_name,
                    'end_res_name': end_res_name,
                    'class': helix_class,
                    'length': length
                }
                
                # Filter out non-alpha helices - only include class 1 (right-handed alpha helices)
                if helix_class == 1:  # Only include class 1 (alpha helices)
                    # Assuming start and end are on the same chain (typical case)
                    if start_chain == end_chain:
                        helices_by_chain[start_chain].append(helix_info)
                    else:
                        print(f"Warning: Helix {helix_id} spans multiple chains ({start_chain}-{end_chain})")
                        helices_by_chain[start_chain].append(helix_info)
                else:
                    # Skip all non-alpha helices (classes 2,3,4,5,6,7,8,9,10)
                    helix_type_names = {
                        2: "right-handed omega",
                        3: "right-handed pi", 
                        4: "right-handed gamma",
                        5: "right-handed 3₁₀",
                        6: "left-handed alpha",
                        7: "left-handed omega", 
                        8: "left-handed gamma",
                        9: "2₇ ribbon",
                        10: "polyproline"
                    }
                    helix_type = helix_type_names.get(helix_class, f"unknown class {helix_class}")
                    print(f"Skipping {helix_type} helix (class {helix_class}): {helix_id} in chain {start_chain}")
    
    return helices_by_chain

def write_helix_dat_files(helices_by_chain, output_dir):
    """Write .dat files with helix residue ranges for each chain"""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    for chain, helices in helices_by_chain.items():
        # Sort helices by start residue number
        helices.sort(key=lambda x: x['start_res'])
        
        # Create output filename
        output_file = os.path.join(output_dir, f"helix_chain_{chain}.dat")
        
        with open(output_file, 'w') as f:
            # Write header comment
            f.write(f"# Alpha helix residue ranges for chain {chain} (class 1 only)\n")
            f.write(f"# Format: first_residue last_residue first_residue last_residue ...\n")
            f.write(f"# Note: Only class 1 (right-handed alpha) helices included - all other helix types filtered out\n")
            
            # Write helix ranges on one line
            helix_ranges = []
            for helix in helices:
                helix_ranges.extend([str(helix['start_res']), str(helix['end_res'])])
            
            f.write(" ".join(helix_ranges) + "\n")
            
            # Write detailed info as comments
            f.write(f"\n# Detailed helix information:\n")
            for helix in helices:
                f.write(f"# Helix {helix['id']}: {helix['start_res_name']}{helix['start_res']}-{helix['end_res_name']}{helix['end_res']} "
                       f"(class {helix['class']}, length {helix['length']})\n")
        
        print(f"Written helix data for chain {chain} to {output_file}")
        print(f"  Found {len(helices)} helices with ranges: {' '.join(helix_ranges)}")

def main():
    # Input and output paths
    pdb_file = "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb"
    output_dir = "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/bendix/alpha_only"
    
    print(f"Reading HELIX records from: {pdb_file}")
    print(f"Output directory: {output_dir}")
    print("Filtering: Only including class 1 (right-handed alpha) helices, ignoring all other helix types (classes 2-10)")
    
    # Parse helix records
    helices_by_chain = parse_helix_records(pdb_file)
    
    if not helices_by_chain:
        print("No HELIX records found in the PDB file!")
        return
    
    print(f"Found helices in {len(helices_by_chain)} chains: {list(helices_by_chain.keys())}")
    
    # Write .dat files
    write_helix_dat_files(helices_by_chain, output_dir)
    
    print("Done!")

if __name__ == "__main__":
    main()
