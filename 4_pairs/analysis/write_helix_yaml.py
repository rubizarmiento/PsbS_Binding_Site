"""
Workflow:
    -Read the HELIX lines from PDB file(s)
    -odir: /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/
    -Write YAML file with helix information organized by chain and helix number:
     Structure:
       chain_n:
         H1:
           start: 1
           end: 20
         H2:
           start: 30
           end: 50
       chain_s:
         H1:
           start: 1
           end: 20
"""

import yaml
import os
from collections import defaultdict
import argparse

def parse_helix_records(pdb_file, include_all_helices=False):
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
                
                # Filter based on helix class
                if include_all_helices:
                    # Include all helix classes
                    if start_chain == end_chain:
                        helices_by_chain[start_chain].append(helix_info)
                    else:
                        print(f"Warning: Helix {helix_id} spans multiple chains ({start_chain}-{end_chain})")
                        helices_by_chain[start_chain].append(helix_info)
                else:
                    # Only include class 1 (alpha helices)
                    if helix_class == 1:
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

def write_helix_yaml(helices_by_chain, output_file):
    """Write YAML file with helix information organized by chain and helix number"""
    
    yaml_data = {}
    
    for chain in sorted(helices_by_chain.keys()):
        helices = helices_by_chain[chain]
        
        # Sort helices by start residue number
        helices.sort(key=lambda x: x['start_res'])
        
        chain_data = {}
        
        # Assign helix numbers (H1, H2, H3, etc.)
        for helix_num, helix in enumerate(helices, start=1):
            helix_key = f"H{helix_num}"
            chain_data[helix_key] = {
                'start': helix['start_res'],
                'end': helix['end_res'],
                'start_res_name': helix['start_res_name'],
                'end_res_name': helix['end_res_name'],
                'length': helix['length'],
                'class': helix['class']
            }
        
        yaml_data[f"chain_{chain}"] = chain_data
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Write YAML file with nice formatting
    with open(output_file, 'w') as f:
        f.write("# Helix residue ranges organized by chain and helix number\n")
        f.write("# Format: chain_ID -> H1, H2, H3, etc -> start/end residue numbers\n")
        f.write("# Note: Only class 1 (right-handed alpha) helices included\n\n")
        
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    
    print(f"Written helix YAML to: {output_file}")
    
    # Print summary
    print("\nSummary:")
    for chain in sorted(yaml_data.keys()):
        helices = yaml_data[chain]
        print(f"  {chain}: {len(helices)} helices")
        for helix_id, helix_data in helices.items():
            print(f"    {helix_id}: residues {helix_data['start']}-{helix_data['end']}")

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description='Parse HELIX records from PDB file and write to YAML format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Use default values
  python write_helix_yaml.py
  
  # Specify custom PDB file and output location
  python write_helix_yaml.py -p /path/to/custom.pdb -o /path/to/output.yaml
  
  # Include all helix types (not just alpha helices)
  python write_helix_yaml.py --all-helices
        '''
    )
    
    parser.add_argument(
        '-f', '--pdb-file',
        type=str,
        default="/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb",
        help='Path to PDB file (default: 5XNL.pdb)'
    )
    
    parser.add_argument(
        '-o', '--output-file',
        type=str,
        default="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/helix_definitions.yaml",
        help='Path to output YAML file (default: helix_definitions.yaml)'
    )
    
    parser.add_argument(
        '--all-helices',
        action='store_true',
        help='Include all helix types, not just class 1 (alpha) helices'
    )
    
    args = parser.parse_args()
    
    print(f"Reading HELIX records from: {args.pdb_file}")
    print(f"Output file: {args.output_file}")
    
    if args.all_helices:
        print("Including ALL helix types (classes 1-10)")
    else:
        print("Filtering: Only including class 1 (right-handed alpha) helices, ignoring all other helix types (classes 2-10)\n")
    
    # Parse helix records
    helices_by_chain = parse_helix_records(args.pdb_file, include_all_helices=args.all_helices)
    
    if not helices_by_chain:
        print("No HELIX records found in the PDB file!")
        return
    
    print(f"Found helices in {len(helices_by_chain)} chains: {list(helices_by_chain.keys())}\n")
    
    # Write YAML file
    write_helix_yaml(helices_by_chain, args.output_file)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
