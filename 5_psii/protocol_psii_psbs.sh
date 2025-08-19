#Description: Generates PSII conformations from previously generated DAFT conformations, excluding all the structures colliding.

#Workflow: 
# 1) Gets previously generated DAFT structures per chain.
# 2) Fits them to the PSII crystal structure.
# 3) Excludes all the structures colliding with the PSII crystal structure.
# 4) 4 chains are being studied, n PSII conformations including PsbS attached to four chains are generated.

#---GLOBAL VARIABLES---
wdir=/martini/rubiz/Github/PsbS_Binding_Site
dir3=${wdir}/3_reference_proteins
dir4=${wdir}/4_pairs
dir5=${wdir}/5_psii
scripts=${wdir}/scripts

#---IMPORTANT FILES IN THIS SCRIPT
pdb0=${dir5}/base_dir/protein_only.pdb # CG PSII with cofactors and chains 
tpr0=${dir5}/base_dir/em.tpr
#---TRANSFERABLE AUXILIARY FUNCTIONS---
function align_structures_mda() {
    # Align a structure to a reference structure using mda
    # Arguments: -ref reference.gro -mobile target.pdb -o output.pdb [-sel "name CA"]
    # Example: align_structures_mda -reference_f reference.gro -reference_s reference.pdb -target target.pdb -o aligned.pdb -sel "name CA"
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -ref) ref="$2"; shift ;;
            -mobile) mobile="$2"; shift ;;
            -o) o="$2"; shift ;;
            -sel) sel="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; exit 1 ;;
        esac
        shift
    done
    selection=${sel:-"name CA"} # Default selection is the atomistic Backbone if not provided
    echo "SELECTION:" $selection
    python3 /martini/rubiz/thylakoid/scripts/align_structures.py -mobile ${mobile} -ref ${ref} -sel "${sel}" -o ${o}
}


#---STEP1---
# 2) Rotate reference structure to align with the x axis.
function rotate_to_x_axis() {
    gmx editconf -f ${pdb0} -o ${dir5}/base_dir/rotated.pdb -rotate 0 0 -25 
}

function fit_daft_files() {
    chains=("4" "s" "r" "c" "4" "s" "r" "c") # Identical chains, but with different names
    chainsIDs=("4" "s" "r" "c" "8" "S" "R" "C")

    rotations=20
    seeds=(242 484)
    mkdir -p 5_psii/daft_initial
    cd ${dir5}/daft_initial
    n_chains=${#chains[@]}
    for ((i=0; i<n_chains; i++)); do
        chain=${chains[$i]}
        chainID=${chainsIDs[$i]}
        mkdir -p chain_${chainID}
        cd chain_${chainID}
        for rot in $(seq 1 1 ${rotations}); do
            for seed in "${seeds[@]}"; do 
                rotation=$(printf "%04d" ${rot})
                dir=${dir4}/chain_${chain}/seed${seed}/chain_${chain}/initial-${rotation}
                pdb_mobile=${dir}/initial.pdb
                echo "Processing file ${pdb_mobile}"
                # If chain ! chainID, change it
                if [[ ${chain} != ${chainID} ]]; then
                    echo "Changing chain from ${chain} to ${chainID}"
                    python3 ${scripts}/modify_structure.py -f ${pdb_mobile} -o chain_${chainID}_seed${seed}_${rotation}.pdb -chain ${chainID} -sel "chainID ${chain}"
                    pdb_mobile=chain_${chainID}_seed${seed}_${rotation}.pdb
                fi                
                align_structures_mda -ref  ${dir5}/base_dir/rotated.pdb -mobile ${pdb_mobile} -o chain_${chainID}_seed${seed}_${rotation}.pdb -sel "name BB and chainID ${chainID}"        
            done
        done
        cd ..
    done
}

#---STEP3---
# 3) Excludes all the structures colliding with the PSII crystal structure.

#---STEP4---
# 4) 4 chains are being studied, n PSII conformations including PsbS attached to four chains are generated.

function main() {
    set -e
    # 1) Fits previously generated DAFT files to the PSII crystal structure.
    #rotate_to_x_axis
    fit_daft_files
    # 2) Fits them to the PSII crystal structure.
    # 3) Excludes all the structures colliding with the PSII crystal structure.
    # 4) 4 chains are being studied, n PSII conformations including PsbS attached to four chains are generated.
}

main

#---IMPORTANT OUTPUTS---
# ${dir5}/base_dir/rotated.pdb # Rotated PSII structure aligned to the x-axis.