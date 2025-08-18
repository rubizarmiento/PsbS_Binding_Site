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
scripts=${scripts}

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
# 1) Fits previously generated DAFT files to the PSII crystal structure.
function fit_daft_files() {
    chains=("4" "s" "r" "c")
    rotations=20
    seeds=(242 484)
    cd ${dir5}
    for chain in "${chains[@]}"; do 
        mkdir -p chain_${chain}
        cd chain_${chain}
        for rot in $(seq 1 1 ${rotations}); do
            for seed in "${seeds[@]}"; do 
                rotation=$(printf "%04d" ${rot})
                dir=${dir4}/chain_${chain}/seed${seed}/chain_${chain}/initial-${rotation}
                pdb_mobile=${dir}/initial.pdb
                echo "Processing file ${pdb_mobile}"
                #align_structures -reference_f ${pdb0} -reference_s ${tpr0} -target ${pdb_mobile} -o chain_${chain}_seed${seed}_${rotation}.pdb -sel "name BB and chainID ${chain}" 
                align_structures_mda -ref ${pdb0} -mobile ${pdb_mobile} -o chain_${chain}_seed${seed}_${rotation}.pdb -sel "name BB and chainID ${chain}"        
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
    fit_daft_files
    # 2) Fits them to the PSII crystal structure.
    # 3) Excludes all the structures colliding with the PSII crystal structure.
    # 4) 4 chains are being studied, n PSII conformations including PsbS attached to four chains are generated.
}

main