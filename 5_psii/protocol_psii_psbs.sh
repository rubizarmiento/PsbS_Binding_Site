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
# 3) Mimimize proteins in vacuum 
function create_directories() {
    dir=${dir5}/psii_psbs/initial_structures
    cd "$dir5/psii_psbs/"
    pdb_files=("$dir"/*.pdb)
    n_files=${#pdb_files[@]}
    for ((i=0; i<n_files; i++)); do
        mkdir -p "sim_$((i+1))"
        cd "sim_$((i+1))"
        rm *pdb
        cp "${pdb_files[$i]}" "initial.pdb"
        cp "$dir5/psii_psbs/system.top" .
        cd ..
    done
}

#----SIMULATION----

function minimize () {
    f=${1}
    p=${2}
    o=${3}
    gmx grompp -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/em.mdp -c ${f} -p ${p} -o ${o}.tpr -maxwarn 1000
    gmx mdrun -v -deffnm ${o} -ntmpi 1 -pin on -pinoffset 0 
}

function minimize_all() {
    gmx editconf -f initial.pdb -o initial.pdb -d 20
    minimize initial.pdb system.top min.tpr
}

function insane () {
    python2 /martini/rubiz/shared/to_Tugba/to_Tugba/insane_M3_lipids.py -l YFPG:2 -l YPPG:1 -l FPGG:1 -l DFGG:5 -l FPMG:1 -l DFMG:7 -l FPSQ:3 -f min.gro -center -sol W -x 55 -y 55 -dz 15 -pbc dodecahedron -solr 0.5 -o insane.gro -p insane.top
    gmx editconf -f insane.gro -o insane.gro -resnr 1
    python3 /martini/rubiz/thylakoid/scripts/correct_insane.py -f insane.gro -sel "resname YFPG YPPG FPGG DFGG FPMG DFMG FPSQ" -p insane.top -of insane2.gro -op insane2.top
}

function merge_top_files () {
    reference_top=system.top 
    insane_top=insane2.top


    #Extract the lines after and not including "Protein" from the insane2.top
    temp_top1=$(mktemp)
    awk '/Protein 1/{found=1; next} found' ${insane_top} > ${temp_top1}
    
    #Lines after and including the intermolecular_bonds.itp line
    temp_top2=$(mktemp) 
    awk '/intermolecular_bonds.itp/,0' ${reference_top} > ${temp_top2}                
    
    #Lines before the intermolecular_bonds.itp line 
    temp_top3=$(mktemp)
    awk '/intermolecular_bonds.itp/ {exit} {print}' ${reference_top} > ${temp_top3}
    
    # Merge the files
    cat ${temp_top3} ${temp_top1} ${temp_top2} > system.top

    echo "Preview of system.top:"
    head -20 system.top
    echo "..."
    tail -10 system.top
}

function minimize_membrane () {
    minimize insane2.gro system.top 1-min
}

function add_ions () {
    cp system.top ions.top
    echo -e "W\n" | gmx genion -s 1-min.tpr -o ions.gro -p ions.top -pname NA -nname CL -neutral -conc 0.10 
}

function minimize_ions () {
    minimize ions.gro ions.top ions_min
}


function equilibrium () {
    python3 /martini/rubiz/thylakoid/scripts/sel_to_ndx.py -f ions_min.gro -sel "not resname W NA CL" "all" "resname W NA CL" -name "Membrane" "System" "Solvent" -o index.ndx
    gmx grompp -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/eq_seed242_dt01.mdp -c ions_min.gro -p ions.top -o 2-eq.tpr -maxwarn 1000 -n index.ndx
    #gmx mdrun -v -deffnm 2-eq -ntmpi 3 -pin on -pinoffset 0 
}

function simulation () {
    dir=/martini/rubiz/Github/PsbS_Binding_Site/initial_structures
    cd "$dir5/psii_psbs/"
    pdb_files=("$dir"/*.pdb)
    n_files=${#pdb_files[@]}
    for ((i=0; i<n_files; i++)); do
        cd "sim_$((i+1))"
        # STEP1: Mimimize proteins in vacuum 
        #gmx editconf -f initial.pdb -o initial.pdb -d 20
        #minimize initial.pdb system.top min

        # STEP2: Add membrane 
        insane

        # STEP3: Merge top files
        merge_top_files

        # STEP4: Minimize the membrane
        minimize_membrane

        # STEP5: Add ions
        add_ions

        # STEP6: Minimize ions
        minimize_ions

        # STEP7: Equilibration
        equilibrium
        cd ..
    done

}

function simulation () {
    dir=/martini/rubiz/Github/PsbS_Binding_Site/initial_structures
    cd "$dir5/psii_psbs/"
    pdb_files=("$dir"/*.pdb)
    n_files=${#pdb_files[@]}
    for ((i=0; i<n_files; i++)); do
        cd "sim_$((i+1))"
        # STEP1: Mimimize proteins in vacuum 
        cp /martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/system.top .
        #gmx editconf -f initial.pdb -o initial.pdb -d 20
        #minimize initial.pdb system.top min

        # STEP2: Add membrane 
        insane

        # STEP3: Merge top files
        merge_top_files

        # STEP4: Minimize the membrane
        minimize_membrane

        # STEP5: Add ions
        add_ions

        # STEP6: Minimize ions
        minimize_ions

        # STEP7: Equilibration
        equilibrium
        cd ..
    done

}




function main() {
    set -e
    # 1) Fits previously generated DAFT files to the PSII crystal structure.
    #rotate_to_x_axis
    #fit_daft_files
    
    # 2) Excludes all the structures colliding with the PSII crystal structure.
    #.   4 chains are being studied, n PSII conformations including PsbS attached to four chains are generated.
    #code generate_inital_configurations.ipynb

    # 3) Creare directories
    #create_directories
}  

main

#---IMPORTANT OUTPUTS---
# ${dir5}/base_dir/rotated.pdb # Rotated PSII structure aligned to the x-axis.
