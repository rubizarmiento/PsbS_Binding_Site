#----------------WORKFLOW----------------
# 1. Sync directories from Snellius to Martini
# 2. Drop every n frames from a trajectory using gmx trjconv
# 3. Remove boundary conditions from a trajectory using mdvwhole



#----------------FUNCTIONS----------------

function store_directories {
    local directories_file=$1
    local -a local_directories=()
    
    while read -r line
    do
        local_directories+=("$line")
    done < "$directories_file"
    
    printf '%s\n' "${local_directories[@]}"

    # Usage - capture output back into array
    #readarray -t directories < <(store_directories "your_file.txt")
}

function sync_snellius_to_martini {
    directory=${1}
    #echo $PWD
    #Change /martini/rubiz to /scratch-shared/rubil
    remote_directory=$(echo $directory | sed 's|/martini/rubiz|/scratch-shared/rubil|g')
    rsync -aP rubil@snellius.surf.nl:${remote_directory}/ ${directory}/
}

function drop_every_n_frames {
    #Drop every n frames from a trajectory using gmx trjconv
    #Arguments: directory preffix start_time end_time n 
    #Example: drop_every_n_frames 3-run.xtc 3-run.tpr 1000_frames.xtc 10 true
    traj=$1
    tpr=$2
    output=$3
    n=$4
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    #If the xtc file exists, drop every n frames
    if [ -f ${traj} ]; then
        if [ ${rewrite} = true ]; then
            echo "Dropping every ${n} frames from ${traj} to ${output}"
            echo -e "System" | gmx trjconv -s ${tpr} -f ${traj} -o ${output} -skip ${n} 
        else
            echo "Output file ${output} already exists. Skipping, unless you set rewrite to true."
        fi
    fi
}

function remove_pbc_mdvwhole(){
    #Remove boundary conditions from a trajectory using gmx trjconv
    #Arguments: directory preffix start_time end_time 
    #Example: remove_pbc 3 directory 3-run 250000 500000 
    gro=$1
    traj=$2
    output=$3
    sel=${4:-"all"} # Default selection is "all" if not provided
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    #If the xtc file exists, remove boundary conditions
    if [ -f ${traj} ]; then
        if [ ${rewrite} = true ]; then
            mdvwhole -x ${traj} -f ${gro} -o ${output} -sel "${sel}" -mol True
        fi
    fi 
}

function remove_pbc_mdvwhole_mol(){
    #Remove boundary conditions from a trajectory using gmx trjconv
    #Arguments: directory preffix start_time end_time 
    #Example: remove_pbc 3 directory 3-run 250000 500000 
    gro=$1
    traj=$2
    output=$3
    sel=${4:-"all"} # Default selection is "all" if not provided
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    #If the xtc file exists, remove boundary conditions
    if [ -f ${traj} ]; then
        if [ ${rewrite} = true ]; then
            mdvwhole -x ${traj} -f ${gro} -o ${output} -sel "${sel}" 


            #mdvwhole -x ${traj} -f ${gro} -o ${output} -sel "${sel}" -mol True
        fi
    fi 
}

function gmx_whole() {
    gro=$1
    tpr=$2
    traj=$2
    output=$3
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    if [ -f ${traj} ]; then
        if [ ${rewrite} = true ]; then
            echo "System" | gmx trjconv -s ${tpr} -f ${gro} -o ${output} -pbc whole

        fi
    fi 
}

function tpr_no_nb() {
    #Create a tpr file without non-bonded interactions
    #Arguments: directory preffix start_time end_time 
    #Example: tpr_no_nb 3 directory 3-run 250000 500000 
    gro=$1
    output=$2
    rewrite=${3:-true} # If true, it will overwrite the output file if it exists
    if [ -f ${gro} ]; then
        if [ ${rewrite} = true ]; then
            echo "Creating tpr file without non-bonded interactions from ${gro} to ${output}"
            cp ions.top no_nb.top
            #Delete the lines that contain:
            #"intermolecular_bonds.itp"
            #"dihedrals.itp"
            sed -i '/intermolecular_bonds.itp/d' no_nb.top
            sed -i '/dihedrals.itp/d' no_nb.top
            sed -i "s/scratch-shared\/rubil/martini\/rubiz/g" no_nb.top
            gmx grompp -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/em.mdp -c ${gro} -p no_nb.top -o ${output} -maxwarn 1000 
        fi
    fi 
}

function tpr_protein() {
    #Create a tpr file without non-bonded interactions
    #Arguments: directory preffix start_time end_time 
    #Example: tpr_no_nb 3 directory 3-run 250000 500000 
    chain=${1}
    gro=$2
    top=${3}
    output=$4
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    dir3="/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/chain_${chain}"
    if [ -f ${gro} ]; then
        if [ ${rewrite} = true ]; then
            echo "Creating tpr file without non-bonded interactions from ${gro} to ${output}"
            cp ${dir3}/system.top ${top}
            #Delete the lines that contain:
            #"intermolecular_bonds.itp"
            #"dihedrals.itp"
            sed -i '/intermolecular_bonds.itp/d' ${top}
            sed -i '/dihedrals.itp/d' ${top}
            gmx grompp -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/em.mdp -c ${gro} -p ${top} -o ${output} -maxwarn 1000 
        fi
    fi 
}

function index_file_proteins() {      
    f=${1} 
    python3 /martini/rubiz/thylakoid/scripts/sel_to_ndx.py -f ${f} -sel "chainID A B 4 s r c" "chainID A B" "chainID 4 s r c" -name "proteins" "psbs" "partner" -o proteins.ndx
    python3 /martini/rubiz/thylakoid/scripts/sel_to_ndx.py -f ${f} -sel "chainID A B 4 s r c" "chainID A B 4 s r c" -name "System" "Protein" -o index3.ndx

}

function extract_proteins() {
    #Drop every n frames from a trajectory using gmx trjconv
    #Arguments: directory preffix start_time end_time n 
    #Example: extract_proteins 3-run.xtc 3-run.tpr proteins.xtc proteins.ndx true
    traj=$1
    tpr=$2
    output=$3
    ndx=$4
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    #If the xtc file exists, drop every n frames
    if [ -f ${traj} ]; then
        if [ ${rewrite} = true ]; then
            echo "Extracting protein from ${traj} to ${output}"
            echo -e "proteins" | gmx trjconv -s ${tpr} -f ${traj} -o ${output} -n ${ndx}
        else
            echo "Output file ${output} already exists. Skipping, unless you set rewrite to true."
        fi
    else
        echo "Trajectory file ${traj} does not exist. Skipping extraction."
    fi
}

function align_structures() {
    #TODO: TEST
    # Align a structure to a reference structure
    # Arguments: -reference reference.gro -target target.gro -o output.gro [-sel "name CA"]
    # Example: align_structures -reference reference.gro -target target.gro -o aligned.gro -sel "name CA"
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -reference_f) reference_f="$2"; shift ;;
            -reference_s) reference_s="$2"; shift ;;
            -target) target="$2"; shift ;;
            -o) o="$2"; shift ;;
            -sel) sel="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; exit 1 ;;
        esac
        shift
    done
    #selection=${selection:-"name CA"} # Default selection is the atomistic Backbone if not provided
    echo $sel
    # Create index file for alignment
    python3 /martini/rubiz/thylakoid/scripts/sel_to_ndx.py -f ${reference_f} -sel "${sel}" "all" -name "Backbone" "System" -o align.ndx
    
    # Align the structure with gmx trjconv
    #echo "gmx trjconv -s ${reference_s} -f ${target} -o ${o} -fit rot+trans -n align.ndx"
    echo -e "Backbone\nSystem" | gmx trjconv -s ${reference_s} -f ${target} -o ${o} -fit rot+trans -n align.ndx
    
    # Clean up
    #rm -f align.ndx 
}

function copy_useful_files() {
    dir_a="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis"
    chains=("4" "s" "r" "c")

    for chain in "${chains[@]}"; do 
        dir="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_${chain}/seed242/chain_${chain}/initial-0001"
        cp ${dir}/*top ${dir_a}/chain_${chain}/
        cp ${dir}/*gro ${dir_a}/chain_${chain}/
        cp ${dir}/*pdb ${dir_a}/chain_${chain}/
        cp ${dir}/*tpr ${dir_a}/chain_${chain}/
    done
}


function align_trajectories() {
    #copy_useful_files
    ref_pdb="initial_fit.pdb"
    ref_tpr="fit.tpr"
    traj_0="proteins_5000ns_concat.xtc"
    dir_a="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis"
    output="aligned_5000ns.xtc"
    chains=("r")
    #chains=("s")

    for chain in "${chains[@]}"; do 
        cd ${dir_a}/chain_${chain}
        align_structures -reference_s ${ref_tpr} -reference_f ${ref_pdb} -target ${traj_0} -o aligned_5000ns.xtc -sel "name BB and not chainID A B"
    done
}

function recover_backup() {
    #Rename the file '#3-run.xtc.1#' to '3-run.xtc'
    chains=("c" "r")
    seeds=(242 484) #Two seeds to generate different conformations
    rotations=20
    files=('#3-run2.part0002.xtc.1#'
        '#3-run2.part0002.trr.1#'
        '#3-run2.part0002.edr.1#'
        '#3-run2.part0002.log.1#'
    )
    outputs=(3-run3.xtc
        3-run3.trr
        3-run3.edr
        3-run3.log
    )
    dir4="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs"
    for chain in "${chains[@]}"; do
        for seed in "${seeds[@]}"; do
            for i in $(seq 1 1 ${rotations}); do
                number=$(printf "%04d" ${i}) #For example, if i=1, it will be 0001, if i=10, it will be 0010
                echo "Recovering chain_${chain} with seed ${seed} and rotation ${number}"
                cd ${dir4}/chain_${chain}/seed${seed}/chain_${chain}/initial-${number}
                for j in "${!files[@]}"; do
                    file=${files[$j]}
                    output=${outputs[$j]}
                    if [ -f ${file} ]; then
                        echo "Renaming ${file} to ${output}"
                        mv ${file} ${output}
                    else
                        echo "File ${file} does not exist. Skipping."
                    fi
                done

            done
        done
    done
}

function extract_frames_special() {
    #Extract the first frames except the first 200 files
    #Arguments: directory preffix start_time end_time 
    #Example: extract_frames_special 3-run.xtc 3-run.tpr 2000_frames.xtc 200 true
    traj=$1
    tpr=$2
    output=$3
    b=$4
    rewrite=${5:-true} # If true, it will overwrite the output file if it exists
    dirs=(/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_c/seed484/chain_c/initial-0002
        /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_c/seed484/chain_c/initial-0013
        /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_c/seed242/chain_c/initial-0007
        /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_r/seed484/chain_r/initial-0008
    )
    #If the xtc file exists, drop every n frames
    for dir in "${dirs[@]}"; do
        cd ${dir}
        if [ -f ${traj} ]; then
            if [ ${rewrite} = true ]; then
                echo "Extracting from ${traj} to ${output}"
                echo -e "System" | gmx trjconv -s ${tpr} -f ${traj} -o ${output} -b ${b} 
            else
                echo "Output file ${output} already exists. Skipping, unless you set rewrite to true."
            fi
        fi
    done

}

function change_name() {
    #Change the name of the file
    #Arguments: old_name new_name
    #Example: change_name old_name.txt new_name.txt
    old_name=$1
    new_name=$2
    if [ -f ${old_name} ]; then
        mv ${old_name} ${new_name}
    else
        echo "File ${old_name} does not exist."
    fi
}

function remove_pbc_dimer_martini(){
    #Remove boundary conditions from a trajectory using gmx trjconv
    #Arguments: directory preffix start_time end_time 
    #Example: remove_pbc 3 directory 3-run 250000 500000 
    xtc=$1
    tpr=$2
    o=$3
    n=$4
    #If the xtc file exists, remove boundary conditions
    echo -e "System\n" | gmx trjconv -s ${tpr}  -f ${xtc} -o pbc_step1.xtc -pbc nojump -n ${n}
    echo -e "Protein\nSystem" | gmx trjconv -s ${tpr} -f pbc_step1.xtc -o pbc_step2.xtc -pbc cluster -n ${n}
    echo -e "System" | gmx trjconv -s ${tpr} -f pbc_step2.xtc -o pbc_step3.xtc -pbc mol -n ${n}
    echo -e "Protein\nSystem" | gmx trjconv -s ${tpr} -f pbc_step3.xtc -o ${o} -pbc cluster -n ${n}
    rm -f pbc_step1.xtc pbc_step2.xtc pbc_step3.xtc
}

function individual_trajectories() {
    # Set the script to exit on error
    set -e

    #-----------ITERATE OVER CHAINS AND SEEDS----------------
    dir4="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs"
    #chains=("4" "s" "r" "c")
    chains=("c")
    #chains=("4" "s" "c")
    rotations=20

    ignore_rotaions=() #Rotations to ignore, for example, if you want to skip the first rotation, add it here
    #rotations=1
    seeds=(242) #Two seeds to generate different conformations
    gro_0="eq.gro"
    pdb_proteins="initial.pdb"
    tpr_0="3-run.tpr"
    traj_0="3-run.xtc"
    traj_0="3-run3.xtc"
    traj_0="all.xtc"
    for chain in "${chains[@]}"; do 
        for seed in "${seeds[@]}"; do
            for i in $(seq 1 1 ${rotations}); do
                if [[ " ${ignore_rotaions[@]} " =~ " ${i} " ]]; then
                    echo "Skipping rotation ${i} for chain ${chain} and seed ${seed}"
                    continue
                fi
            
                number=$(printf "%04d" ${i}) #For example, if i=1, it will be 0001, if i=10, it will be 0010
                echo "Analyzing chain_${chain} with seed ${seed} and rotation ${number}"
                cd ${dir4}/chain_${chain}/seed${seed}/chain_${chain}/initial-${number}
                #1) Sync directories from Snellius to Martini
                #sync_snellius_to_martini $(pwd)

                #2) Drop every n frames from a trajectory using gmx trjconv
                #drop_every_n_frames 3-run.xtc 3-run.tpr 200ns.xtc 5 true

                #3) Remove boundary conditions from a trajectory using mdvwhole
                #TOO EXPENSIVE, will do proteins onlyremove_pbc_mdvwhole ${gro_0} "200ns.xtc" "3-run_whole.xtc" "all" true 
                
                #4) TPR for easier visualization
                #tpr_no_nb ${gro_0} ${tpr_0} "no_nb.tpr" true #Vis all
                #tpr_protein ${chain} ${pdb_proteins} proteins.top "protein.tpr" true # Vis only protein
                #tpr_protein ${chain} initial_fit.pdb fit.top "fit.tpr" true #Fit protein to reference structure

                #5) Analysis Index file
                index_file_proteins ${pdb_proteins}
                #6) Extract protein (THIS IS THE ONLY PART THAT IS NEEDED FOR  THE NEW FRAMES)
                extract_proteins ${traj_0} ${tpr_0} proteins_5000.xtc proteins.ndx true
                remove_pbc_dimer_martini proteins_5000.xtc protein.tpr whole_proteins_5000.xtc index3.ndx
                #NOTUSEDremove_pbc_mdvwhole protein.tpr "proteins_5000.xtc" "whole_proteins_5000.xtc" "all" true # Produces corrupted files that do not allow concatenation
                #7 Extract the first n nanoseconds
                echo -e "0\n" | gmx trjconv -s protein.tpr -f whole_proteins_5000.xtc -o whole_proteins_5000.xtc -b 50000 -t0 0
            done
        done
    done
}


function concatenate_trajectories() {
    # Set the script to exit on error
    set -e

    #-----------ITERATE OVER CHAINS AND SEEDS----------------
    dir4="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs"
    #chains=("4" "s" "r" "c")
    #chains=("s" "r" "c")

    chains=("r")

    rotations=20
    #rotations=1
    seeds=(242 484) #Two seeds to generate different conformations
    pdb_proteins="initial.pdb"
    tpr_0="protein.tpr"
    #traj_0="proteins_5000.xtc"
    #traj_0="aligned_5000ns.xtc"
    traj_0="whole_proteins_5000.xtc"
    output="proteins_5000ns_concat.xtc"
    #sim_length_ps=1000000 # 1 microsecond simulation length
    sim_length_ps=4950000 # 10 nanoseconds simulation length
    #dir3="/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/chain_${chains[0]}"
    b=0
    e=${sim_length_ps} #End time in ps

    for chain in "${chains[@]}"; do 
        array_chain=()
        time=0
        time_str="0\n"
        for seed in "${seeds[@]}"; do
            for i in $(seq 1 1 ${rotations}); do
                number=$(printf "%04d" ${i}) #For example, if i=1, it will be 0001, if i=10, it will be 0010
                array_chain+=("${dir4}/chain_${chain}/seed${seed}/chain_${chain}/initial-${number}/${traj_0}") #Add the trajectory to the array
                time=$((${time} + ${sim_length_ps})) #Increment time by sim_length_ps ps for each trajectory
                #time_str+="${time}\n" #Store time as a string
                time_str+="c\n" #Store time as a string

            done
        done
        #echo "${time_str}"
        #gmx trjcat -f "${array_chain[@]}" -o ${dir4}/chain_${chain}/concat.xtc -settime true -b ${b} -e ${e}

        echo -e "${time_str}" | gmx trjcat -f "${array_chain[@]}" -o ${dir4}/analysis/chain_${chain}/${output} -settime "true" 
    done
}



function chain_trajectories() {
    #copy_useful_files
    align_trajectories

}

function special_rename(){
    dirs=(/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_s/seed484/chain_s/initial-0008
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_r/seed484/chain_r/initial-0006
    /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/chain_c/seed242/chain_c/initial-0005
    )

    for dir in "${dirs[@]}"; do
        cd ${dir}
        echo "Changing name in ${dir}"
        change_name 3-run2.part0002.xtc 3-run2.part0003.xtc
        #gmx check -f 3-run2.part0003.xtc
    done
}

function check_trajectories() {
    chains=("s")
    seeds=(242 484) #Two seeds to generate different conformations
    rotations=20
    dir4="/martini/rubiz/Github/PsbS_Binding_Site/4_pairs"
    traj=all.xtc
    for chain in "${chains[@]}"; do 
        for seed in "${seeds[@]}"; do
            for i in $(seq 1 1 ${rotations}); do
                number=$(printf "%04d" ${i}) #For example, if i=1, it will be 0001, if i=10, it will be 0010
                echo "Checking chain_${chain} with seed ${seed} and rotation ${number}"
                cd ${dir4}/chain_${chain}/seed${seed}/chain_${chain}/initial-${number}
                #Check if the trajectory file exists
                #if [ ! -f ${traj} ]; then
                #    echo "Trajectory file ${traj} does not exist. Skipping."
                #    continue
                #fi
                
                gmx check -f ${traj} 
            done
        done
    done
}


function special_cases() {
    echo "Special cases handled."
    #8) Troubleshooting
    #recover_backup

    #9) #Extract the first frames except the first 200 frames
    #extract_frames_special 3-run.xtc 3-run.tpr "3-run2.part0002.xtc" 200000 true

    #10) Rename files
    #special_rename
}

# Call the main function
function main() {
    set -e
    #special_cases

    individual_trajectories
    #concatenate_trajectories
    #chain_trajectories

    #check_trajectories

}
main 
                               
