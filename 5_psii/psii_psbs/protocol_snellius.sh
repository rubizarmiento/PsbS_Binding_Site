#---INSIDE MARTINI--
function sync_martini_to_snellius {
    directories=(
        "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs"
        "/martini/rubiz/shared/to_Tugba/"
        "/martini/rubiz/thylakoid/proteins_thylakoid/M3_cofactors/"
        "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/"
        "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/sh_snel"
        "/martini/rubiz/thylakoid/templates/"
        "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/chains_opm/"
        "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/"
        "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/"
    )
    for directory in "${directories[@]}"
    do
        cd $directory
        #echo $PWD
        #Change /martini/rubiz to /scratch-shared/rubil
        remote_directory=$(echo $directory | sed 's|/martini/rubiz|/scratch-shared/rubil|g')
        rsync -aP ${directory}/ rubil@snellius.surf.nl:${remote_directory}/
        #echo "${directory}/ ${remote_directory}/"
    done
}

function sync_snellius_to_martini {
    directories=(
        "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs"
        "/martini/rubiz/shared/to_Tugba/"
        "/martini/rubiz/thylakoid/proteins_thylakoid/M3_cofactors/"
        "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/"
        "/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/sh_snel"
        "/martini/rubiz/thylakoid/templates/"
        "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/chains_opm/"
        "/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/"
        "/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/"
    )
    for directory in "${directories[@]}"
    do
        cd $directory
        #echo $PWD
        #Change /martini/rubiz to /scratch-shared/rubil
        remote_directory=$(echo $directory | sed 's|/martini/rubiz|/scratch-shared/rubil|g')
        rsync -aP rubil@snellius.surf.nl:${remote_directory}/ ${directory}/
        #echo "${directory}/ ${remote_directory}/"
    done
}



#---INSIDE SNELLIUS--
function mkdir_snellius() {
    directories=(
        "/scratch-shared/rubil/Github/PsbS_Binding_Site/5_psii/psii_psbs"
        "/scratch-shared/rubil/shared/to_Tugba/"
        "/scratch-shared/rubil/thylakoid/proteins_thylakoid/M3_cofactors/"
        "/scratch-shared/rubil/PsbS_Binding_Site/4_pairs/mdps/"
        "/scratch-shared/Github/PsbS_Binding_Site/4_pairs/sh_snel"
        "/scratch-shared/Github/thylakoid/templates/"
        "/scratch-shared/rubil/Github/PsbS_Binding_Site/3_reference_proteins/chains_opm/"
        "/scratch-shared/rubil/Github/PsbS_Binding_Site/3_reference_proteins/psbs/"
        "/scratch-shared/rubil/Github/PsbS_Binding_Site/5_psii/base_dir/"
    )
    for directory in "${directories[@]}"
    do
        mkdir -p $directory
    done
}

function snellius_top {
    cd "/scratch-shared/rubil/PsbS_Binding_Site/psii_psbs/"
    n_files=8
    for ((i=0; i<n_files; i++)); do
        cd "sim_$((i+1))"
        bash /home/rubil/thylakoid/scripts/make_snell.sh ions.top
    done
}

function load_gromacs {
    module purge
    module load 2024
    module load GROMACS/2024.3-foss-2024a
}

function equllibrate_snel {
    cd "/scratch-shared/rubil/PsbS_Binding_Site/psii_psbs/"
    seed=242
    gro=ions_min.gro
    tpr=2-eq.tpr
    n_files=8
    for ((i=0; i<n_files; i++)); do
        cd "sim_$((i+1))"    
        gmx grompp -f /scratch-shared/rubil/PsbS_Binding_Site/4_pairs/mdps/eq_seed${seed}.mdp -c ${gro} -p ions.top -o ${tpr} -maxwarn 1000 -n index.ndx
    done

}

function submit_simulations {
    file=run_long2.sh
    cd "/scratch-shared/rubil/PsbS_Binding_Site/psii_psbs/"
    n_files=8
    for ((i=0; i<n_files; i++)); do
        cd "sim_$((i+1))"    
        sbatch ${file}
    done
}

#sync_martini_to_snellius
sync_snellius_to_martini