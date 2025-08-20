n_sim=${1}

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
        n_sim=${1}
        cd /martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/psii_psbs/
        cd "sim_$n_sim"
        # STEP1: Mimimize proteins in vacuum 
        cp /martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/system.top .
        gmx editconf -f initial.pdb -o initial.pdb -d 20
        minimize initial.pdb system.top min

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
}

simulation ${n_sim}

