
function lifetime_analysis_protein_protein(){
  odir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis_psii
  dir5="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs"
  gro=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/sim_1/tem4.pdb
  
  #xtc=test_u  # Special case for 5_8_c_e_f_j_k_p_z - swap chainIDs B and f
  #cp ${odir2}/clust_c075/fix_5_8_c_e_f_j_k_p_z.pdb ${odir2}/clust_c075/temp_step1.pdb
  ## Step 1: Change chainID B to temporary segid 'tmpB'
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step1.pdb -o ${odir2}/clust_c075/temp_step2.pdb -sel "chainID B" -segid "tmpB"
  ## Step 2: Change chainID f to chainID B
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step2.pdb -o ${odir2}/clust_c075/temp_step3.pdb -sel "chainID f" -chain "B"
  ## Step 3: Change segid tmpB to chainID f
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step3.pdb -o ${odir2}/clust_c075/fix_5_8_c_e_f_j_k_p_z.pdb -sel "segid tmpB" -chain "f"
  # Clean up temp files
  rm -f ${odir2}/clust_c075/temp_step*.pdb.xtc
  xtc=analysis.xtc
  sel1="not segid A1 A2 A3 A4 and not resname *GG* *SQ* *MG* *HEME* *PG* W* HOH *HG* *MG* PLQ PL9 LUT VIO XAT NEO NEX BCR"
  sel3="segid A1 A2 A3 A4" # Only chlorophylls and proteins
  cutoff=8
  dt=1 # time step between frames
  min_event_ns=1000
  sim=("sim_1" "sim_2" "sim_3" "sim_4" "sim_5" "sim_6" "sim_7" "sim_8")
  for s in "${sim[@]}"; do 
    cd ${dir5}/${s} 
    echo "Starting contact analysis for chain ${chain}..."
    python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel3}" -sel2 "${sel1}" -o . -prefix lifetime_protein -group_by1 "segids" -group_by2 "segids" > lifetime.log 2>&1 &
  done
}

function extract_binding(){
  dir5="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs"
  sim=("sim_1" "sim_2" "sim_3" "sim_4" "sim_5" "sim_6" "sim_7" "sim_8")
  rm -rf /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/*
  for s in "${sim[@]}"; do 
  cd ${dir5}/${s}
    python /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/extract_binding_trajectories_psii.py \
      -f ${dir5}/sim_1/tem4.pdb \
      -trj analysis.xtc \
      -df lifetime_protein_events_df.csv \
      -odir /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj \
      -p analysis.top \
      -mdp /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/em.mdp \
      -preffix ${s}
  done
}

function binding_pose_pdb(){
  cd chain_${chain}
  # the trajectories are named chain${chain}_${id}_${start}_${end}.xtc
  sel1="chainID A B"
  sel2="chainID ${chain} and (not resname *MG* *HEME* *GG* *SQ* *PG* W* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)" # Only chlorophylls and proteins
  f=${an1}/chain_${chain}/initial_fit_merged.pdb
  trj_arr=($(ls ${an1}/top_lifetimes_pdbs/chain_${chain}_*.xtc))
  # 'counter' is used to generate unique output file names for each trajectory in the loop.
  counter=0
  rm -f chain_${chain}_*.log

  for trj in "${trj_arr[@]}"; do
    counter=$((counter + 1))
    python3 ${an1}/binding_pose.py -f ${f} -trj ${trj} -tpr ${an1}/chain_${chain}/protein.tpr -sel1 "${sel2}" -sel2 "${sel1}" -o chain_${chain}_${counter} --cutoff 0.5 0.75 1 1.5 2 >> chain_${chain}_binding.log 2>&1 &
  done
}

function write_equivalent_binding_sites(){
  # The binding events are in the dir: /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  # The files are named sim_{n_sim}_{PsbSID}_{chains}.* e.g. sim_4_A3_6_7.pdb
  # We generate a csv file in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv
  # With the columns:
  # tag original     new          old_chains count tag_number
  # n_s sim_4_A4_N_S sim_4_A4_n_s N_S        6     1

  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/write_equivalent_binding_sites.py
  
  # Then simply copies the first {original}.*pdb and {original}.*tpr files in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  # With the name 
  # {tag_number}_{tag}
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj 
  csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/grouped_binding_pdbs.py ${dir} ${csv} ${odir}
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/concatenated_grouped_trajectories.py ${dir} ${csv} ${odir}


  ocsv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/write_unique_basenames.py ${csv} ${ocsv}
}

function align_trajectories(){
  # Read *grouped.pdb files, select not segids A1... and BB to align the trajectories
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv
  idir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned

  # Make the output directory if it doesn't exist
  mkdir -p ${odir}

  #Remove *log and *xtc files from the output directory
  rm -f ${odir}/*log ${odir}/*xtc

  n_lines=$(wc -l < ${csv})
  for (( i=1; i<=n_lines; i++ )); do
    if [[ $i -eq 1 ]]; then
      continue # Skip the header line
    fi
    #Print the first and second column of the i-th line
    line=$(sed -n "${i}p" ${csv})
    if [[ -n "${line}" ]]; then
      chains=$(echo ${line} | awk '{print $2}')
      basename=$(echo ${line} | awk '{print $1}')

      # Split by underscores
      chains_arr=(${chains//_/ })

      echo "Aligning ${basename} using chains: ${chains_arr[*]}"
      python ${script}/align_structures.py -mobile ${idir}/${basename}.xtc -mobiletop ${idir}/${basename}.pdb -ref ${ref_pdb} -sel "name BB and chainID ${chains_arr[*]}" -o ${odir}/${basename}.xtc > ${odir}/${basename}_align.log 2>&1 &
    fi
    cp ${idir}/*pdb ${odir}/
    cp ${idir}/*tpr ${odir}/

  done
}

function write_occupancy(){
  #Returns: /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned/occupancy.csv

  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  total_frames=$((4950*32)) # 32 PsbS copies, 4950 frames each
  python3 ${script}/write_occupancy.py ${dir} ${total_frames}
}

function binding_pose_grouped(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  idir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_cluster
  trj_arr=($(ls ${idir}/*.xtc))

  mkdir -p ${odir}
  rm -rf ${odir}/*

  for trj in "${trj_arr[@]}"; do
    basename=$(basename ${trj} .xtc)
    f=${idir}/${basename}.pdb
    tpr=${idir}/${basename}.tpr
    special_basename=6_8_c_e_f_j_k_p_z
    special_selection="chainID 8 and name BB"
    # the trajectories are named chain${chain}_${id}_${start}_${end}.xtc
    sel1="segid A1 A2 A3 A4 and name BB" # PsbS
    sel2="(not segid A1 A2 A3 A4 and (not resname *MG* *HEM* *GG* DGD *SQ* *PG* W* *HG* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)) and name BB" # Only chlorophylls and proteins
  
    # 'counter' is used to generate unique output file names for each trajectory in the loop.
    mkdir -p ${odir}/${basename}
    if [[ "${basename}" == *"${special_basename}"* ]]; then
      sel2="${special_selection}"
      echo "Processing ${basename} with special selection: ${sel2}"
      python3 ${script}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -osel "all" -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/${basename} --cutoff 0.75 >> ${odir}/${basename}/cluster.log 2>&1 &
    else
      echo "Processing ${basename} with standard selection..."
      python3 ${script}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -osel "all" -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/${basename} --cutoff 0.75 >> ${odir}/${basename}/cluster.log 2>&1 &
    fi
  done
}

function lifetime_analysis_grouped (){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  rm -f ${odir}/*
  tpr_files=($(ls ${dir}/*.tpr))
  for tpr_file in "${tpr_files[@]}"; do
    basename=$(basename ${tpr_file} .tpr)
    
    f=${dir}/${basename}.pdb
    trj=${dir}/${basename}.xtc
    chains_part=$(echo ${basename} | cut -d'_' -f2-)
    chains_arr=(${chains_part//_/ })
    sel1="segid A*" # PsbS
    sel2="not segid A* and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG*)" # Only chlorophylls, HEME and proteins, carotenoids
    cutoff=8
    dt=2 # time step between frames
    min_event_ns=100


    #python3 ${script}/lifetime_analysis.py -prefix test -n_frames 100 -dt ${dt} -min_event_ns 0 -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix test_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    #Contacts PsbS and chains
    python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix psbs_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    for chain in "${chains_arr[@]}"; do
      # Contacts chains and PsbS
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "chainID ${chain}" -sel2 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "resids" -group_by2 "segids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
    done
  done
}

function extract_cluster(){
  # Returns: /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/biggest_clusters_c075 
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  idir1=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned
  idir2=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_cluster
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  trj_arr=($(ls ${idir1}/*.xtc))

  mkdir -p ${odir}
  rm -rf ${odir}/*

  for trj in "${trj_arr[@]}"; do
    basename=$(basename ${trj} .xtc)
    f=${idir1}/${basename}.pdb
    trj=${idir1}/${basename}.xtc
    log=${idir2}/${basename}/clust_c075/cluster.log
    python3 ${script}/extract_cluster.py -f ${f} -trj ${trj} -g ${log} -o ${odir}/${basename}.pdb
  done
}

function cg2at(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  cg2at_path=/martini/rubiz/thylakoid/scripts/para/bin/cg2at
  scripts=/martini/rubiz/thylakoid/scripts
  rewrite=${1:-"false"}  # Default to false if not provided
  if [ "$rewrite" == "true" ]; then
    rm -rf ${odir}/*
  fi
  
  files=("$dir"/*.pdb)

  cd $odir
  for file in "${files[@]}"; do
      basename=$(basename ${file} .pdb)  # Get filename without path and extension

      o=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at/${basename}/FINAL/final_cg2at_de_novo.pdb
      
      # Skip if output exists
      if [ -f ${o} ]; then
        echo "Skipping ${file}, output already exists."
        continue
      else
        echo "Processing ${file}..."
        rm -rf ${odir}/${basename}/*


        sel="not resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR"
        cofactors="resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR"

        python3 ${script}/sel_to_ndx.py -f ${file} -sel "${sel}" -name "Protein" -o ${odir}/${basename}_protein_cg.ndx
        python3 ${script}/sel_to_ndx.py -f ${file} -sel "${cofactors}" -name "Cofactors" -o ${odir}/${basename}_cofactors_cg.ndx

        gmx editconf -f ${file} -n ${odir}/${basename}_protein_cg.ndx -o ${odir}/${basename}_protein_cg.pdb
        
        # Some proteins have no cofactors, check if ndx is empty
        if [ -s ${odir}/${basename}_cofactors_cg.ndx ]; then
          gmx editconf -f ${file} -n ${odir}/${basename}_cofactors_cg.ndx -o ${odir}/${basename}_cofactors_cg.pdb
        fi

        # cg2at
        ${cg2at_path} -c ${basename}_protein_cg.pdb -ff charmm36-jul2020-updated -fg martini_3-0_charmm36 -w tip3p -loc ${basename} >> ${odir}/${basename}.log 2>&1 &
      fi
  done
}

function add_bonds_pdb_cofactors() {
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  dir1=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at

  files=("$dir1"/*.pdb)
  cd $odir 

  for file in "${files[@]}"; do
    basename=$(basename ${file} .pdb)
    #TODO Copy top file and comment the chains.
    f=${dir}/${basename}.pdb
    trj=${dir}/${basename}.xtc
    chains_part=$(echo ${basename} | cut -d'_' -f2-)
    top_file=${dir1}/*${chains_part}.top
    echo ${top_file}
    #chains_arr=(${chains_part//_/ })


    #gmx genconf -f ${file} -o ${file} -s ${odir}/${basename}_cofactors.tpr
  done

  


}

function check_sucess_cg2at(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at

  files=("$dir"/*.pdb)
  rm -rf ${odir}/cg2at.log
  for file in "${files[@]}"; do
    basename=$(basename ${file} .pdb)  # Get filename without path and extension
    file=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at/${basename}/FINAL/final_cg2at_de_novo.pdb
    # Check if the file exists and is not empty
    if [[ -s ${file} ]]; then
      echo "${file} SUCCESS" >> ${odir}/cg2at.log
    else
      echo "${file} FAILED" >> ${odir}/cg2at.log
    fi
  done
  grep "FAILED" ${odir}/cg2at.log
}

function reassign_chains(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  script=/martini/rubiz/thylakoid/scripts
  an1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  
  cd ${dir}
  for i in `find -name "final_cg2at_de_novo.pdb"`;do 
    subdirpath=$(dirname "$(dirname "$i")")
    subdirpath="${subdirpath#./}"
    #python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${dir}/${subdirpath}/FINAL/final.pdb  -ref ${dir}/${subdirpath}_protein_cg.pdb -i ${dir}/${subdirpath}/FINAL/final_cg2at_de_novo.pdb
    # Align PDB to input CG structure
    chains_arr=(${subdirpath//_/ })
    # Remove the first element (the tag number)
    chains_arr=("${chains_arr[@]:1}")
    echo "Aligning ${subdirpath} using chains: ${chains_arr[*]}"

    ref_pdb=${dir}/${subdirpath}_protein_cg.pdb
    #python ${script}/align_structures.py -mobile ${dir}/${subdirpath}/FINAL/final.pdb -ref ${ref_pdb} -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o ${dir}/${subdirpath}/FINAL/final_aligned.pdb > ${dir}/${subdirpath}_align.log 2>&1 &
    python ${an1}/align_structures.py -mobile ${dir}/${subdirpath}/FINAL/final.pdb -ref ${ref_pdb} -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o ${dir}/${subdirpath}/FINAL/final_aligned.pdb     


  done
}

function lifetimes_to_pdb_psii(){
  pdb_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  lifetimes_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/pdbs_lifetimes
  sel_protein="not resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR"
  sel_cofactors="resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR"
  sel_psbs="chainID 9"
  rm -rf ${odir}/*
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_to_pdb_psii.py ${pdb_dir} ${lifetimes_dir} ${odir} "${sel_protein}" "${sel_cofactors}" "${sel_psbs}"
}

function main(){
  set -e  

  #lifetime_analysis_protein_protein  # Get the binding events a csv file.
  #sleep 80m

  #extract_binding                    # Extract binding events (pdb, xtc, tpr)
  #sleep 10m
  
  #write_equivalent_binding_sites     # Group binding sites
  #write_occupancy                    # !!!Change "total_frames" if the trajectory is extended

  #lifetime_analysis_grouped          # Calculate contacts for each subtrajectory
  #sleep 80m
  #plot_lifetimes                     #TODO
  
  #align_trajectories             
  #sleep 20m
  
  #binding_pose_grouped               # Clustering analysis. !!!Change "special selection" if the trajectory is extended
  #sleep 30m
  #extract_cluster                    # Extract middle cluster as gmx cluster generates corrupted PDBs

  #cg2at
  #sleep 60m
  add_bonds_pdb_cofactors
  #check_sucess_cg2at
  #reassign_chains
  #lifetimes_to_pdb_psii
}

main
