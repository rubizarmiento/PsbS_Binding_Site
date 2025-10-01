
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

function gen_test_trajectories(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_test
  mkdir -p ${odir}
  trj_files=($(ls ${dir}/*.xtc))
  for trj_file in "${trj_files[@]}"; do
    basename=$(basename ${trj_file} .xtc)
    echo -e "0\n" | gmx trjconv -f ${trj_file} -s ${dir}/${basename}.tpr -o ${odir}/${basename}.xtc -e 50000 
  done
}

function lifetime_analysis(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  tpr_files=($(ls ${dir}/*.tpr))
  for tpr_file in "${tpr_files[@]}"; do
    basename=$(basename ${tpr_file} .tpr)
    f=${dir}/${basename}.pdb
    trj=${dir}/${basename}.xtc
    chains_part=$(echo ${basename} | cut -d'_' -f4-)
    chains_arr=(${chains_part//_/ })
    sel1="segid A*" # PsbS
    sel2="segid ${chains_part} and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG*)" # Only chlorophylls, HEME and proteins, carotenoids

    cutoff=8
    dt=1 # time step between frames
    min_event_ns=100
    python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir} -prefix psbs_${basename} -group_by1 "segids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    for chain in "${chains_arr[@]}"; do
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel2 "segid ${chain}" -sel1 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "segids" -group_by2 "resids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel2 "segid ${chain}" -sel1 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "segids" -group_by2 "resids" > chain_${chain}_${basename}.log 2>&1 &
    done
  done
}

function clean_csv(){
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes   
  sed -i '/A1/d' ${odir}/*residue_summary_df.csv
  sed -i '/A2/d' ${odir}/*residue_summary_df.csv
  sed -i '/A3/d' ${odir}/*residue_summary_df.csv
  sed -i '/A4/d' ${odir}/*residue_summary_df.csv
}

function clustering(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  if [[ "${odir}" != "/" && -d "${odir}" ]]; then
    rm -rf ${odir}/*
  else
    echo "Error: Refusing to delete contents of root or non-existent directory: ${odir}"
    exit 1
  fi
  tpr_files=($(ls ${dir}/*.tpr))
  echo "Found ${#tpr_files[@]} TPR files"
  counter=0
  for tpr_file in "${tpr_files[@]}"; do
    basename=$(basename ${tpr_file} .tpr)
    counter=$((counter + 1))
    mkdir -p ${odir}/${counter}_${basename}
    cd ${odir}/${counter}_${basename} || { echo "Failed to cd to ${odir}/${counter}_${basename}"; continue; }
    echo "Processing ${basename} in $(pwd)"
    f=${dir}/${basename}.pdb
    trj=${dir}/${basename}.xtc
    tpr=${dir}/${basename}.tpr
    chains_part=$(echo ${basename} | cut -d'_' -f4-)
    chains_arr=(${chains_part//_/ })
    echo "Chains: ${chains_arr[@]}"
    sel1="segid A* and name BB"
    sel2="segid ${chains_arr[@]} and (not resname *GG* *SQ* *PG* W* HOH *HG* *MG* *HEME* PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB CHL)"
    osel="segid A* ${chains_arr[@]}"
    python3 ${script}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -sel1 "${sel2}" -sel2 "${sel1}" -osel "${osel}" -o . --cutoff 0.5 0.75 >> ${odir}/${counter}_${basename}.log 2>&1 &
    echo "Started Python for ${basename} (PID: $!)"
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

function group_centers(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering 
  odir2=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering_grouped   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir2}
  rm -rf ${odir2}/*
  dirs=( $(find "${odir}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;) )
  array=("clust_c05" "clust_c075")
  for arr in "${array[@]}"; do
    mkdir -p ${odir2}/${arr}
    for dir_ in "${dirs[@]}"; do
      basename=$(basename ${dir_})
      cp ${odir}/${dir_}/${arr}/centers.pdb ${odir2}/${arr}/${basename}.pdb
    done
  done
  python3 ${script}/n_clusters_psii.py ${odir2} 
  python3 ${script}/dict_5xnl.py ${odir2}/clust_c075 # Tries to correct centers.pdb but they files are too corrupted, later are just simply extracted
}

function align_structures(){
  #TODO: Implement
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/clustering_grouped/clust_c075 
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  pdb_files=($(ls ${dir}/fix*.pdb))
  cd $dir
  for file in "${pdb_files[@]}"; do
    rm -f ${basename}_aligned.pdb
    basename=$(basename ${file} .pdb)
    chains_part=$(echo ${basename} | cut -d'_' -f5-)
    chains_arr=(${chains_part//_/ })
    #Remove the elements that start with A1, A2, A3, A4 
    chains_arr=($(echo "${chains_arr[@]}" | tr ' ' '\n' | grep -vE '^(A1|A2|A3|A4)' | tr '\n' ' '))
    #echo ${file}
    #echo ${chains_arr[@]} 
    #Remove fix_ from the basename
    basename=${basename#fix_}
    python ${script}/align_structures.py -mobile ${file} -ref ${ref_pdb} -sel "name BB and chainID ${chains_arr[*]}" -o ${basename}_aligned.pdb > ${basename}_align.log 2>&1 &
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
}

function align_trajectories(){
  # Read *grouped.pdb files, select not segids A1... and BB to align the trajectories
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv
  odir
  n_lines=$(wc -l < ${csv})
  for (( i=1; i<=n_lines; i++ )); do
    if [[ $i -eq 1 ]]; then
      continue # Skip the header line
    fi
    #Print the first and second column of the i-th line
    line=$(sed -n "${i}p" ${csv})
    if [[ -n "${line}" ]]; then
      #echo ${line}
      chains=$(echo ${line} | awk '{print $1}')
      # Split by underscores
      chains_arr=(${chains//_/ })
      basename=$(echo ${line} | awk '{print $2}') # Second column (original name)
      echo ${basename}
      #python ${script}/align_structures.py -mobile ${basename}.xtc -mobiletop ${basename}.pdb -ref ${ref_pdb} -sel "name BB and chainID ${chains_arr[*]}" -o ${basename}_aligned.xtc > ${basename}_align.log 2>&1 &
    fi
  done
}

function concatenate_trajectories(){
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/concatenated_grouped_trajectories.py
}

function binding_pose_grouped(){
  an1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  trj_arr=($(ls /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/*.xtc))
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering
  for trj in "${trj_arr[@]}"; do
    basename=$(basename ${trj} .xtc)
    rm -f ${odir}/${basename}/cluster.log
    f=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/${basename}.pdb
    tpr=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/${basename}.tpr
    special_basename=5_8_c_e_f_j_k_p_z
    special_selection="chainID 5 c 8 and name BB"
    # the trajectories are named chain${chain}_${id}_${start}_${end}.xtc
    sel1="segid A1 A2 A3 A4 and name BB" # PsbS
    sel2="(not segid A1 A2 A3 A4 and (not resname *MG* *HEME* *GG* *SQ* *PG* W* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)) and name BB" # Only chlorophylls and proteins
  
    # 'counter' is used to generate unique output file names for each trajectory in the loop.
    mkdir -p ${odir}/${basename}
    if [[ "${basename}" == *"${special_basename}"* ]]; then
      sel2="${special_selection}"
      echo "Processing ${basename} with special selection: ${sel2}"
      python3 ${an1}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -osel "all" -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/${basename} --cutoff 0.75 >> ${odir}/${basename}/cluster.log 2>&1 &
    else
      echo "Processing ${basename} with standard selection..."
    fi
  done
}

function lifetime_analysis_grouped (){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/lifetimes   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
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
    dt=1 # time step between frames
    min_event_ns=100

    #python3 ${script}/lifetime_analysis.py -prefix test -n_frames 100 -dt ${dt} -min_event_ns 0 -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix test_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &

    #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix psbs_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    for chain in "${chains_arr[@]}"; do
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel2 "segid ${chain}" -sel1 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "segids" -group_by2 "resids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel2 "segid ${chain}" -sel1 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "segids" -group_by2 "resids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
    done
  done
}

function write_occupancy(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  python3 ${script}/write_occupancy.py
}

function group_centers_all(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering 
  odir2=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering_grouped   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir2}
  rm -rf ${odir2}/*
  dirs=( $(find "${odir}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;) )
  array=("clust_c05" "clust_c075")
  for arr in "${array[@]}"; do
    mkdir -p ${odir2}/${arr}
    for dir_ in "${dirs[@]}"; do
      basename=$(basename ${dir_})
      cp ${odir}/${dir_}/${arr}/centers.pdb ${odir2}/${arr}/${basename}.pdb
    done
  done
  python3 ${script}/n_clusters_psii.py ${odir2}
  python3 ${script}/dict_5xnl_v2.py ${odir2}/clust_c075 # Correct the chainIDs in the pdb files.
  # Special case for 5_8_c_e_f_j_k_p_z - swap chainIDs B and f
  #cp ${odir2}/clust_c075/fix_5_8_c_e_f_j_k_p_z.pdb ${odir2}/clust_c075/temp_step1.pdb
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step1.pdb -o ${odir2}/clust_c075/temp_step2.pdb -sel "chainID B" -segid "tmpB"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step2.pdb -o ${odir2}/clust_c075/temp_step3.pdb -sel "chainID f" -chain "B"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step3.pdb -o ${odir2}/clust_c075/fix_5_8_c_e_f_j_k_p_z.pdb -sel "segid tmpB" -chain "f"
  ## Special case for 4_7_8 - swap chainIDs 8 and A
  #cp ${odir2}/clust_c075/fix_4_7_8.pdb ${odir2}/clust_c075/temp_step1.pdb
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step1.pdb -o ${odir2}/clust_c075/temp_step2.pdb -sel "chainID 8" -segid "tmp"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step2.pdb -o ${odir2}/clust_c075/temp_step3.pdb -sel "chainID A" -chain "A"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step3.pdb -o ${odir2}/clust_c075/temp_step4.pdb -sel "chainID B" -chain "8"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step4.pdb -o ${odir2}/clust_c075/fix_4_7_8.pdb -sel "segid tmp" -chain "B"
#
  #cp ${odir2}/fix_5_8_c_e_f_j_k_p_z.pdb ${odir2}/clust_c075/temp_step1.pdb
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step1.pdb -o ${odir2}/clust_c075/temp_step2.pdb -sel "chainID 8" -segid "tmp"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step2.pdb -o ${odir2}/clust_c075/temp_step3.pdb -sel "chainID B" -chain "8"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step3.pdb -o ${odir2}/clust_c075/temp_step4.pdb -sel "chainID A" -chain "B"
  #python3 ${script}/modify_structure.py -f ${odir2}/clust_c075/temp_step4.pdb -o ${odir2}/clust_c075/fix_5_8_c_e_f_j_k_p_z.pdb -sel "segid tmp" -chain "A"
  ## Clean up temp files
  #rm -f ${odir2}/clust_c075/temp_step*.pdb 

}




function align_structures_grouped(){
  #TODO: Implement
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering_grouped/clust_c075 
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  special_basename=5_8_c_e_f_j_k_p_z
  special_selection="8"
  pdb_files=($(ls ${dir}/fix*.pdb))
  cd $dir
  for file in "${pdb_files[@]}"; do
    rm -f ${basename}_aligned.pdb
    basename=$(basename ${file} .pdb)
    chains_part=$(echo ${basename} | cut -d'_' -f3-)
    chains_arr=(${chains_part//_/ })
    #echo ${file}
    #Remove fix_ from the basename
    basename=${basename#fix_}
    if [[ "${basename}" == "${special_basename}" ]]; then
      chains_arr=(${special_selection})
    fi
    echo "Aligning ${basename} using chains: ${chains_arr[*]}"
    python ${script}/align_structures.py -mobile ${file} -ref ${ref_pdb} -sel "name BB and chainID ${chains_arr[*]}" -o ${dir}/${basename}_aligned.pdb > ${dir}/${basename}_align.log 2>&1 &
  done
}

function extract_cluster(){
  python3 extract_middle_structures.py --clustering-dir /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all/clustering --sim-base-dir /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/all --output-dir /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/biggest_clusters_c075 --clust-c075-only --biggest-only
}


function reassign_chains(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/biggest_clusters_c075
  cd ${dir}
  for i in `find -name "final_cg2at_de_novo.pdb"`;do 
    subdirpath=$(dirname "$(dirname "$i")")
    subdirpath="${subdirpath#./}"
    python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${dir}/${subdirpath}/FINAL/final.pdb  -ref ${dir}/${subdirpath}_protein.pdb -i ${dir}/${subdirpath}/FINAL/final_cg2at_de_novo.pdb 
  done
}

function lifetimes_to_pdb_psii(){
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_to_pdb_psii.py
}


function main(){
  #set -e  
  # 1) Get the binding events a csv file.
  #lifetime_analysis_protein_protein 
  
  # 2) Extract binding events (pdb, xtc, tpr)
  #extract_binding

  #DELgen_test_trajectories # Optional, debugging
  
  #(Writes in: /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes)
  #DELlifetime_analysis 
  #DELclean_csv

  #DEL - clustering 
  #DEL - group_centers
  #DEL - align_structures

  # Group binding sites
  #write_equivalent_binding_sites # Check here
  #align_trajectories
  #concatenate_trajectories
  #write_occupancy
  #binding_pose_grouped
  
  #HERE
  #group_centers_all
  #align_structures_grouped
  #extract_cluster
  #lifetime_analysis_grouped

  #reassign_chains
  #lifetimes_to_pdb_psii
  
}

main
