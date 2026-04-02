root_dir=/martini/rubiz/Github/PsbS_Binding_Site
analysis_dir=${root_dir}/4_pairs/analysis/analysis_pairs
scripts_dir=${root_dir}/4_pairs/analysis

source ${root_dir}/analysis_dataset/scripts/paths.sh


function lifetime_analysis_protein_protein(){
  odir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis_psii
  dir5="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs"
  gro=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/sim_1/tem4.pdb
  
  xtc=analysis.xtc
  sel1="not segid A1 A2 A3 A4 and not resname *GG* DGD *SQ* *MG* *HEM* *PG* W* HOH *HG* *MG* PLQ PL9 LUT VIO XAT NEO NEX BCR" # Only chlorophylls and proteins
  cutoff=8
  dt=1 # time step between frames
  min_event_ns=1000
  sim=("sim_1" "sim_2" "sim_3" "sim_4" "sim_5" "sim_6" "sim_7" "sim_8")
  for s in "${sim[@]}"; do 
    cd ${dir5}/${s} 
    echo "Starting contact analysis for ${s}..."
    python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel3}" -sel2 "${sel1}" -o . -prefix lifetime_protein -group_by1 "segids" -group_by2 "segids" > lifetime.log 2>&1 &
  done
}

function extract_binding(){
  # Extract all the subtrajectories and generate a tpr file
  # By default it comments the itps after the [ molecules ] section as it genarates problems with gmx trjconv -conect 
  dir5="/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs"
  sim=("sim_1" "sim_2" "sim_3" "sim_4" "sim_5" "sim_6" "sim_7" "sim_8")
  #sim=("sim_1")

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
      -preffix ${s} > ${s}_extract.log 2>&1 &
  done
}

function binding_pose_pdb(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
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
    python3 ${script}/binding_pose.py -f ${f} -trj ${trj} -tpr ${an1}/chain_${chain}/protein.tpr -sel1 "${sel2}" -sel2 "${sel1}" -o chain_${chain}_${counter} --cutoff 0.5 0.75 1 1.5 2 >> chain_${chain}_binding.log 2>&1 &
  done
}

function get_sim_length_in_dir(){
  idir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  trjs=($(ls ${idir}/*.xtc))
  for trj in "${trjs[@]}"; do
    basename=$(basename ${trj} .xtc)
    tpr=${idir}/${basename}.tpr
    #If the last part of the basename is _cofactors, skip
    if [[ ${basename} == *"_cofactors" ]]; then
      echo "Skipping ${basename} as it is a cofactors trajectory"
      continue
    fi
    gmx check -f ${trj} > ${idir}/${basename}_length.log 2>&1
  done
}

function postprocessing_sim_length_in_dir(){
  idir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  ocsv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/simulation_length_ns.csv
  echo "basename,total_length_ns" > ${ocsv}
  log_files=($(ls ${idir}/*_length.log))
  for log in "${log_files[@]}"; do
    time=$(grep "Time          " ${log} | awk '{print $2}')
    basename=$(basename ${log} _length.log)
    echo "${basename},${time}" >> ${ocsv}
  done
  echo "Simulation lengths written to ${ocsv}"
  cat ${ocsv}
}




function write_equivalent_binding_sites(){
  # The binding events are in the dir: /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  # The files are named sim_{n_sim}_{PsbSID}_{chains}.* e.g. sim_4_A3_6_7.pdb
  # We generate a csv file in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv
  # With the columns:
  # tag original     new          old_chains count tag_number
  # n_s sim_4_A4_N_S sim_4_A4_n_s N_S        6     1

  #python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/write_equivalent_binding_sites.py
  
  # Then simply copies the first {original}.*pdb and {original}.*tpr files in /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  # With the name 
  # {tag_number}_{tag}
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj 
  csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_equivalent_chains.csv
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  #python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/grouped_binding_pdbs.py ${dir} ${csv} ${odir}
  #python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/concatenated_grouped_trajectories.py ${dir} ${csv} ${odir}
  ocsv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv # This is the Table used in the Manuscript
  yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/equivalent_chains.yaml
  #python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/write_unique_basenames.py ${csv} ${ocsv} ${yaml}

  #Copy cofactors tpr from trj to trj_grouped using the equivalent_chains mapping
  # Use the first occurrence of each tag in basenames_equivalent_chains.csv
  # to find the cofactors.tpr file from trj/ and copy it with the grouped name
  idir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  declare -A seen_tags
  while read -r tag original new old_chains count tag_number rest; do
    [[ "${tag}" == "tag" ]] && continue  # Skip header
    unique_basename="${tag_number}_${tag}"
    if [[ -z "${seen_tags[${tag}]}" ]]; then
      seen_tags[${tag}]=1
      src_top=${idir}/${original}_cofactors.top
      src_pdb=${idir}/${original}_cofactors.pdb
      dst_top=${odir}/${unique_basename}_cofactors_all.top
      dst_top_2=${odir}/${unique_basename}_cofactors.top
      dst_pdb=${odir}/${unique_basename}_cofactors_all.pdb
      dst_pdb_2=${odir}/${unique_basename}_cofactors.pdb
      if [ -f "${src_top}" ]; then
        # If ${dst_pdb} does not exist continue
        if [ ! -f "${dst_pdb}" ]; then
          echo "Warning: Destination PDB not found: ${dst_pdb}. Skipping ${src_top} and ${src_pdb}"
          continue
        fi
        echo "Copying ${src_top} -> ${dst_top}"
        cp "${src_top}" "${dst_top}"
        # Comment all the lipids in the dest_top_2 file: lines that contain "*MG" "*GG" "*SQ" "*PG"
        lipids="DPPG|DSMG|DSGG|DFGG|FPGG|YPPG|YFPG"
        sed -E "s/^(.*\b(${lipids})\b.*)$/;\1/g" "${dst_top}" > "${dst_top_2}"        
        python3 ${scripts_dir}/sel_to_ndx.py -f ${src_pdb} -sel "not resname *MG* *GG* *SQ* *PG* *HG* *DG*" -name "System" -o ${odir}/${unique_basename}_cofactors.ndx

        
        gmx editconf -f ${src_pdb} -o ${dst_pdb_2} -n ${odir}/${unique_basename}_cofactors.ndx  

      else
        echo "Warning: cofactors TOP not found: ${src_top}"
      fi
      if [ -f "${src_pdb}" ]; then
        echo "Copying ${src_pdb} -> ${dst_pdb}"
        cp "${src_pdb}" "${dst_pdb}"
        gmx grompp -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/em.mdp -c ${dst_pdb} -p ${dst_top} -o ${odir}/${unique_basename}_cofactors_all.tpr -maxwarn 10
        gmx grompp -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/mdps/em.mdp -c ${dst_pdb_2} -p ${dst_top_2} -o ${odir}/${unique_basename}_cofactors.tpr -maxwarn 10
      else
        echo "Warning: cofactors PDB not found: ${src_pdb}"
      fi
    fi
  done < ${csv}
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
  rm -f ${odir}/*log ${odir}/*xtc ${odir}/*pdb ${odir}/*tpr

  special_basenames=("6_8_c_e_f_j_k_p_z")

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

      #If baseanme includes any of the special basenames cange selection to chainID c
      if [[ " ${special_basenames[@]} " =~ " ${chains} " ]]; then
        chains_arr=("c")
        echo "Special case for ${chains}, using chainID C only"
      fi

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

  # Test selections
  trj_arr=("${idir}/6_8_c_e_f_j_k_p_z.xtc")
  special_basenames=("6_8_c_e_f_j_k_p_z")
  special_selection="chainID c and name BB"
  for trj in "${trj_arr[@]}"; do
    basename=$(basename ${trj} .xtc)
    f=${idir}/${basename}.pdb
    tpr=${idir}/${basename}.tpr
    # the trajectories are named chain${chain}_${id}_${start}_${end}.xtc
    sel1="segid A1 A2 A3 A4 and name BB" # PsbS
    sel2="(not segid A1 A2 A3 A4 and (not resname *MG* *HEM* *GG* DGD *SQ* *PG* W* *HG* HOH *MG* *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)) and name BB" # Only chlorophylls and proteins
  
    # 'counter' is used to generate unique output file names for each trajectory in the loop.
    mkdir -p ${odir}/${basename}
    basename_no_id=$(echo ${basename} | cut -d'_' -f2-)
    #echo ${basename_no_id}
    if [[ " ${special_basenames[@]} " =~ " ${basename} " ]]; then
      sel2="${special_selection}"
      echo "Processing ${basename} with special selection: ${sel2}"
      python3 ${script}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -osel "all" -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/${basename} --cutoff 0.45 >> ${odir}/${basename}/cluster.log 2>&1 &
    else
      echo "Processing ${basename} with standard selection..."
      python3 ${script}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -osel "all" -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/${basename} --cutoff 0.45 >> ${odir}/${basename}/cluster.log 2>&1 &
    fi
  done
}

function plot_binding_modes () {
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/PSII_LHCII/psii_with_cofactors_aa.pdb
  binding_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  binding_modes_occupancy_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv
  chains_occupancy_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary/lifetimes_summary_df_bychain_chains_all.csv
  chain_labels_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_labels.yaml
  output1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/psii_binding_sites_overview.eps
  output2=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/psii_binding_sites_overview_all_labels.eps
  output3=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/psii_binding_reaction_center.eps
  output4=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures/psii_binding_reaction_center_all_labels.eps

  vmin=100
  vmax=30000
  cmap='colorcet:CET_L17'
  python3 ${script}/plot_psii_binding_modes.py -min_lifetime 1000 -top_label 4 -log_transform  -alpha_value 0.1 -ref ${ref_pdb} -binding_dir ${binding_dir} -binding_modes_occupancy_csv ${binding_modes_occupancy_csv} -chains_occupancy_csv ${chains_occupancy_csv} -chain_labels_yaml ${chain_labels_yaml} -output ${output1} -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} 
  python3 ${script}/plot_psii_binding_modes.py -min_lifetime 1000 -log_transform  -alpha_value 0.1 -ref ${ref_pdb} -binding_dir ${binding_dir} -binding_modes_occupancy_csv ${binding_modes_occupancy_csv} -chains_occupancy_csv ${chains_occupancy_csv} -chain_labels_yaml ${chain_labels_yaml} -output ${output2} -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} 
  python3 ${script}/plot_psii_binding_modes.py -no_labels -min_lifetime 1000 -exclude_sites "S1" "S2" "S3" "S4" "S11" "S5" "S9" -log_transform  -alpha_value 0.1 -ref ${ref_pdb} -binding_dir ${binding_dir} -binding_modes_occupancy_csv ${binding_modes_occupancy_csv} -chains_occupancy_csv ${chains_occupancy_csv} -chain_labels_yaml ${chain_labels_yaml} -output ${output3} -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} 
  python3 ${script}/plot_psii_binding_modes.py -min_lifetime 1000 -exclude_sites "S1" "S2" "S3" "S4" "S11" "S5" "S9" -log_transform  -alpha_value 0.1 -ref ${ref_pdb} -binding_dir ${binding_dir} -binding_modes_occupancy_csv ${binding_modes_occupancy_csv} -chains_occupancy_csv ${chains_occupancy_csv} -chain_labels_yaml ${chain_labels_yaml} -output ${output4} -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} 
  
}

function plot_psii_venn_diagram () {
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  binding_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  occupancy_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned/occupancy.csv
  output_dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures
  python3 ${script}/plot_VennDiagram.py -ref ${ref} -occupancy_csv ${occupancy_csv} -output_dir ${output_dir} -ref ${ref}
}

function check_selections(){
    scripts_dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
    dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
    tpr_files=($(ls ${dir}/*.tpr))
    for tpr_file in "${tpr_files[@]}"; do
      basename=$(basename ${tpr_file} .tpr)
      f=${dir}/${basename}.pdb
      sel="not segid A* and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG* *DS* *DP* *DG*)"
      python3 ${scripts_dir}/return_non_protein_residues.py -f ${f} -sel "${sel}" 
    done
}

function add_segids_to_pdb(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  tpr_files=($(ls ${dir}/*.tpr))
  for tpr_file in "${tpr_files[@]}"; do
    #Ignore the cofactors.tpr
    if [[ ${tpr_file} == *"cofactors"* ]]; then
      echo "Skipping ${tpr_file} as it is a cofactors trajectory"
      continue
    fi
    basename=$(basename ${tpr_file} .tpr)
    f=${dir}/${basename}.pdb
    o=${dir}/${basename}_labelled.pdb
    echo "Adding segids to ${f} and writing to ${o}"
    python3 ${scripts_dir}/add_segids_to_pdb.py \
      --input_pdb ${f} \
      --helix_yaml ${HELIX_DEFINITIONS_YAML_COMBINED_simtype2} \
      --output_pdb ${o}
    sed -i "s/_A//g" ${o} # Remove _A from segids
    sed -i "s/_B//g" ${o} # Remove _B from segids
  done
}

function lifetime_analysis_grouped (){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  #rm -f ${odir}/*
  tpr_files=($(ls ${dir}/*.tpr))
  for tpr_file in "${tpr_files[@]}"; do
    basename=$(basename ${tpr_file} .tpr)
    
    f=${dir}/${basename}_labelled.pdb
    trj=${dir}/${basename}.xtc
    chains_part=$(echo ${basename} | cut -d'_' -f2-)
    chains_arr=(${chains_part//_/ })
    sel1="chainID 9" # PsbS

    sel2="not segid A* and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG* *DS* *DP* *DG*)" # Only chlorophylls, HEME and proteins, carotenoids
    sel7="chainID 9 and segid T* H*" # PsbS and helix regions
    cutoff=8
    dt=1 # time step between frames
    min_event_ns=100


    #TEST_line python3 ${script}/lifetime_analysis.py -prefix test -n_frames 100 -dt ${dt} -min_event_ns 0 -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix test_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    #Contacts PsbS and chains
    #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix psbs_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix bychain_psbs_${basename} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/bychain_psbs_${basename}.log 2>&1 &

    for chain in "${chains_arr[@]}"; do
      # Contacts chains and PsbS
      sel3="chainID ${chain} and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG* *DS* *DP* *DG*)" # Only chlorophylls, carotenoids and proteins
      sel4="chainID ${chain} and (resname CLA CLB CHL, PLQ PL9 LUT VIO XAT NEO NEX BCR HEME)" # Only chlorophylls and carotenoids 
      sel5="chainID ${chain} and segid B E C A D" # Only helices
      sel6="chainID ${chain} and segid L" # Only loops
      sel8="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* LMG DGD CLA CLB CHL, PLQ PL9 LUT VIO XAT NEO NEX BCR HEM*)" # Only proteins

      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir} -prefix bychain_chain_${chain}_${basename} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/bychain_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel8}" -sel2 "${sel1}" -o ${odir} -prefix bychain_protein_chain_${chain}_${basename} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/bychain_protein_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel3}" -o ${odir} -prefix psbs_chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel4}" -sel2 "${sel1}" -o ${odir} -prefix cofactors_chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/cofactors_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel4}" -sel2 "${sel1}" -o ${odir} -prefix bychain_cofactors_chain_${chain}_${basename} -group_by1 "resnames" -group_by2 "chainIDs" > ${odir}/bychain_cofactors_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel5}" -sel2 "${sel1}" -o ${odir} -prefix bychain_helix_chain_${chain}_${basename} -group_by1 "segids" -group_by2 "chainIDs" > ${odir}/bychain_helix_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel6}" -sel2 "${sel1}" -o ${odir} -prefix bychain_loop_chain_${chain}_${basename} -group_by1 "segids" -group_by2 "chainIDs" > ${odir}/bychain_loop_chain_${chain}_${basename}.log 2>&1 &
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel5}" -sel2 "${sel7}" -o ${odir} -prefix bychain_helix_helix_chain_${chain}_${basename} -group_by1 "segids" -group_by2 "chainIDs" > ${odir}/bychain_helix_helix_chain_${chain}_${basename}.log 2>&1 &
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel7}" -sel2 "${sel5}" -o ${odir} -prefix bychain_helix_helix_psbs_chain_${chain}_${basename} -group_by1 "segids" -group_by2 "chainIDs" > ${odir}/bychain_helix_helix_psbs_chain_${chain}_${basename}.log 2>&1 &
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
    log=${idir2}/${basename}/clust_c045/cluster.log
    python3 ${script}/extract_cluster.py -f ${f} -trj ${trj} -g ${log} -o ${odir}/${basename}.pdb
  done
}

function write_yaml_chainid_and_size(){
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/PSII_LHCII/psii_with_cofactors_aa.pdb
  oyaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_sizes.yaml
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/write_yaml_chainid_and_size.py -f ${ref_pdb} -o ${oyaml} 
}

function get_biggest_chain_in_binding_site(){
  yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_sizes.yaml
  csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv
  ocsv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_biggest_chain.csv
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/get_biggest_chain_in_binding_site.py -yaml ${yaml} -csv ${csv} -ocsv ${ocsv}
}


function extract_structures(){
  # Extract CG protein and cofactor PDBs from middle_cluster snapshots.
  # Can be re-run independently (e.g. if ref_pdb changes) without repeating cg2at.
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster   # pdb files: 9_c_s_z.pdb -> ${basename}.pdb 
  tpr_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at          
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  biggest_chain_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_biggest_chain.csv
  files=("$dir"/*.pdb)
  sel_cofactors="resname CLA CLB CHL *HEM* PLQ PL9 LUT VIO XAT NEO NEX BCR"

  echo "Extracting protein/cofactor PDBs for ${#files[@]} files in ${dir}..."

  cd $odir
  for file in "${files[@]}"; do
      basename=$(basename ${file} .pdb)  # Get filename without path and extension
      echo "Processing ${file}..."
      chains_arr=($(echo ${basename} | cut -d'_' -f2- | tr '_' ' '))  # Extract chains from filename and split into array

      # Get biggest chain in the binding site from csv with grep and awk
      biggest_chain=$(grep "${basename}" ${biggest_chain_csv} | awk -F',' '{print $3}')

      echo "Aligning ${basename} to reference using biggest chain ${biggest_chain} for selection..."  
      python ${script}/align_structures.py -mobile ${dir}/${basename}.pdb -ref ${ref_pdb} -sel_ref "name BB and chainID ${biggest_chain}" -sel_mobile "name BB and chainID ${biggest_chain}" -o ${dir}/${basename}.pdb    


      pdb_aa=${odir}/${basename}/FINAL/final_aligned.pdb
      cp ${pdb_aa} ${pdb_aa}.bak
      #Fit to reference using the biggest chain in the binding site for selection
      python ${script}/align_structures.py -mobile ${pdb_aa} -ref ${ref_pdb} -sel_ref "name BB and chainID ${biggest_chain}" -sel_mobile "name CA and chainID ${biggest_chain}" -o ${pdb_aa}
      # Write ndxs
      sel="not resname CLA CLB CHL *HG* *HEM* PLQ PL9 *GG* *SQ* *PG* DGD LMG DSMG LUT VIO XAT NEO NEX W2 HOH BCR"
      cofactors="resname CLA CLB CHL *HG* *HEM* PLQ PL9 *GG* *SQ* *PG* DGD LMG DSMG LUT VIO XAT NEO NEX W2 HOH BCR"

      echo "-----Writing ndx files for ${basename}-----"
      python3 ${script}/sel_to_ndx.py -f ${dir}/${basename}.pdb -sel "all" -name "System" -o ${odir}/${basename}_protein_cg.ndx
      python3 ${script}/sel_to_ndx.py -f ${dir}/${basename}.pdb -sel "${sel}" -name "Protein" -o ${odir}/${basename}_protein_cg.ndx
      python3 ${script}/sel_to_ndx.py -f ${dir}/${basename}.pdb -sel "${cofactors}" -name "Cofactors" -o ${odir}/${basename}_cofactors_cg.ndx
      chains_str=$(echo ${basename} | cut -d'_' -f2- | tr '_' ' ')
      #python3 ${script}/sel_to_ndx.py -f ${ref_pdb} -sel "chainID ${chains_str} and ${sel}" -name "Protein" -o ${odir}/${basename}_protein_aa.ndx

      # Extract proteins
      echo "--- Extracting protein structure for ${basename}..."
      gmx editconf -f ${file} -n ${odir}/${basename}_protein_cg.ndx -o ${odir}/${basename}_protein_cg.pdb
      #gmx editconf -f ${ref_pdb} -n ${odir}/${basename}_protein_aa.ndx -o ${odir}/${basename}_protein_aa.pdb  

      echo "---Extracting cofactor structure for ${basename}..."
      # Generate ndx for the cofactors
      python3 ${script}/sel_to_ndx.py -f ${file} -sel "${sel_cofactors}" -name "Cofactors" -o ${odir}/${basename}_cofactors_cg.ndx

      if [ -s ${odir}/${basename}_cofactors_cg.ndx ]; then
        echo "Cofactors found for ${basename}, extracting..."
        gmx editconf -f ${file} -n ${odir}/${basename}_cofactors_cg.ndx -o ${odir}/${basename}_cofactors_cg_1.pdb
        echo -e "0\n" | gmx trjconv -f ${odir}/${basename}_cofactors_cg_1.pdb -s ${tpr_dir}/${basename}_cofactors.tpr -conect -o ${odir}/${basename}_cofactors_cg.pdb
        python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${odir}/${basename}_cofactors_cg.pdb  -ref ${odir}/${basename}_cofactors_cg_1.pdb -i ${odir}/${basename}_cofactors_cg.pdb
        
        # Some proteins have no cofactors, check if ndx is empty
        if [ -s ${odir}/${basename}_cofactors_cg.ndx ]; then
          gmx editconf -f ${file} -n ${odir}/${basename}_cofactors_cg.ndx -o ${odir}/${basename}_cofactors_cg_nb.pdb # No bonds info, but ok chains
          echo "Cofactors\n" | gmx trjconv -f ${file} -s ${tpr_dir}/${basename}_cofactors.tpr -n ${odir}/${basename}_cofactors_cg.ndx -conect -o ${odir}/${basename}_cofactors_cg_1.pdb # Bonds info but wrong chains
          python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${odir}/${basename}_cofactors_cg.pdb  -ref ${odir}/${basename}_cofactors_cg_nb.pdb -i ${odir}/${basename}_cofactors_cg_1.pdb
        fi
      else
        echo "No cofactors found for ${basename}, skipping cofactor extraction."
      fi
    
  done
}

function cg2at(){
  # Run cg2at conversion. Requires align_chains() to have been run first.
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at          
  cg2at_path=/martini/rubiz/thylakoid/scripts/para/bin/cg2at
  rewrite=${1:-"false"}  # Default to false if not provided
  if [ "$rewrite" == "true" ]; then
    rm -rf ${odir}/*
  fi

  cd $odir
  for file in ${odir}/*_protein_cg.pdb; do
      basename=$(basename ${file} _protein_cg.pdb)
      o=${odir}/${basename}/FINAL/final_cg2at_de_novo.pdb
      
    # Skip if output exists
    if [ -f ${o} ]; then
      echo "Skipping ${basename}, output already exists."
      continue
    else
      echo "Running cg2at for ${basename}..."
      rm -rf ${odir}/${basename}/*
      ${cg2at_path} -c ${basename}_protein_cg.pdb -ff charmm36-jul2020-updated -fg martini_3-0_charmm36 -w tip3p -loc ${basename} >> ${odir}/${basename}.log 2>&1 &
      fi
  done
}

function reassign_chains(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  script=/martini/rubiz/thylakoid/scripts
  an1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  
  cd ${dir}
  for i in `find -name "final_cg2at_de_novo.pdb"`;do 
    subdirpath=$(dirname "$(dirname "$i")")
    subdirpath="${subdirpath#./}"
    python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${dir}/${subdirpath}/FINAL/final.pdb  -ref ${dir}/${subdirpath}_protein_cg.pdb -i ${dir}/${subdirpath}/FINAL/final_cg2at_de_novo.pdb
    
    chains_arr=(${subdirpath//_/ })
    # Remove the first element (the tag number)
    chains_arr=("${chains_arr[@]:1}")
    echo "Aligning ${subdirpath} using chains: ${chains_arr[*]}"
    ref_pdb="${dir}/${subdirpath}_protein_cg.pdb"
    #python ${script}/align_structures.py -mobile ${dir}/${subdirpath}/FINAL/final.pdb -ref ${ref_pdb} -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o ${dir}/${subdirpath}/FINAL/final_aligned.pdb > ${dir}/${subdirpath}_align.log 2>&1 &
    #python ${an1}/align_structures.py -mobile "${dir}/${subdirpath}/FINAL/final.pdb" -ref "${ref_pdb}" -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o "${dir}/${subdirpath}/FINAL/final_aligned.pdb"     
  done
}

function lifetimes_to_cif_psii_old(){
  pdb_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  lifetimes_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes
  basenames_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv
  sel_protein="not resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR and not chainID 9"
  sel_cofactors="resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR"
  sel_psbs="chainID 9"
  rm -rf ${odir}/*pdb
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_to_cif_psii.py ${pdb_dir} ${lifetimes_dir} ${odir} "${sel_protein}" "${sel_cofactors}" "${sel_psbs}" ${basenames_csv}
}

function lifetimes_to_cif_psii(){
  # Adapted from the pairs lifetimes_to_cif function.
  # Uses add_lifetimes_to_cif.py / add_lifetimes_to_pdb.py per binding mode
  # instead of the monolithic lifetime_to_cif_psii.py script.
  source /usr/local/gromacs-2023.4/bin/GMXRC

  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  pdb_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  lifetimes_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  tpr_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes
  basenames_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv

  mkdir -p ${odir}
  rm -rf ${odir}/*

  # Read basenames from first column of CSV (skip header)
  mapfile -t basenames < <(tail -n +2 "${basenames_csv}" | cut -d',' -f1)

  sel_psbs="chainID 9"

  for basename in "${basenames[@]}"; do
    pdb=${pdb_dir}/${basename}/FINAL/final_aligned.pdb
    cofactors_pdb=${pdb_dir}/${basename}_cofactors_cg.pdb
    cofactors_tpr=${tpr_dir}/${basename}_cofactors.tpr

    if [ ! -f "${pdb}" ]; then
      echo "Warning: PDB not found for ${basename}, skipping: ${pdb}"
      continue
    fi

    # Extract chains from basename (strip the leading number)
    chains_part=$(echo ${basename} | cut -d'_' -f2-)
    chains_arr=(${chains_part//_/ })

    # --- Protein: single CIF per case (all chains) ---
    sel_protein="chainID ${chains_arr[*]} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR LMG DGD CHL CLA CLB HEM*)"
    protein_csvs=()
    for chain in "${chains_arr[@]}"; do
      pcsv=${lifetimes_dir}/chain_${chain}_${basename}_respairs_residue_summary_df.csv
      if [ -f "${pcsv}" ]; then
        protein_csvs+=("${pcsv}")
      else
        echo "Warning: No lifetime CSV for chain ${chain} of ${basename}."
      fi
    done

    if [ ${#protein_csvs[@]} -gt 0 ]; then
      echo "Processing ${basename} protein (${#protein_csvs[@]} chain CSVs)..."
      python3 ${script}/add_lifetimes_to_cif.py \
        -f ${pdb} -sel "${sel_protein}" \
        -csv ${protein_csvs[@]} \
        -o ${odir}/${basename}_protein.cif --log_transform
    else
      echo "Warning: No protein lifetime CSVs found for ${basename}, skipping."
    fi

    # --- Cofactors: single PDB with CONECT records (all chains) ---
    cofactors_csvs=()
    for chain in "${chains_arr[@]}"; do
      ccsv=${lifetimes_dir}/cofactors_chain_${chain}_${basename}_residue_summary_df.csv
      if [ -f "${ccsv}" ]; then
        cofactors_csvs+=("${ccsv}")
      fi
    done

    if [ -f "${cofactors_pdb}" ] && [ -f "${cofactors_tpr}" ] && [ ${#cofactors_csvs[@]} -gt 0 ]; then
      echo "Processing ${basename} cofactors (${#cofactors_csvs[@]} chain CSVs)..."
      echo "0" | gmx trjconv -f ${cofactors_pdb} -s ${cofactors_tpr} -conect -o ${odir}/${basename}_cofactors.pdb
      python3 ${script}/add_lifetimes_to_pdb.py \
        -f ${cofactors_pdb} \
        -csv ${cofactors_csvs[@]} \
        -o ${odir}/${basename}_cofactors.pdb --log_transform
    fi

    # --- PsbS: symmetrize dimer lifetimes before CIF generation ---
    psbs_csv=${lifetimes_dir}/psbs_${basename}_residue_summary_df.csv
    psbs_csv_sym=${odir}/${basename}_psbs_symmetrized.csv

    if [ -f "${psbs_csv}" ]; then
      echo "Symmetrizing ${basename} psbs..."
      python3 ${script}/sum_lifetimes_psbs_chains.py \
        -csv ${psbs_csv} -o ${psbs_csv_sym}

      echo "Processing ${basename} psbs..."
      python3 ${script}/add_lifetimes_to_cif.py \
        -f ${pdb} -sel "${sel_psbs}" \
        -csv ${psbs_csv_sym} \
        -o ${odir}/${basename}_psbs.cif --log_transform
    else
      echo "Warning: No PsbS lifetime CSV for ${basename}, skipping."
    fi
  done
}

function lifetimes_statistics_psii(){
  lifetimes_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_statistics
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  rm -f ${odir}/*
  python3 ${script}/lifetimes_statistics.py -lifetimes_dir ${lifetimes_dir} -o ${odir}/lifetimes_statistics.dat 
}

function plot_lifetimes(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  cifs_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes
  basenames_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj/basenames_binding.csv
  chain_labels_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_labels.yaml
  color_config_yaml=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/color_definitions.yaml
  psii_helix_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix_labels.yaml
  psbs_helix_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_helix_labels.yaml
  output_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/8_cifs_lifetimes
  
  # PDB with helices defined (for reference)
  psii_pdbdatabase=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb
  psbs_pdbdatabase=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/4ri2.pdb
  
  # Generate helix YAML files from PDB (if needed)
  #python3 ${script}/write_helix_yaml.py -o ${psii_helix_yaml} -f ${psii_pdbdatabase}
  #python3 ${script}/write_helix_yaml.py -o ${psbs_helix_yaml} -f ${psbs_pdbdatabase}
  
  # Generate protein sequence plots with B-factor coloring and helix annotations
  echo "Generating protein sequence plots with lifetimes visualization..."
  python3 ${script}/plot_lifetimes_sequences.py \
    -d ${cifs_dir} \
    -b ${basenames_csv} \
    -l ${chain_labels_yaml} \
    -c ${color_config_yaml} \
    -p ${psii_helix_yaml} \
    -s ${psbs_helix_yaml} \
    -o ${output_dir} \
    --split-sequences 106
  
  echo "Sequence plots saved to: ${output_dir}"
}

function write_databases(){
HELIX_DEFINITIONS_YAML_GROUP2="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_helix_labels_merged_psii.yaml"
HELIX_DEFINITIONS_YAML_GROUP1="/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix_labels.yaml"
analysis_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites
scripts_dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
wdir=${analysis_dir}
odir=${wdir}/10_database
idir7=${wdir}/lifetimes
mkdir -p ${odir}
python3 ${scripts_dir}/write_databases.py \
  --helix_def_yaml_group1 ${HELIX_DEFINITIONS_YAML_GROUP1} \
  --helix_def_yaml_group2 ${HELIX_DEFINITIONS_YAML_GROUP2} \
  --output ${odir}/database.csv \
  --csv_files "${idir7}/*respairs_events*.csv" \
  --add_labels_colname "sim_type" \
  --add_labels_values "psii_lhcii"

}

function join_databases(){
  file1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/analysis_pairs/chain_4/10_database/database.csv
  file2=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/analysis_pairs/chain_r/10_database/database.csv
  file3=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/analysis_pairs/chain_s/10_database/database.csv
  file4=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/10_database/database.csv
  
  odir=/martini/rubiz/Github/PsbS_Binding_Site/combined_database
  mkdir -p ${odir}
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/join_databases.py \
    --input_files "${file1} ${file2} ${file3} ${file4}" \
    --output ${odir}/combined_database.csv
}

function main(){
  #set -e  

  #lifetime_analysis_protein_protein  # Get the binding events a csv file.

  #extract_binding                    # Extract binding events (pdb, xtc, tpr)
  #get_sim_length_in_dir
  #postprocessing_sim_length_in_dir
  #write_equivalent_binding_sites     # Group binding sites
  #write_occupancy                    # !!!Change "total_frames" if the trajectory is extended
  #check_selections
  #add_segids_to_pdb
  #lifetime_analysis_grouped           # Calculate contacts for each subtrajectory

  #align_trajectories             
  
  #binding_pose_grouped               # Clustering analysis. 
  #extract_cluster                    # Extract middle structure from largest cluster as gmx cluster generates corrupted PDBs
  #write_yaml_chainid_and_size          # Write YAML with chain IDs and sizes for PSII reference structure
  #get_biggest_chain_in_binding_site    # Get the biggest chain in each binding site and write to CSV for plotting purposes (e.g. to color by chain in the binding modes overview)
  extract_structures                 # Extract CG protein/cofactor PDBs (re-runnable if ref_pdb changes)


  #plot_binding_modes
  #plot_psii_venn_diagram
  #cg2at                              # CG-to -AA conversion (uses align_chains output)
  #reassign_chains                    # Reassign chain IDs and align AA to CG backbone
  lifetimes_to_cif_psii              # CIF files allow bfactors > 999 while PDB files do not.
  #lifetimes_to_cif_psii_old         # Old monolithic version (uses lifetime_to_cif_psii.py)
  #lifetimes_statistics_psii          # Max occupancy
  #plot_lifetimes                     # Generate protein sequence plots with B-factor coloring
  #write_databases
  #join_databases
}

main
