root_dir=/martini/rubiz/Github/PsbS_Binding_Site

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
      -preffix ${s}
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

function plot_psii_binding_modes () {
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  ref=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
  binding_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  binding_modes_occupancy_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped/binding_modes_lifetimes.csv
  chains_occupancy_csv=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary/lifetimes_summary_df_bychain_chains_all.csv
  output_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/figures
  vmin=100
  vmax=50000
  cmap='colorcet:CET_L17'
  python3 ${script}/plot_psii_binding_modes.py -log_transform -ref ${ref} -binding_dir ${binding_dir} -binding_modes_occupancy_csv ${binding_modes_occupancy_csv} -chains_occupancy_csv ${chains_occupancy_csv} -output_dir ${output_dir} -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} -move_site_label "S10" -move_offset "0 20"
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

function lifetime_analysis_grouped (){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes   
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  #rm -f ${odir}/*
  tpr_files=($(ls ${dir}/*.tpr))
  for tpr_file in "${tpr_files[@]}"; do
    basename=$(basename ${tpr_file} .tpr)
    
    f=${dir}/${basename}.pdb
    trj=${dir}/${basename}.xtc
    chains_part=$(echo ${basename} | cut -d'_' -f2-)
    chains_arr=(${chains_part//_/ })
    sel1="segid A*" # PsbS
    sel2="not segid A* and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG* *DS* *DP* *DG*)" # Only chlorophylls, HEME and proteins, carotenoids
    cutoff=8
    dt=1 # time step between frames
    min_event_ns=100


    #TEST_line python3 ${script}/lifetime_analysis.py -prefix test -n_frames 100 -dt ${dt} -min_event_ns 0 -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix test_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    #Contacts PsbS and chains
    #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix psbs_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
    #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix bychain_psbs_${basename} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/bychain_psbs_${basename}.log 2>&1 &

    for chain in "${chains_arr[@]}"; do
      # Contacts chains and PsbS
      sel4="chainID ${chain} and (resname CLA CLB CHL, PLQ PL9 LUT VIO XAT NEO NEX BCR HEME)" # Only chlorophylls and carotenoids 
      sel3="chainID ${chain} and (not resname *GG* *SQ* *PG* *MG* W* HOH *HG* *DS* *DP* *DG*)" # Only chlorophylls, carotenoids and proteins
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir} -prefix bychain_chain_${chain}_${basename} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/bychain_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel3}" -o ${odir} -prefix psbs_chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_chain_${chain}_${basename}.log 2>&1 &
      #python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel4}" -sel2 "${sel1}" -o ${odir} -prefix cofactors_chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/cofactors_chain_${chain}_${basename}.log 2>&1 &
      python3 ${script}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel4}" -sel2 "${sel1}" -o ${odir} -prefix bychain_cofactors_chain_${chain}_${basename} -group_by1 "resnames" -group_by2 "chainIDs" > ${odir}/bychain_cofactors_chain_${chain}_${basename}.log 2>&1 &
    done
  done
}

function sum_csv_lifetimes(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  #rm -f ${odir}/*
  chains=("s" "n" "8" "7" "k" "z" "6" "g")
  
  outputs_arr_psbs=()
  outputs_arr_chains=()
  outputs_arr_bychain_chains=()
  outputs_arr_cofactors=()
  outputs_arr_bychain_cofactors=()

  for chain in "${chains[@]}"; do
    echo "Processing chain: ${chain}"
    outputs_arr_psbs+=("${odir}/lifetimes_summary_df_psbs_${chain}.csv")
    outputs_arr_chains+=("${odir}/lifetimes_summary_df_chain_${chain}.csv")
    outputs_arr_grouped_chains+=("${odir}/lifetimes_summary_df_grouped_chain_${chain}.csv")
    # Check if the cofactors file is not empty before adding to the array
    if [ -s "${odir}/lifetimes_summary_df_cofactors_chain_${chain}.csv" ]; then
      outputs_arr_cofactors+=("${odir}/lifetimes_summary_df_cofactors_chain_${chain}.csv")
    fi

    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_chain_${chain}.csv -prefix chain_${chain}_ # Suffix is *residue_summary_df.csv
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_bychain_chain_${chain}.csv -prefix bychain_chain_${chain}_ # Suffix is *chain_summary_df.csv python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_grouped_psbs_chain_${chain}.csv -prefix bychain_psbs_chain_${chain}_ # Suffix is *chain_summary_df.csv
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_psbs_chain_${chain}.csv -prefix psbs_chain_${chain}_ # Suffix is *chain_summary_df.csv
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_cofactors_chain_${chain}.csv -prefix cofactors_chain_${chain}_  --ignore-missing  # Suffix is *residue_summary_df.csv, ignores empty files for chains without cofactors
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_bychain_cofactors_chain_${chain}.csv -prefix bychain_cofactors_chain_${chain}_ --ignore-missing # Suffix is *chain_summary_df.csv, ignores empty files for chains without cofactors
  done
  python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_psbs.csv -prefix bychain_psbs_ # Suffix is *residue_summary_df.csv

  python3 ${script}/join_csvs.py -c "${outputs_arr_psbs[@]}" -o ${odir}/lifetimes_summary_df_psbs_all.csv
  python3 ${script}/join_csvs.py -c "${outputs_arr_chains[@]}" -o ${odir}/lifetimes_summary_df_chains_all.csv
  python3 ${script}/join_csvs.py -c "${outputs_arr_bychain_chains[@]}" -o ${odir}/lifetimes_summary_df_bychain_chains_all.csv
  python3 ${script}/join_csvs.py -c "${outputs_arr_cofactors[@]}" -o ${odir}/lifetimes_summary_df_cofactors_all.csv
  python3 ${script}/join_csvs.py -c "${outputs_arr_bychain_cofactors[@]}" -o ${odir}/lifetimes_summary_df_bychain_cofactors_all.csv

}

function symmetric_sum_lifetimes_psbs_chains(){
  idir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis

  chains=("s" "n" "8" "7" "k" "z" "6" "g")


  for chain in "${chains[@]}"; do
    echo "Processing chain: ${chain}"
    python  ${script}/sum_lifetimes_psbs_chains.py -csv ${idir}/lifetimes_summary_df_psbs_chain_${chain}.csv -o ${odir}/lifetimes_summary_df_psbs_chain_${chain}_symmetrized.csv

  done
}

function get_lifetimes_per_binding_mode(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  output=${dir}/binding_modes_lifetimes.csv
  dt=1 # time step between frames
  python3 ${script}/get_lifetimes_per_binding_mode.py -d ${dir} -o ${output} -dt ${dt}
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

function cg2at(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster   # pdb files: 9_c_s_z.pdb -> ${basename}.pdb 
  tpr_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_grouped
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at          
  cg2at_path=/martini/rubiz/Github/PsbS_Binding_Site/analysis_dataset/venv/bin/cg2at
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

      python3 ${script}/sel_to_ndx.py -f ${dir}/${basename}.pdb -sel "${sel}" -name "Protein" -o ${odir}/${basename}_protein_cg.ndx
      python3 ${script}/sel_to_ndx.py -f ${dir}/${basename}.pdb -sel "${cofactors}" -name "Cofactors" -o ${odir}/${basename}_cofactors_cg.ndx
      gmx editconf -f ${file} -n ${odir}/${basename}_protein_cg.ndx -o ${odir}/${basename}_protein_cg.pdb
      
      # Check if ndx is empty
      if [ -s ${odir}/${basename}_cofactors_cg.ndx ]; then
        gmx editconf -f ${file} -n ${odir}/${basename}_cofactors_cg.ndx -o ${odir}/${basename}_cofactors_cg_1.pdb
        echo -e "0\n" | gmx trjconv -f ${odir}/${basename}_cofactors_cg_1.pdb -s ${tpr_dir}/${basename}_cofactors.tpr -conect -o ${odir}/${basename}_cofactors_cg.pdb
        python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${odir}/${basename}_cofactors_cg.pdb  -ref ${odir}/${basename}_cofactors_cg_1.pdb -i ${odir}/${basename}_cofactors_cg.pdb
      fi
      
      # Some proteins have no cofactors, check if ndx is empty
      if [ -s ${odir}/${basename}_cofactors_cg.ndx ]; then
        gmx editconf -f ${file} -n ${odir}/${basename}_cofactors_cg.ndx -o ${odir}/${basename}_cofactors_cg_nb.pdb # No bonds info, but ok chains
        echo "Cofactors\n" | gmx trjconv -f ${file} -s ${tpr_dir}/${basename}.tpr -n ${odir}/${basename}_cofactors_cg.ndx -conect -o ${odir}/${basename}_cofactors_cg_1.pdb # Bonds info but wrong chains
        python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${odir}/${basename}_cofactors_cg.pdb  -ref ${odir}/${basename}_cofactors_cg_nb.pdb -i ${odir}/${basename}_cofactors_cg_1.pdb
      fi
      #cg2at
      ${cg2at_path} -ncpus 1 -c ${basename}_protein_cg.pdb -ff charmm36-jul2020-updated -fg martini_3-0_charmm36 -w tip3p -loc ${basename} >> ${odir}/${basename}.log 2>&1 &
    fi
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
    python3 /martini/rubiz/thylakoid/scripts/assing_resid_chain_from_pdb.py -o ${dir}/${subdirpath}/FINAL/final.pdb  -ref ${dir}/${subdirpath}_protein_cg.pdb -i ${dir}/${subdirpath}/FINAL/final_cg2at_de_novo.pdb
    
    chains_arr=(${subdirpath//_/ })
    # Remove the first element (the tag number)
    chains_arr=("${chains_arr[@]:1}")
    echo "Aligning ${subdirpath} using chains: ${chains_arr[*]}"
    ref_pdb="${dir}/${subdirpath}_protein_cg.pdb"
    #python ${script}/align_structures.py -mobile ${dir}/${subdirpath}/FINAL/final.pdb -ref ${ref_pdb} -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o ${dir}/${subdirpath}/FINAL/final_aligned.pdb > ${dir}/${subdirpath}_align.log 2>&1 &
    python "${an1}/align_structures.py" -mobile "${dir}/${subdirpath}/FINAL/final.pdb" -ref "${ref_pdb}" -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o "${dir}/${subdirpath}/FINAL/final_aligned.pdb"     
  done
}

function lifetimes_to_cif_psii(){
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

function change_lhcbm_chain_lifetime_binding_mode(){
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_modified_chain_ids
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  mkdir -p ${odir}
  rm -f ${odir}/*


}

function add_lifetimes_to_cif(){
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_psbs_summary
  lifetimes_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis

  chains=("s" "n" "8" "7" "k" "z")
  wdir=/martini/rubiz/Github/PsbS_Binding_Site
  dir3=${wdir}/3_reference_proteins
  pdb0=${dir3}/chains_aa/5XNL_chains.pdb
  pdb_psbs=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_4_0_dimer_aligned.pdb

  for chain in "${chains[@]}"; do
    echo "Processing chain: ${chain}"
    python  ${script}/add_lifetimes_to_cif.py --log_transform  -f ${pdb0} -sel "chainID ${chain}" \
      -csv ${lifetimes_dir}/lifetimes_summary_df_chain_${chain}.csv -o ${odir}/sum_chain_${chain}.cif
    python  ${script}/add_lifetimes_to_cif.py --log_transform  -f ${pdb_psbs} -sel "all" \
      -csv ${lifetimes_dir}/lifetimes_summary_df_psbs_chain_${chain}_symmetrized.csv -o ${odir}/sum_psbs_chain_${chain}.cif
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

function plot_aligned_sequences(){
  script=${root_dir}/4_pairs/analysis
  edge_color='None'
  edge_linewidth=0.3
  vmin=100
  vmax=30000
  cmap='colorcet:CET_L17' # Uses the cmap library
  dir1="${root_dir}/5_psii/binding_sites/lifetimes_summary"

  split_every=82
  fasta="${root_dir}/fasta/5xnl_lhc_set_aligned_paper.fasta"
  output="${root_dir}/4_pairs/analysis/figures/aligned_sequences_chains_psii.eps"

  file1=${dir1}/lifetimes_summary_df_chain_s.csv
  file2=${dir1}/lifetimes_summary_df_chain_8.csv
  file3=${dir1}/lifetimes_summary_df_chain_7.csv
  file4=${dir1}/lifetimes_summary_df_chain_n.csv

  headers=("CP26" "CP24" "LHCB3" "LHCBM")
  files=($file1 $file2 $file3 $file4)
  python ${script}/plot_aligned_sequences.py --log_transform --edge_color ${edge_color} --edge_linewidth ${edge_linewidth} --fasta ${fasta} --files ${files[@]} --headers ${headers[@]} --output ${output} --vmin ${vmin} --vmax ${vmax} --cmap ${cmap} --split_every ${split_every}

  split_every=110
  cmap='colorcet:CET_L17' # Uses the cmap library
  fasta="${root_dir}/fasta/psbs.fasta"
  output="${root_dir}/4_pairs/analysis/figures/aligned_sequences_psbs_psii.eps"

  file1=${dir1}/lifetimes_summary_df_psbs_chain_s_symmetrized.csv
  file2=${dir1}/lifetimes_summary_df_psbs_chain_8_symmetrized.csv
  file3=${dir1}/lifetimes_summary_df_psbs_chain_7_symmetrized.csv
  file4=${dir1}/lifetimes_summary_df_psbs_chain_n_symmetrized.csv

  headers=("PsbS_CP26" "PsbS_CP24" "PsbS_LHCB3" "PsbS_LHCBM")
  files=($file1 $file2 $file3 $file4)

  python ${script}/plot_aligned_sequences.py --log_transform --edge_color ${edge_color} --edge_linewidth ${edge_linewidth} --max_resids 212 --fasta ${fasta} --files ${files[@]} --headers ${headers[@]} --output ${output} --vmin ${vmin} --vmax ${vmax} --cmap ${cmap} --split_every ${split_every}

}

function plot_cofactors_lifetimes(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes_summary
  output_dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/PSII_LHCII/psii_with_cofactors_aa.pdb
  equivalent_chains_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/equivalent_chainids.yaml
  change_resnames_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/change_resnames.yaml
  helix_labels_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix_labels.yaml
  vmin=100
  vmax=30000
  cmap='colorcet:CET_L17' # Uses the cmap library
  mkdir -p ${output_dir}
  python3 ${script}/plot_cofactors_lifetimes.py -log_transform -ref ${ref_pdb} -csv ${dir}/lifetimes_summary_df_cofactors_all.csv -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} -output ${output_dir}/cofactors_lifetimes_psii.eps -change_resnames_yaml ${change_resnames_yaml} 
}

function main(){
  set -e  

  #lifetime_analysis_protein_protein  # Get the binding events a csv file.

  #extract_binding                    # Extract binding events (pdb, xtc, tpr)
  #write_equivalent_binding_sites     # Group binding sites


  #check_selections
  #lifetime_analysis_grouped          # Calculate contacts for each subtrajectory
  #sum_csv_lifetimes                  # Sum per chain 
  #symmetric_sum_lifetimes_psbs_chains
  #get_lifetimes_per_binding_mode      # Get binding time per binding mode
  #align_trajectories             
  
  #binding_pose_grouped               # Clustering analysis. 
  #extract_cluster                    # Extract middle structure from largest cluster as gmx cluster generates corrupted PDBs
  plot_psii_binding_modes
  #plot_psii_venn_diagram
  #check_sucess_cg2at
  #cg2at 
  #check_sucess_cg2at
  #reassign_chains 
  #lifetimes_to_cif_psii              # CIF files allow bfactors > 999 while PDB files do not.
  #add_lifetimes_to_cif
  #lifetimes_statistics_psii          # Max occupancy
  #plot_lifetimes                     # Generate protein sequence plots with B-factor coloring
  #write_databases
  #join_databases
  #plot_aligned_sequences
  #plot_cofactors_lifetimes
}

main
