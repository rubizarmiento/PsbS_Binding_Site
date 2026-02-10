root_dir=/martini/rubiz/Github/PsbS_Binding_Site

analysis_dir=${root_dir}/4_pairs/analysis/analysis_pairs
scripts_dir=${root_dir}/4_pairs/analysis
chains=("4" "c" "r" "s")
chains_analyze=("4" "r" "s") # chain c has no binding events longer than 1000 ns
cg2at_path=/martini/rubiz/thylakoid/scripts/para/bin/cg2at
mdp_tpr=${root_dir}/4_pairs/mdps/em.mdp # Dummy mdp to generate tpr files

function check_selections(){
  for chain in "${chains[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir1=${wdir}/1_trj
    pdb=${idir1}/initial_fit_merged.pdb
    sel="chainID ${chain} and (not resname *GG* *DG* *SQ* *PG* *MG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR HEM*)" # Only chlorophylls and proteins
    python3 ${scripts_dir}/return_non_protein_residues.py -f ${pdb} -sel "${sel}" 
  done
}

function lifetimes_classification(){
  for chain in "${chains[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    odir=${wdir}/2_lifetimes_classification
    idir1=${wdir}/1_trj
    pdb=${idir1}/initial_fit_merged.pdb
    #xtc=${idir}/test_ultrashort.xtc
    xtc=${idir1}/proteins_5000ns_concat.xtc
    sel1="chainID A B"
    sel2="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR LMG DGD)" # Only chlorophylls and proteins
    cutoff=8
    dt=2 # time step between frames
    min_event_ns=1000
    sub_trj_nframes=2475 # 5000 ns, 2 ns per frame, minus 25 frames discarded as equilibration
    echo "Starting contact analysis for chain ${chain}..."
    echo "LOG: ${odir}/lifetime_protein.log"
    mkdir -p ${odir}  
    rm -f ${odir}/*  
    # -split every n_frames ensure the analysis is done indivially for each sub-trajectory
    python3 ${scripts_dir}/lifetime_analysis.py -split_every_n_frames ${sub_trj_nframes} -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${pdb} -traj ${xtc} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir} -prefix lifetime_protein -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/2_lifetimes_classification.log 2>&1 &
  done
}

function extract_binding(){
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir1=${wdir}/1_trj
    idir2=${wdir}/2_lifetimes_classification
    odir=${wdir}/3_trj_grouped
    
    df=${idir2}/lifetime_protein_events_df.csv
    f=${idir1}/initial_fit_merged.pdb
    trj=${idir1}/aligned_5000ns.xtc
    sel="all"

    mkdir -p ${odir}
    echo "Extracting binding events for chain ${chain}..."
    echo "LOG: ${odir}/extract_binding_chain_${chain}.log"
    python3 ${scripts_dir}/extract_trajectories_from_dataframe.py -sort "lifetime_ns" -sel ${sel} -i ${df} -f ${f} -trj ${trj} -o ${odir} -prefix chain_${chain} > ${odir}/extract_binding_chain_${chain}.log 2>&1 &
  done
}

function binding_pose_grouped(){
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir1=${wdir}/1_trj
    idir3=${wdir}/3_trj_grouped
    odir=${wdir}/4_trj_cluster

    trj_arr=($(ls ${idir3}/chain_${chain}_*.xtc))
    f=${idir1}/initial_fit_merged.pdb
    tpr=${idir1}/protein.tpr
    ocsv=${odir}/basenames_binding.csv
    
    sel1="chainID A B and name BB"  # PsbS backbone
    sel2="(chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR LMG DGD)) and name BB" # Protein backbone
    osel="all"
    cutoff=0.45

    # 'counter' is used to generate unique output file names for each trajectory in the loop.
    counter=0
    mkdir -p ${odir}
    rm -rf ${odir}/*
          
    echo "Clustering binding poses for chain ${chain}"
    echo "LOG: ${odir}/binding_pose_grouped.log"
    for trj in "${trj_arr[@]}"; do
      counter=$((counter + 1))
      echo chain_${chain}_${counter} >> ${odir}/binding_basenames.csv
      python3 ${scripts_dir}/binding_pose.py -f ${f} -trj ${trj} -tpr ${tpr} -sel1 "${sel2}" -sel2 "${sel1}"  -osel "${osel}" -o ${odir}/chain_${chain}_${counter} --cutoff ${cutoff} >> ${odir}/binding_pose_grouped.log 2>&1 &
    done
  done
}

function extract_cluster(){
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
  
    idir1=${wdir}/1_trj
    idir3=${wdir}/3_trj_grouped
    idir4=${wdir}/4_trj_cluster
    odir=${wdir}/5_middle_cluster
    
    csv=${idir4}/binding_basenames.csv
    f=${idir1}/initial_fit_merged.pdb
    mapfile -t basenames < "${csv}"

    mkdir -p ${odir}
    rm -rf ${odir}/*
    for basename in "${basenames[@]}"; do
      trj=${idir3}/${basename}*.xtc
      log=${idir4}/${basename}/clust_c045/cluster.log
      python3 ${scripts_dir}/extract_cluster.py -f ${f} -trj ${trj} -g ${log} -o ${odir}/${basename}.pdb -n 5 # Extract 5 frames before the middle cluster time, useful when cg2at fails
    done
  done
}

function extract_cluster_patch(){
  # Extract only the first cluster member for chain_4_3 case
  chain="4"
  basename="chain_4_3"
  
  wdir=${analysis_dir}/chain_${chain}
  idir1=${wdir}/1_trj
  idir3=${wdir}/3_trj_grouped
  idir4=${wdir}/4_trj_cluster
  odir=${wdir}/5_middle_cluster
  
  f=${idir1}/initial_fit_merged.pdb
  trj=${idir3}/${basename}*.xtc
  log=${idir4}/${basename}/clust_c045/cluster.log
  
  mkdir -p ${odir}
  
  echo "Extracting first cluster member for ${basename}..."
  python3 ${scripts_dir}/extract_cluster.py -f ${f} -trj ${trj} -g ${log} -o ${odir}/${basename}.pdb --first-member # Extract the first cluster member instead of the middle cluster
  rm -f /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/analysis_pairs/chain_4/6_cg2at/chain_4_3/FINAL/final_cg2at_de_novo.pdb
}

function tpr_cofactors(){
  for chain in "${chains[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir1=${wdir}/1_trj
    odir=${wdir}/1_trj
    f=${idir1}/initial_fit_merged.pdb
    top=${idir1}/no_nb.top

    # Copy .top and comment the lines that start with chain and psbs
    cp ${top} ${odir}/cofactors.top
    non_cofactors=("chain" "psbs" "W" "NA" "CL " "YFPF" "YFPG" "YPPG" "FPGG" "DFGG" "FPMG" "DFMG" "FPSQ" "DPPG" "DSGG" "DSMG" "DPSQ")
    for res in "${non_cofactors[@]}"; do
      sed -i "/^${res}/s/^/;/" ${odir}/cofactors.top
    done

    new_path=Github/PsbS_Binding_Site
    sed -i "s|interaction_partner_paper|${new_path}|g" ${odir}/cofactors.top
    sel="resname PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB HEME CHL"
    python3 ${scripts_dir}/sel_to_ndx.py -f ${f} -sel "${sel}" -name "Cofactors" -o ${odir}/cofactors.ndx
    echo "Cofactors\n" | gmx editconf -f ${f} -n ${odir}/cofactors.ndx -o ${odir}/cofactors.pdb
    gmx grompp -f ${mdp_tpr} -c ${odir}/cofactors.pdb -p ${odir}/cofactors.top -o ${odir}/cofactors.tpr -maxwarn 100
  done
}

function cg2at(){
  rewrite=${1:-"false"}  # Default to false if not provided
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir5=${wdir}/5_middle_cluster
    dir1=${wdir}/1_trj 

    ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/base_dir/rotated.pdb
    files=("$idir5"/*.pdb)

    csv=${wdir}/4_trj_cluster/binding_basenames.csv
    mapfile -t basenames < "${csv}"

    for basename in "${basenames[@]}"; do

      odir=${wdir}/6_cg2at/${basename}/
      mkdir -p ${odir}
      if [ "$rewrite" == "true" ]; then
        rm -rf ${odir}/*
      fi

      o=${odir}/FINAL/final_cg2at_de_novo.pdb
        
      # Skip if output exists
      if [ -f "${o}" ]; then
        echo "Skipping ${o}, output already exists."
        continue
      else
        echo "Processing ${file}..."
        rm -rf ${odir}/*  
        sel="not resname *GG* *DG* *SQ* *PG* *MG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR HEM* CLA CHL CLB"
        cofactors="resname PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB HEME CHL"
        file=${idir5}/${basename}.pdb
        #align to reference pdb
        python3 ${scripts_dir}/align_structures.py -mobile ${file} -ref ${ref_pdb} -sel_ref "name BB and chainID ${chain}" -sel_mobile "name BB and chainID ${chain}" -o ${odir}/${basename}_aligned.pdb

        python3 ${scripts_dir}/sel_to_ndx.py -f ${odir}/${basename}_aligned.pdb -sel "${sel}" -name "Protein" -o ${odir}/${basename}_protein_cg.ndx
        python3 ${scripts_dir}/sel_to_ndx.py -f ${odir}/${basename}_aligned.pdb -sel "${cofactors}" -name "Cofactors" -o ${odir}/${basename}_cofactors_cg.ndx
        gmx editconf -f ${odir}/${basename}_aligned.pdb -n ${odir}/${basename}_protein_cg.ndx -o ${odir}/${basename}_protein_cg.pdb
        
        # Check if ndx is empty
        if [ -s ${odir}/${basename}_cofactors_cg.ndx ]; then
          gmx editconf -f ${odir}/${basename}_aligned.pdb -n ${odir}/${basename}_cofactors_cg.ndx -o ${odir}/${basename}_cofactors_cg.pdb
          python3 ${scripts_dir}/assing_resid_chain_from_pdb.py -o ${odir}/${basename}_cofactors_cg.pdb  -ref ${odir}/${basename}_cofactors_cg.pdb -i ${odir}/${basename}_cofactors_cg.pdb
        fi
        #cg2at - use relative path to avoid path duplication issues
        (cd ${odir} && ${cg2at_path} -c ${basename}_protein_cg.pdb -ff charmm36-jul2020-updated -fg martini_3-0_charmm36 -w tip3p -loc . >> cg2at.log 2>&1) &
      fi
    done
  done
}

function lifetime_analysis_grouped (){
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir1=${wdir}/1_trj
    idir4=${wdir}/4_trj_cluster
    idir3=${wdir}/3_trj_grouped
    odir=${wdir}/7_lifetimes_grouped
    
    csv=${idir4}/binding_basenames.csv
    f=${idir1}/initial_fit_merged.pdb
    mapfile -t basenames < "${csv}"
    tpr=${idir1}/protein.tpr
    
    mkdir -p ${odir}
    #rm -f ${odir}/*    
    for basename in "${basenames[@]}"; do

      trj=${idir3}/${basename}*.xtc

      chains_part=$(echo ${basename} | cut -d'_' -f2-)
      chains_arr=(${chains_part//_/ })
      sel1="chainID A" # PsbS
      sel2="resname PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB HEME CHL" # Only chlorophylls, HEME, and carotenoids
      cutoff=8
      dt=2 # time step between frames
      min_event_ns=100


      #python3 ${scripts_dir}/lifetime_analysis.py -prefix test -n_frames 100 -dt ${dt} -min_event_ns 0 -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix test_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
      #Contacts PsbS and chains
      #python3 ${scripts_dir}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix psbs_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_${basename}.log 2>&1 &
      for chain in "${chains_arr[@]}"; do
        # Contacts chains and PsbS
        sel3="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* LMG DGD)" # Only chlorophylls, carotenoids and proteins
        sel4="chainID ${chain} and (resname CLA CLB CHL, PLQ PL9 LUT VIO XAT NEO NEX BCR HEME)" # Only chlorophylls and carotenoids 

        #python3 ${scripts_dir}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel3}" -o ${odir} -prefix psbs_chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/psbs_chain_${chain}_${basename}.log 2>&1 &
        #python3 ${scripts_dir}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir} -prefix chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/chain_${chain}_${basename}.log 2>&1 &
        python3 ${scripts_dir}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${f} -traj ${trj} -sel1 "${sel4}" -sel2 "${sel1}" -o ${odir} -prefix cofactors_chain_${chain}_${basename} -group_by1 "resids" -group_by2 "resids" > ${odir}/chlorophylls_chain_${chain}_${basename}.log 2>&1 &
      done
    done
  done
}



function reassign_chains(){
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    odir6=${wdir}/6_cg2at
    idir5=${wdir}/5_middle_cluster
    idir4=${wdir}/4_trj_cluster
    csv=${idir4}/binding_basenames.csv
    mapfile -t basenames < "${csv}"
    files=("$idir5"/*.pdb)
  
    for basename in "${basenames[@]}"; do
      odir=${wdir}/6_cg2at/${basename}/

      o=${odir}/FINAL/final_cg2at_de_novo.pdb

      python3 ${scripts_dir}/assing_resid_chain_from_pdb.py -o ${odir}/FINAL/final.pdb  -ref ${odir}/${basename}_protein_cg.pdb -i ${odir}/FINAL/final_cg2at_de_novo.pdb
      
      chains_arr=(${basename//_/ })
      # Remove the first element (the tag number) and the last
      chains_arr=("${chains_arr[@]:1}")
      chains_arr=("${chains_arr[@]::${#chains_arr[@]}-1}")
      echo "Aligning ${basename} using chains: ${chains_arr[*]}"
      ref_pdb="${odir}/${basename}_protein_cg.pdb"
      #python ${script}/align_structures.py -mobile ${dir}/${subdirpath}/FINAL/final.pdb -ref ${ref_pdb} -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o ${dir}/${subdirpath}/FINAL/final_aligned.pdb > ${dir}/${subdirpath}_align.log 2>&1 &
      python "${scripts_dir}/align_structures.py" -mobile "${odir}/FINAL/final.pdb" -ref "${ref_pdb}" -sel_ref "name BB and chainID ${chains_arr[*]}" -sel_mobile "name CA and chainID ${chains_arr[*]}" -o "${odir}/FINAL/final_aligned.pdb"
    done
  done
}

function lifetimes_to_cif(){
  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    idir7=${wdir}/7_lifetimes_grouped
    idir4=${wdir}/4_trj_cluster
    csv=${idir4}/binding_basenames.csv
    odir=${wdir}/8_cifs_lifetimes
    pdb_dir=${wdir}/6_cg2at

    mapfile -t basenames < "${csv}"
    
    # Copy cofactors files from subdirectories to pdb_dir for Python script to find them
    for basename in "${basenames[@]}"; do
      cofactors_src=${pdb_dir}/${basename}/${basename}_cofactors_cg.pdb
      cofactors_dst=${pdb_dir}/${basename}_cofactors_cg.pdb
      if [ -f "${cofactors_src}" ]; then
        cp "${cofactors_src}" "${cofactors_dst}"
      fi
    done
    
    sel_protein="(chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR LMG DGD CHL CLA CLB HEM*))"
    sel_cofactors="resname PLQ PL9 LUT VIO XAT NEO NEX BCR CLA CLB HEME CHL"
    sel_psbs="chainID A"
    rm -rf ${odir}/*pdb
    python3 ${scripts_dir}/lifetime_to_cif.py ${pdb_dir} ${idir7} ${odir} "${sel_protein}" "${sel_cofactors}" "${sel_psbs}" ${csv}
  done
}

function sum_csv_lifetimes(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  odir=${analysis_dir}/10_lifetimes_summary

  mkdir -p ${odir}
  rm -f ${odir}/*
  
  chains=("4" "r" "s")
  outputs_arr_psbs=()
  outputs_arr_chains=()
  outputs_arr_cofactors=()
  for chain in "${chains[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    dir=${wdir}/7_lifetimes_grouped
    echo "Processing chain: ${chain}"
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_chain_${chain}.csv -prefix chain_${chain}_ # Suffix is *residue_summary_df.csv
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_cofactors_chain_${chain}.csv -prefix cofactors_chain_${chain}_  --ignore-missing  # Suffix is *residue_summary_df.csv, ignores empty files for chains without cofactors
    python3 ${script}/sum_csv_lifetimes.py -d ${dir} -o ${odir}/lifetimes_summary_df_psbs_chain_${chain}.csv -prefix psbs_chain_${chain}_ # Suffix is *residue_summary_df.csv
    outputs_arr_psbs+=("${odir}/lifetimes_summary_df_psbs_chain_${chain}.csv")
    outputs_arr_chains+=("${odir}/lifetimes_summary_df_chain_${chain}.csv")
    # Check if the cofactors file is not empty before adding to the array
    if [ -s "${odir}/lifetimes_summary_df_cofactors_chain_${chain}.csv" ]; then
      outputs_arr_cofactors+=("${odir}/lifetimes_summary_df_cofactors_chain_${chain}.csv")
    fi
  done
  python3 ${script}/join_csvs.py -c "${outputs_arr_psbs[@]}" -o ${odir}/lifetimes_summary_df_psbs_all.csv
  python3 ${script}/join_csvs.py -c "${outputs_arr_chains[@]}" -o ${odir}/lifetimes_summary_df_chains_all.csv
  python3 ${script}/join_csvs.py -c "${outputs_arr_cofactors[@]}" -o ${odir}/lifetimes_summary_df_cofactors_all.csv
}

function sum_lifetimes_psbs_chains(){
  odir=${analysis_dir}/10_lifetimes_summary
  idir=${analysis_dir}/10_lifetimes_summary
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis

  chains=("4" "r" "s")

  for chain in "${chains[@]}"; do
    python  ${script}/sum_lifetimes_psbs_chains.py -csv ${idir}/lifetimes_summary_df_psbs_chain_${chain}.csv -o ${odir}/lifetimes_summary_df_psbs_chain_${chain}_symmetrized.csv
  done
}

function add_lifetimes_to_cif(){
  odir=${analysis_dir}/11_cifs_psbs_summary
  lifetimes_dir=${analysis_dir}/10_lifetimes_summary
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis



  chains=("4" "r" "s")
  wdir=/martini/rubiz/Github/PsbS_Binding_Site
  dir3=${wdir}/3_reference_proteins
  pdb0=${dir3}/chains_aa/5XNL_chains.pdb
  pdb_psbs=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/psbs/psbs_4_0_dimer_aligned.pdb

  for chain in "${chains[@]}"; do
    python  ${script}/add_lifetimes_to_cif.py --log_transform  -f ${pdb0} -sel "chainID ${chain}" \
      -csv ${lifetimes_dir}/lifetimes_summary_df_chain_${chain}.csv -o ${odir}/sum_chain_${chain}.cif
    python  ${script}/add_lifetimes_to_cif.py --log_transform  -f ${pdb_psbs} -sel "all" \
      -csv ${lifetimes_dir}/lifetimes_summary_df_psbs_chain_${chain}_symmetrized.csv -o ${odir}/sum_psbs_chain_${chain}.cif
  done
}

function plot_lifetimes(){
  #Same as PSII analysis but for each pair chain
  
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  
  # Share with PSII analysis, plot costummization
  chain_labels_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/chain_labels.yaml
  color_config_yaml=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/psii_psbs/color_definitions.yaml

  # Output YAML files defining helices
  psii_helix_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix.yaml
  psbs_helix_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psbs_helix.yaml
  # PDB with helices defined (for reference)
  psii_pdbdatabase=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/5XNL.pdb
  psbs_pdbdatabase=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/4ri2.pdb

  # Generate helix YAML files from PDB (helix definitions)
  python3 ${script}/write_helix_yaml.py -o ${psii_helix_yaml} -f ${psii_pdbdatabase}
  python3 ${script}/write_helix_yaml.py -o ${psbs_helix_yaml} -f ${psbs_pdbdatabase}


  for chain in "${chains_analyze[@]}"; do
    wdir=${analysis_dir}/chain_${chain}
    cifs_dir=${wdir}/8_cifs_lifetimes
    basenames_csv=${wdir}/4_trj_cluster/binding_basenames.csv
    output_dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/analysis_pairs/chain_${chain}/9_lifetimes_sequences

    mkdir -p ${output_dir}

    # Generate protein sequence plots with B-factor coloring and helix annotations
    echo "Generating protein sequence plots with lifetimes visualization for chain ${chain}..."
    python3 ${script}/plot_lifetimes_sequences.py \
      -d ${cifs_dir} \
      -b ${basenames_csv} \
      -l ${chain_labels_yaml} \
      -c ${color_config_yaml} \
      -p ${psii_helix_yaml} \
      -s ${psbs_helix_yaml} \
      -o ${output_dir} \
      --split-sequences 106

    echo "Sequence plots for chain ${chain} saved to: ${output_dir}"
  done
}


function plot_aligned_sequences(){
  script=${root_dir}/4_pairs/analysis

  vmin=100
  vmax=30000
  cmap='colorcet:CET_L17' # Uses the cmap library
  split_every=82
  edge_color='None'
  edge_linewidth=0.3
  col_width=0.1
  row_height=0.4
  dir1=${analysis_dir}/10_lifetimes_summary



  fasta="${root_dir}/fasta/5xnl_lhc_set_aligned_paper.fasta"
  output="${root_dir}/4_pairs/analysis/figures/aligned_sequences_chains_pairs.eps"
  

  file1=${dir1}/lifetimes_summary_df_chain_s.csv
  file2=${dir1}/lifetimes_summary_df_chain_4.csv
  file3=${dir1}/lifetimes_summary_df_chain_r.csv

  headers=("CP26" "CP24" "CP29")
  files=($file1 $file2 $file3)

  python ${script}/plot_aligned_sequences.py --log_transform --col_width ${col_width} --row_height ${row_height} --edge_color ${edge_color} --edge_linewidth ${edge_linewidth} --fasta ${fasta} --files ${files[@]} --headers ${headers[@]} --output ${output} --vmin ${vmin} --vmax ${vmax} --cmap ${cmap} --split_every ${split_every}


  fasta="${root_dir}/fasta/psbs.fasta"
  output="${root_dir}/4_pairs/analysis/figures/aligned_sequences_psbs_pairs.eps"

  file1=${dir1}/lifetimes_summary_df_psbs_chain_s_symmetrized.csv
  file2=${dir1}/lifetimes_summary_df_psbs_chain_4_symmetrized.csv
  file3=${dir1}/lifetimes_summary_df_psbs_chain_r_symmetrized.csv
  split_every=110
  cmap='colorcet:CET_L17' # Uses the cmap library

  headers=("PsbS_CP26" "PsbS_CP24" "PsbS_CP29")
  files=($file1 $file2 $file3)

  python ${script}/plot_aligned_sequences.py --log_transform --col_width ${col_width} --row_height ${row_height} --edge_color ${edge_color} --edge_linewidth ${edge_linewidth} --max_resids 212 --fasta ${fasta} --files ${files[@]} --headers ${headers[@]} --output ${output} --vmin ${vmin} --vmax ${vmax} --cmap ${cmap} --split_every ${split_every}
}

function plot_cofactors_lifetimes(){
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  dir=${analysis_dir}/10_lifetimes_summary
  output_dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/figures
  ref_pdb=/martini/rubiz/Github/PsbS_Binding_Site/3_reference_proteins/PSII_LHCII/psii_with_cofactors_aa.pdb
  equivalent_chains_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/equivalent_chainids.yaml
  change_resnames_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/change_resnames.yaml
  change_chains_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/change_chains.yaml
  helix_labels_yaml=/martini/rubiz/Github/PsbS_Binding_Site/definitions_yaml/psii_helix_labels.yaml
  vmin=100
  vmax=30000
  cmap='colorcet:CET_L17' # Uses the cmap library
  mkdir -p ${output_dir}
  python3 ${script}/plot_cofactors_lifetimes.py -log_transform -ref ${ref_pdb} -csv ${dir}/lifetimes_summary_df_cofactors_all.csv -vmin ${vmin} -vmax ${vmax} -cmap ${cmap} -output ${output_dir}/cofactors_lifetimes_pairs.eps -change_resnames_yaml ${change_resnames_yaml} -change_chains_yaml ${change_chains_yaml} 

}



function main(){
  set -e
  #check_selections
  #lifetimes_classification

  
  #extract_binding
  #binding_pose_grouped   
  #lifetime_analysis_grouped

  #extract_cluster
  #extract_cluster_special          # Backmapping failed for chain_4_3, extract first cluster member instead
  #tpr_cofactors

  #cg2at 
  #reassign_chains

    
  #lifetimes_to_cif                   # CIF files allow bfactors > 999 while PDB files do not.
  #sum_csv_lifetimes
  #sum_lifetimes_psbs_chains
  #add_lifetimes_to_cif
  
  #lifetimes_statistics_psii          # Max occupancy
  #plot_lifetimes
  #plot_aligned_sequences
  plot_cofactors_lifetimes
}

main
