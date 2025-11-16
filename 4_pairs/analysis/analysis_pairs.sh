an1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
odir=${an1}/csv_files
chains=("4" "c" "r" "s")

function contact_analysis_protein(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit.pdb
    #xtc=test_ultrashort.xtc
    xtc=aligned_5000ns.xtc
    sel1="chainID A B"
    sel2="chainID ${chain} +and not resname CLA CLB CHL *HG* PLQ PL9 *GG* *SQ* *PG* LUT VIO XAT NEO NEX W2 HOH BCR"
    cutoff=8
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/contact_analysis.py -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_protein_chain${chain}_5000ns.csv -group_by1 "resids" -group_by2 "resids" > ${odir}/lifetime/contact_matrix_protein_chain${chain}_5000ns.log 2>&1 &
  done
}

function contact_analysis_cofactors(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit.pdb
    #xtc=test_ultrashort.xtc
    xtc=aligned_5000ns.xtc
    sel1="chainID A B"
    sel2="resname CLA CLB CHL PLQ PL9 LUT VIO XAT NEO NEX BCR"
    cutoff=8
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/contact_analysis.py -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_cofactors_chain${chain}_5000ns.csv -group_by1 "resids" -group_by2 "resnames" > ${odir}/lifetime/contact_matrix_cofactors_chain${chain}_5000ns.log 2>&1 &
  done
}

function contact_analysis_cofactors_resids(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit.pdb
    #xtc=test_ultrashort.xtc
    xtc=aligned_5000ns.xtc
    sel1="chainID A B"
    sel2="resname CLA CLB CHL PLQ PL9 LUT VIO XAT NEO NEX BCR"
    cutoff=8
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/contact_analysis.py -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_cofactors_resids_chain${chain}_5000ns.csv -group_by1 "resids" -group_by2 "resids" > ${odir}/lifetime/contact_matrix_cofactors_resids_chain${chain}_5000ns.log 2>&1 &
  done
}

function lifetime_analysis_protein(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit.pdb
    #xtc=test_ultrashort.xtc
    xtc=aligned_5000ns.xtc
    sel1="chainID A B"
    sel2="chainID ${chain} and not resname CLA CLB CHL *HG* PLQ PL9 *GG* *SQ* *PG* LUT VIO XAT NEO NEX BCR"
    sel3="chainID ${chain}"
    cutoff=8
    dt=2
    min_event_ns=1000
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel3}" -o ${odir}/lifetime -prefix lifetime_psbs_chain_${chain} -group_by1 "resids" -group_by2 "chainIDs" > ${odir}/lifetime/lifetime_psbs_chain_${chain}.log 2>&1 &
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/lifetime -prefix lifetime_chain${chain} -group_by1 "resids" -group_by2 "chainIDs" > ${odir}/lifetime/lifetime_chain_${chain}.log 2>&1 &
  done
}

function lifetime_analysis_cofactors(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit.pdb
    #xtc=test_ultrashort.xtc
    xtc=aligned_5000ns.xtc
    sel1="chainID A B"
    sel2="resname CLA CLB CHL PLQ PL9 LUT VIO XAT NEO NEX BCR"
    sel3="chainID ${chain}"
    cutoff=8
    dt=2
    min_event_ns=100
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel3}" -o ${odir}/lifetime -prefix lifetime_cofactors_psbs_chain_${chain} -group_by1 "resids" -group_by2 "chainIDs" > ${odir}/lifetime/lifetime_psbs_chain_${chain}.log 2>&1 &
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/lifetime -prefix lifetime_cofactors_chain_${chain} -group_by1 "resids" -group_by2 "chainIDs" > ${odir}/lifetime/lifetime_chain_${chain}.log 2>&1 &
  done
}

function lifetime_analysis_protein_protein(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit_merged.pdb
    #xtc=test_ultrashort.xtc
    xtc=aligned_5000ns.xtc
    sel1="chainID A"
    sel3="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)" # Only chlorophylls and proteins
    cutoff=8
    dt=2 # time step between frames
    min_event_ns=1000
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir}/lifetime -prefix lifetime_protein_chain_${chain} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/lifetime/lifetime_protein_chain_${chain}.log 2>&1 &
  done
}

function save_trj_lifetimes(){
  mkdir -p ${an1}/top_lifetimes_pdbs
  cd ${an1}/top_lifetimes_pdbs
  sel="all"
  for chain in "${chains[@]}"; do
    df=${odir}/lifetime/lifetime_protein_chain_${chain}_events_df.csv
    f=${an1}/chain_${chain}/initial_fit_merged.pdb
    trj=${an1}/chain_${chain}/aligned_5000ns.xtc
    python3 ${an1}/extract_trajectories_from_dataframe.py -sort "lifetime_ns" -sel ${sel} -i ${df} -f ${f} -trj ${trj} -o . -prefix chain_${chain} > chain_${chain}.log 2>&1 &
  done
}

function binding_pose_pdb(){
  chain=${1}
  odir=${an1}/binding_poses_pairs
  
  
  mkdir -p ${odir}
  cd ${odir}
  mkdir -p chain_${chain}
  
  cd chain_${chain}
  # the trajectories are named chain${chain}_${id}_${start}_${end}.xtc
  sel1="chainID A B"
  sel2="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR DP* DS*)" # Only chlorophylls and proteins
  f=${an1}/chain_${chain}/initial_fit_merged.pdb
  trj_arr=($(ls ${an1}/top_lifetimes_pdbs/chain_${chain}_*.xtc))
  # 'counter' is used to generate unique output file names for each trajectory in the loop.
  counter=0
  rm -f chain_${chain}_*.log
  for trj in "${trj_arr[@]}"; do
    counter=$((counter + 1))
    python3 ${an1}/binding_pose.py -f ${f} -trj ${trj} -tpr ${an1}/chain_${chain}/protein.tpr -sel1 "${sel2}" -sel2 "${sel1}" -o chain_${chain}_${counter} --cutoff 0.45 >> chain_${chain}_binding.log 2>&1 &
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


function extract_cluster(){
  # Returns: /martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/biggest_clusters_c075 
  script=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
  idir1=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/trj_aligned
  idir2=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/binding_poses_pairs
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/middle_cluster
  trj_arr=($(ls ${an1}/top_lifetimes_pdbs/chain_${chain}_*.xtc))

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
      ${cg2at_path} -c ${basename}_protein_cg.pdb -ff charmm36-jul2020-updated -fg martini_3-0_charmm36 -w tip3p -loc ${basename} >> ${odir}/${basename}.log 2>&1 &
      fi
  done
}

function lifetimes_to_cif_psii(){
  pdb_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cg2at
  lifetimes_dir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/lifetimes
  odir=/martini/rubiz/Github/PsbS_Binding_Site/5_psii/binding_sites/cifs_lifetimes
  sel_protein="not resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR and not chainID 9"
  sel_cofactors="resname CLA CLB CHL *HG* HEM PLQ PL9 *GG* *SQ* *PG* DGD LMG LUT VIO XAT NEO NEX W2 HOH BCR"
  sel_psbs="chainID 9"
  rm -rf ${odir}/*pdb
  python3 /martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/lifetime_to_cif_psii.py ${pdb_dir} ${lifetimes_dir} ${odir} "${sel_protein}" "${sel_cofactors}" "${sel_psbs}"
}


function lifetime_analysis_chain_resname(){
  chains=("4" "r" "s")
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit_merged.pdb
    sel1="chainID A"
    sel2="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)" # Only chlorophylls and proteins
    cutoff=8
    dt=2 # time step between frames
    min_event_ns=0
    dir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs
    n_trj=$(find ${dir} -name "chain_${chain}*.xtc" | wc -l)
    mkdir -p ${dir}/lifetimes_chain_${chain}
    echo "Number of binding events for chain ${chain}: ${n_trj}"
    log_file=${dir}/lifetimes_chain_${chain}/lifetime_analysis_chain_${chain}_resname.log
    rm -f "${log_file}"
    for i in $(seq 1 $((n_trj-1))); do
      trj=$(find ${dir} -name "chain_${chain}_${i}*.xtc")
      if [ -n "$trj" ]; then
        echo "Using trajectory: $trj"
        python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${dir}/lifetimes_chain_${chain} -prefix lifetime_analysis_chain_${chain}_${i}_resname -group_by1 "chainIDs" -group_by2 "resnames" >> "${log_file}" 2>&1 &
      else
        echo "No trajectory found for chain ${chain}, event ${i}" >> "${log_file}"
      fi
    done
  done
}

function lifetime_analysis_classified(){
  chains=("4" "r" "s")
  for chain in "${chains[@]}"; do
    gro=${an1}/chain_${chain}/initial_fit_merged.pdb
    #xtc=test_ultrashort.xtc
    sel1="chainID A" # PsbS
    sel2="chainID ${chain} and (not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR)" # Only chlorophylls and proteins
    cutoff=8
    dt=2 # time step between frames
    min_event_ns=100
    odir=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/top_lifetimes_pdbs/lifetimes_chain_${chain}    
    trj=${an1}/top_lifetimes_pdbs/joined/chain_${chain}_non_chl.xtc
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix chain_${chain}_psbs_non_chl -group_by1 "chainIDs" -group_by2 "resids" > ${odir}/chain_${chain}_psbs_non_chl 2>&1 &
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${trj} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir} -prefix psbs_chain_${chain}_non_chl -group_by1 "chainIDs" -group_by2 "resids" > ${odir}/psbs_chain_${chain}_non_chl 2>&1 &
    trj=${an1}/top_lifetimes_pdbs/joined/chain_${chain}_chl.xtc
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${trj} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir} -prefix chain_${chain}_psbs_chl -group_by1 "chainIDs" -group_by2 "resids" > ${odir}/chain_${chain}_psbs_chl.log 2>&1 &
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${trj} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir} -prefix psbs_chain_${chain}_chl -group_by1 "chainIDs" -group_by2 "resids" > ${odir}/psbs_chain_${chain}_chl.log 2>&1 &
  done
}

function cluster_trajectories(){
  chains=("4" "r" "s")
  for chain in "${chains[@]}"; do
    python3 ${an1}/cluster_trajectories.py ${chain}
  done
}

#python binding_pose.py -f chain_4/initial_fit.pdb -trj chain_4/test.xtc -tpr chain_4/protein.tpr -sel1 "chainID 4" -sel2 "chainID A B" -o binding_poses_main_test -filter 0.8 --cutoff 0.15 0.30 0.45 0.60

function main(){
  #contact_analysis_protein
  #contact_analysis_cofactors
  #contact_analysis_cofactors_resids
  #python3 ${an1}/contacts_to_pdb.py

  #lifetime_analysis_cofactors 
  #lifetime_analysis_protein
  #lifetime_analysis_protein_protein
  

  
  #save_trj_lifetimes
  #Chain c is not being analyzed as it has no binding events longer than 1000 ns.
  #binding_pose_pdb "4"
  #binding_pose_pdb "r"
  #binding_pose_pdb "s"
  #python3 ${an1}/summary_clusters_per_chain.py

  extract_cluster


  #lifetime_analysis_chain_resname # Use to classify interactions between Chl and non-chls
  #lifetime_analysis_classified
  
  #cluster_binding_modes
  #cluster_trajectories
  #lifetime_analysis_filtered                                                                                                                                                        
}

main
