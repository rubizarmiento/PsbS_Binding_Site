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
    sel2="chainID ${chain} and not resname CLA CLB CHL *HG* PLQ PL9 *GG* *SQ* *PG* LUT VIO XAT NEO NEX W2 HOH BCR"
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
    min_event_ns=100
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
    sel3="chainID ${chain} and not resname *GG* *SQ* *PG* W* HOH *HG*" # Only carotenes, chlorophylls and proteins
    cutoff=8
    dt=2 # time step between frames
    min_event_ns=1000
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel3}" -sel2 "${sel1}" -o ${odir}/lifetime -prefix lifetime_protein_chain_${chain} -group_by1 "chainIDs" -group_by2 "chainIDs" > ${odir}/lifetime/lifetime_protein_chain_${chain}.log 2>&1 &
  done
}

function save_top_lifetimes_old(){
  mkdir -p ${an1}/top_lifetimes_pdbs
  cd ${an1}/top_lifetimes_pdbs
  sel="not resname W2 DPPG"
  for chain in "${chains[@]}"; do
    df=${odir}/lifetime/lifetime_protein_chain_${chain}_events_df.csv
    f=${an1}/chain_${chain}/initial_fit.pdb
    xtc=${an1}/chain_${chain}/aligned_5000ns.xtc
    python3 ${an1}/extract_trajectories_from_dataframe.py -join -sel ${sel}-i ${df} -f ${f} -trj ${xtc} -sort lifetime_ns -n_trj 5 -o . -prefix chain_${chain} > chain_${chain}.log 2>&1 &
  done
}

function binding_pose_pdb(){
  mkdir -p ${an1}/top_lifetimes_pdbs/binding_poses_pdbs
  cd ${an1}/top_lifetimes_pdbs/binding_poses_pdbs

  for chain in "${chains[@]}"; do
    f=${an1}/chain_${chain}/initial_fit.pdb
    trj=${an1}/chain_${chain}/aligned_5000ns.xtc
    #The trajectories are named: chain_${chain}_${n_top_liftime}_${start_frame}_${end_frame}.xtc
    # Array trajectories with the *${chain}*xtc
    python ${an1}/binding_pose.py -f ${f} -trj ${trj} -tpr ${an1}/chain_${chain}/protein.tpr -sel1 "${sel1}" -sel2 "${sel2}" -o binding_poses_chain_${chain}_${counter} -filter 0.8 --cutoff 0.5 0.75 1 1.5 2 >> chain_${chain}_binding.log 2>&1 &
  done
}

function binding_pose_pdb(){
  mkdir -p ${an1}/binding_poses
  cd ${an1}/binding_poses
  gro=initial_fit_merged.pdb
  #xtc=test_ultrashort.xtc
  xtc=aligned_5000ns.xtc
  sel1="chainID A"
  for chain in "${chains[@]}"; do
    sel2="chainID ${chain} and not resname *GG* *SQ* *PG* W* HOH *HG* PLQ PL9 LUT VIO XAT NEO NEX BCR" # Only chlorophylls and proteins
    f=${an1}/chain_${chain}/initial_fit.pdb
    trj=${an1}/chain_${chain}/aligned_5000ns.xtc
    timestamp=$(date +"%Y%m%d_%H%M%S")
    python ${an1}/binding_pose.py -f ${f} -trj ${trj} -tpr ${an1}/chain_${chain}/protein.tpr -sel1 "${sel1}" -sel2 "${sel2}" -o binding_poses_chain_${chain} -filter 0.8 --cutoff 0.5 0.75 1 1.5 2 >> chain_${chain}_binding_${timestamp}.log 2>&1 &
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
  #save_top_lifetimes
  binding_pose_pdb
}

main
