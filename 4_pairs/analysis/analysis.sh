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
    python3 ${an1}/contact_analysis.py -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_protein_chain${chain}_5000ns.csv -group_by1 "resids" -group_by2 "resids" > ${odir}/contact_matrix_protein_chain${chain}_5000ns.log 2>&1 &
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
    python3 ${an1}/contact_analysis.py -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_cofactors_chain${chain}_5000ns.csv -group_by1 "resids" -group_by2 "resnames" > ${odir}/contact_matrix_cofactors_chain${chain}_5000ns.log 2>&1 &
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
    python3 ${an1}/contact_analysis.py -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_cofactors_resids_chain${chain}_5000ns.csv -group_by1 "resids" -group_by2 "resids" > ${odir}/contact_matrix_cofactors_resids_chain${chain}_5000ns.log 2>&1 &
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
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel3}" -o ${odir}/lifetime -prefix lifetime_cofactors_psbs_chain_${chain} -group_by1 "resids" -group_by2 "chainIDs" > ${odir}/lifetime_psbs_chain_${chain}.log 2>&1 &
    python3 ${an1}/lifetime_analysis.py -dt ${dt} -min_event_ns ${min_event_ns} -cutoff "${cutoff}" -f ${gro} -traj ${xtc} -sel1 "${sel2}" -sel2 "${sel1}" -o ${odir}/lifetime -prefix lifetime_cofactors_chain_${chain} -group_by1 "resids" -group_by2 "chainIDs" > ${odir}/lifetime_chain_${chain}.log 2>&1 &
  done
}


#python binding_pose.py -f chain_4/initial_fit.pdb -trj chain_4/test.xtc -tpr chain_4/protein.tpr -sel1 "chainID 4" -sel2 "chainID A B" -o binding_poses_main_test -filter 0.8 --cutoff 0.15 0.30 0.45

function main(){
  #contact_analysis_protein
  #contact_analysis_cofactors
  #contact_analysis_cofactors_resids
  #python3 ${an1}/contacts_to_pdb.py

  lifetime_analysis_cofactors 
  lifetime_analysis_protein
}

main