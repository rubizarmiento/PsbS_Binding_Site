an1=/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis
odir=${an1}/csv_files
chains=("4" "c" "r" "s")

function contact_analysis_all(){
  for chain in "${chains[@]}"; do
    cd ${an1}/chain_${chain}
    gro=initial_fit.pdb
    xtc=aligned_5000ns.xtc
    sel1="chainID A B"
    sel2="chainID ${chain} and not resname *HG* *GG* *SQ* *PG* W2 HOH"
    cutoff=8
    echo "Starting contact analysis for chain ${chain}..."
    python3 ${an1}/contact_analysis.py -f ${gro} -traj ${xtc} -sel1 "${sel1}" -sel2 "${sel2}" -o ${odir}/contact_matrix_chain_${chain}_all.csv -group_by1 "resids" -group_by2 "resids" > /dev/null 2>&1 &

  done
}



#python binding_pose.py -f chain_4/initial_fit.pdb -trj chain_4/test.xtc -tpr chain_4/protein.tpr -sel1 "chainID 4" -sel2 "chainID A B" -o binding_poses_main_test -filter 0.8 --cutoff 0.15 0.30 0.45

function main(){
  contact_analysis_all
}

main