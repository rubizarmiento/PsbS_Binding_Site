
repository_dir=/martini/rubiz/Github/PsbS_Binding_Site/analysis_dataset
scripts_dir=${repository_dir}/scripts

HELIX_DEFINITIONS_YAML_GROUP1=${repository_dir}/definitions_yaml/binding_proteins_helix_labels.yaml
HELIX_DEFINITIONS_YAML_GROUP2=${repository_dir}/definitions_yaml/target_protein_helix_labels_merged.yaml

LABELS_CHAIN_YAML_GROUP1=${repository_dir}/definitions_yaml/chain_labels.yaml
LABELS_CHAIN_YAML_GROUP2=${repository_dir}/definitions_yaml/psbs_labels.yaml

SQL_MODIFIERS=${repository_dir}/sql_operations/01_add_psbs_modifiers.sql
SQL_CLASSIFICATION=${repository_dir}/sql_operations/02_add_residue_classifications.sql


function write_databases_pairs(){
  chains_to_analyze=("A" "B" "C") # chain D has no binding events longer than 1000 ns
    
  for chain in "${chains_to_analyze[@]}"; do
      idir=${repository_dir}/0_datasets/pairs/chain_${chain}
      odir=${repository_dir}/1_concatenated_databases/pairs/chain_${chain}
      
      mkdir -p ${odir}
      python3 ${scripts_dir}/write_databases.py \
          --helix_def_yaml_group1 ${HELIX_DEFINITIONS_YAML_GROUP1} \
          --helix_def_yaml_group2 ${HELIX_DEFINITIONS_YAML_GROUP2} \
          --output ${odir}/database.csv \
          --csv_files "${idir}/*respairs_events*.csv" \
          --add_labels_colname "sim_type" \
          --add_labels_values "pairs" \
          --labels_chain_yaml_group1 ${LABELS_CHAIN_YAML_GROUP1} \
          --labels_chain_yaml_group2 ${LABELS_CHAIN_YAML_GROUP2} 

      python ${scripts_dir}/apply_sql_operations.py \
          --input_csv ${odir}/database.csv \
          --output_csv ${odir}/database.csv \
          --operations ${SQL_MODIFIERS}

      python ${scripts_dir}/apply_sql_operations.py \
          --input_csv ${odir}/database.csv \
          --output_csv ${odir}/database.csv \
          --operations ${SQL_CLASSIFICATION}
  done
}

function write_databases_supercomplex(){
  idir=${repository_dir}/0_datasets/supercomplex
  odir=${repository_dir}/1_concatenated_databases/supercomplex
  mkdir -p ${odir}
  python3 ${scripts_dir}/write_databases.py \
    --helix_def_yaml_group1 ${HELIX_DEFINITIONS_YAML_GROUP1} \
    --helix_def_yaml_group2 ${HELIX_DEFINITIONS_YAML_GROUP2} \
    --output ${odir}/database.csv \
    --csv_files "${idir}/*respairs_events*.csv" \
    --add_labels_colname "sim_type" \
    --add_labels_values "psii_lhcii" \
    --labels_chain_yaml_group1 ${LABELS_CHAIN_YAML_GROUP1} \
    --labels_chain_yaml_group2 ${LABELS_CHAIN_YAML_GROUP2} \
    
    python ${scripts_dir}/apply_sql_operations.py \
      --input_csv ${odir}/database.csv \
      --output_csv ${odir}/database.csv \
      --operations ${SQL_MODIFIERS}

    python ${scripts_dir}/apply_sql_operations.py \
      --input_csv ${odir}/database.csv \
      --output_csv ${odir}/database.csv \
      --operations ${SQL_CLASSIFICATION}
}

function join_databases(){
  file1=${repository_dir}/1_concatenated_databases/pairs/chain_A/database.csv
  file2=${repository_dir}/1_concatenated_databases/pairs/chain_B/database.csv
  file3=${repository_dir}/1_concatenated_databases/pairs/chain_C/database.csv
  file4=${repository_dir}/1_concatenated_databases/supercomplex/database.csv

  odir=${repository_dir}/2_joined_database
  mkdir -p ${odir}
  python3 ${scripts_dir}/join_databases.py \
    --input_files ${file1} ${file2} ${file3} ${file4} \
    --output_file ${odir}/joined_database.csv
}

function main(){
  set -e
  write_databases_pairs
  write_databases_supercomplex
  join_databases
}

main