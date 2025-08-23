# Important Scripts

## PsbS and its Potential Partner Proteins

- **4_pairs/analysis/merge_trajectories.ipynb**: Merge simulations with complex names generated after appending was not possible.
- **4_pairs/analysis/1_preprocessing_copy1.sh**: Concatenates and aligns simulations
- **4_pairs/analysis/contact_analysis.ipynb**: Contact and lifetime analysis.
- **4_pairs/analysis/analysis.sh**: Binding poses with clustering analysis.

## PsbS and the PSII-LHCII Supercomplex

- **5_psii/protocol_psii_psbs.sh**: Generate initial structures PSII-PsbS
- **5_psii/psii_psbs/run_all.sh**: Simulation steps
- **5_psii/psii_psbs/protocol_snellius.sh**: Function to run simulation in HPC

# Important Files

## PsbS and its Potential Partner Proteins

- `4_pairs/analysis/chain_${chain}/initial.pdb`
- `4_pairs/analysis/chain_${chain}/protein.tpr`
- `4_pairs/analysis/chain_${chain}/aligned_5000ns.xtc`