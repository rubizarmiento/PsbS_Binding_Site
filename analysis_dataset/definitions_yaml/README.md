# Definitions YAML Directory

This directory contains YAML configuration files that define structural and labeling information for Photosystem II (PSII) and PsbS protein structures. These files are used in the analysis of binding sites and structural alignments in the PsbS Binding Site project.

## Files Overview

### Chain and Protein Definitions
- **`equivalent_chains.yaml`**: Groups equivalent chains by protein name, listing multiple chain IDs that represent the same protein type (e.g., CP24: ["8", "4"]).
- **`equivalent_chainids.yaml`**: Similar to `equivalent_chains.yaml`, but keyed by chain ID instead of protein name.
- **`chain_labels.yaml`**: Maps chain IDs (e.g., "8", "s", "c") to their corresponding protein names in the PSII complex (e.g., CP24, CP26, CP43, LHCB3, LHCBM, PSBX, etc.).
- **`psbs_labels.yaml`**: Maps chain ID "9" to the PsbS protein.
- **`psbs_labels_psii.yaml`**: Alternative or context-specific mapping of chain ID "9" to the PsbS protein, possibly for PSII-integrated analyses.

### Helix Definitions
- **`binding_proteins_helix_labels.yaml`**: Defines helix residue ranges for binding proteins in the PSII complex, organized by chain and using alphabetical helix labels (A, B, C, etc.), including start/end residues, residue names, length, and class.
- **`target_protein_helix_labels_merged.yaml`**: Merged helix definitions for the target protein (PsbS), organized by chain (A) and using helix types (TM1, H1, TM2, etc.), including start/end residues, residue names, length, and class.
- **`target_protein_helix_labels_merged_psii.yaml`**: Merged helix definitions for the target protein (PsbS) integrated with PSII data, organized by chain (9) and using helix types (TM1, H1, TM2, etc.).

## YAML File Formats

All YAML files follow standard YAML syntax. Below are the common data structures used:

- **Simple Key-Value Mappings**: Used in label files (e.g., `chain_labels.yaml`, `psbs_labels.yaml`).
  ```yaml
  "chain_id": "protein_name"
  ```

- **Helix Definitions**: Nested structure for helix ranges (e.g., all helix files).
  ```yaml
  chain_ID:
    helix_name:
      start: residue_number  # int
      end: residue_number    # int
      start_res_name: "RES"  # str, three-letter residue name
      end_res_name: "RES"    # str
      length: helix_length   # int
      class: helix_class     # int (typically 1 for alpha helices)
  ```

## Usage
These YAML files serve as configuration inputs for scripts that analyze protein structures, perform alignments, and identify binding sites. They provide standardized definitions for chains, residues, and structural elements across the PSII-PsbS complex.

## Notes
- All files use YAML format for easy parsing and editing.
- Helix definitions focus on class 1 (right-handed alpha) helices.</content>
<parameter name="filePath">/martini/rubiz/Github/PsbS_Binding_Site/analysis_dataset/definitions_yaml/README.md