# The repository contains the methods to generate the simulations in the directories 4_pairs and 5_psii (Pending final documentation).

# The analysis_dataset directory contains scripts and instructions to process the simulation data and generate the analysis databases.

# Database Creation Pipeline

This pipeline processes protein binding event data from CSV files, enriches it with structural information, and creates combined databases for analysis.

## Quick Start

```bash
# Create a virtual environment
python -m venv venv
source venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt

# Run the pipeline inside analysis_dataset directory
bash write_databases.sh
```

## What it does

The `write_databases.sh` script performs three main steps:

### 1. Process individual chain data (`write_databases_pairs`)
- Reads CSV files containing binding events for chains A, B, and C
- Adds helix labels and chain identifiers from YAML definitions
- Applies SQL operations to add derived columns and classifications
- Outputs: `1_concatenated_databases/pairs/chain_{A,B,C}/database.csv`

### 2. Process supercomplex data (`write_databases_supercomplex`)
- Reads CSV files from PSII-LHCII supercomplex simulations
- Adds structural annotations and labels
- Applies the same SQL transformations
- Outputs: `1_concatenated_databases/supercomplex/database.csv`

### 3. Join all databases (`join_databases`)
- Combines the four databases (chains A, B, C, and supercomplex) into one
- Outputs: `2_joined_database/joined_database.csv`

## Directory Structure

```
analysis_dataset/
├── 0_datasets/               # Input CSV files
│   ├── pairs/
│   └── supercomplex/
├── 1_concatenated_databases/ # Processed individual databases
│   ├── pairs/
│   └── supercomplex/
├── 2_joined_database/        # Final combined database
├── definitions_yaml/         # YAML configuration files
├── sql_operations/           # SQL transformation files
└── scripts/                  # Python processing scripts
```

## Requirements

### Python Dependencies

Install required Python packages:

```bash
pip install -r requirements.txt
```

Required packages:
- `pandas==2.3.3` - Data manipulation and CSV processing
- `PyYAML==6.0.3` - YAML configuration file parsing

### Input Files

- YAML definition files in `definitions_yaml/`
- SQL operation files in `sql_operations/`
- Input CSV files in `0_datasets/`

## Output

Final database: `2_joined_database/joined_database.csv`

Contains enriched binding event data with:
- Lifetime statistics (mean, std, median, min, max, count, sum)
- Helix assignments for both binding partners
- Chain labels and identifiers
- Residue type classifications
- Derived columns from SQL operations
