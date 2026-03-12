Run ProteinCartography

## Installation
Install ProteinCartography from: https://github.com/Arcadia-Science/ProteinCartography
Follow the instructions there to install the associated environments

## Steps to run the code

### 1. Place the `cluster_settings` directory into the `demo` directory in your ProteinCartography installation

### 2. Place all PDB files under `cluster_settings/in_dir`

### 3. Run the program with 
- conda activate cartography_tidy
- snakemake --snakefile Snakefile --configfile demo/cluster_settings/config.yml --use-conda --cores 12

