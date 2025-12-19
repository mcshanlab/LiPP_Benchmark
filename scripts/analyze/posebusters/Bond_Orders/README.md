# recover.py Usage Guide

## Overview
`recover.py` is a helper script for recovering and assigning bond orders to PDB ligand files, adapted from Pat Walters' blog post "Assigning Bond Orders to PDB Ligands." It processes a set of lipid PDB files, assigns correct bond orders using Ligand Expo SMILES templates, and outputs SDF files for further analysis.


## Usage
Run the script from the command line:

```bash
python recover.py --output_pdb_dir <path_to_output_dir> --output_sdf_dir <path_to_sdf_dir>
```

- `--output_pdb_dir`: Directory containing the output lipid files for processing.
- `--output_sdf_dir`: Directory where the output SDF files will be saved.
