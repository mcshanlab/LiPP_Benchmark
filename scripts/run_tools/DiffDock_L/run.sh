#!/usr/bin/env

# activate diffdock environment
micromamba activate diffdock
# Run diffdock inference
python -m inference --config default_inference_args.yaml --protein_ligand_csv all_input.csv  --out_dir diff_outputs