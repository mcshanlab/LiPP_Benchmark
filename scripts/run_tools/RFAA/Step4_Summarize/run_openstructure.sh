#!/bin/bash


mkdir -p ./open_structure_output

for mdl in ../Step3_Evaluation/rfaa_aligned_output/*_alignProtein.pdb; do
  # Extract ID (remove folder path and suffix)
  id=$(basename "$mdl" _alignProtein.pdb)

  # Build reference path
  ref="../../../data/pdb_structures_split/${id}/protein.pdb"

  echo "Comparing model: $mdl"
  echo "With reference:  $ref"

  ost compare-structures -m "$mdl" -r "$ref" \
      --lddt --tm-score --rigid-scores \
      -o "./open_structure_output/${id}.json"

done