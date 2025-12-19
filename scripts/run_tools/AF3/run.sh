#!/bin/bash

# Example of running AlphaFold3 on the preprocessed input files.
# Should be run by scripts for your cluster/configuration with gpu nodes.


CONFIG_FILE="../../../config.yaml"
AF3_path=$(grep '^AF3_path:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"')
AF3_image_path=$(grep '^AF3_image_path:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"')
AF3_database_path=$(grep '^AF3_database_path:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"')
AF3_param_path=$(grep '^AF3_param_path:' "$CONFIG_FILE" | awk '{print $2}' | tr -d '"')

echo 'Using AlphaFold3 path: ' $AF3_path
echo 'Using AlphaFold3 image file: ' $AF3_image_path
echo 'Using AlphaFold3 database path: ' $AF3_database_path
echo 'Using AlphaFold3 param path: ' $AF3_param_path


filenames=`ls ./af_input/BD*`

for eachfile in $filenames
do
   echo $eachfile
   prefix='./af_input/'
   suffix='.json'
   BDID=${eachfile//$prefix/}
   BDID=${BDID//$suffix/}
   echo $BDID
   mkdir -p ./af_output/${BDID}
   apptainer exec --bind ./af_input:/root/af_input --bind ./af_output/${BDID}:/root/af_output --bind ${AF3_param_path}:/root/models --bind ${AF3_database_path}:/root/public_databases --nv ${AF3_image_path} python ${AF3_path}/run_alphafold.py --json_path=/root/af_input/${BDID}.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output
done
