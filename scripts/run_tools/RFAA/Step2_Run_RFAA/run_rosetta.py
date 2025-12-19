import sys
import os
import subprocess
from typing import Dict, Optional
import pandas as pd

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *

def create_test_yaml(config_dir: str, job_name: str, fasta_path: str, smiles_string: str, test_rfaa_dir: str) -> str:
    """
    Creates a YAML configuration file for RoseTTAFold-All-Atom (RFAA) inference.
    
    This function generates a test.yaml file with the necessary configuration parameters
    for running RFAA on protein-small molecule complexes. The configuration includes
    protein input (FASTA file) and small molecule input (SMILES string).
    
    Args:
        config_dir (str): Directory path where the test.yaml file will be created.
        job_name (str): Name of the job for identification in RFAA output.
        fasta_path (str): Path to the protein FASTA file.
        smiles_string (str): SMILES representation of the small molecule.
        test_rfaa_dir (str): Directory where RFAA will write its output files.
    
    Returns:
        str: Full path to the created YAML configuration file.
    
    """

    yaml_content = f"""
defaults:
  - base
job_name: "{job_name}"
protein_inputs:
  A:
    fasta_file: {fasta_path}
sm_inputs:
  B:
    input: "{smiles_string}"
    input_type: "smiles"
hydra:
  run:
    dir: {test_rfaa_dir}
"""

    yaml_path = os.path.join(config_dir, "test.yaml")
    with open(yaml_path, "w") as yaml_file:
        yaml_file.write(yaml_content)
    print(f"Created: {yaml_path}")
    return yaml_path

def prepare_and_run_rfaa(main_directory: str, rfaa_software_dir: str, base_yaml_path: str, dolphin_smiles: Dict[str, str]) -> None:
    """
    Prepares input files, generates configuration files, and runs RoseTTAFold-All-Atom for protein-ligand complexes.
    
    This function processes multiple protein-ligand systems by:
    1. Scanning the main directory for subdirectories containing protein.fasta files
    2. Matching each directory with corresponding SMILES strings from the dolphin_smiles dictionary
    3. Creating RFAA configuration files for each valid protein-ligand pair
    4. Executing RFAA inference to predict protein-ligand complex structures
    
    Args:
        main_directory (str): Root directory containing subdirectories with protein.fasta files.
                             Each subdirectory represents one protein-ligand system.
        rfaa_software_dir (str): Path to the RoseTTAFold-All-Atom software installation directory.
        base_yaml_path (str): Path to the base YAML configuration file that contains default RFAA settings.
        dolphin_smiles (Dict[str, str]): Dictionary mapping folder names (BioDolphin IDs) to
                                        their corresponding SMILES strings.
    
    Returns:
        None: Function executes RFAA for each valid system and prints progress messages.
              Creates output directories and files for each successful run.
    """

    for folder in os.listdir(main_directory):
        folder_path = os.path.join(main_directory, folder)
        if os.path.isdir(folder_path):
            fasta_path = os.path.join(folder_path, "protein.fasta")

            # Check if required files exist
            if not os.path.isfile(fasta_path):
                print(f"Missing protein.fasta in {folder_path}, skipping...")
                continue

            # Get the SMILES string for this folder
            smiles_string = dolphin_smiles.get(folder)
            if not smiles_string:
                print(f"No SMILES string found for {folder}, skipping...")
                continue

            # Prepare inputs and config directory

            test_rfaa_dir = os.path.join(folder_path, "test_RFAA")
            os.makedirs(test_rfaa_dir, exist_ok=True)

            config_dir = os.path.join(test_rfaa_dir, "my_config")
            os.makedirs(config_dir, exist_ok=True)

            # Copy base.yaml file to the config directory
            base_yaml_dest = os.path.join(config_dir, "base.yaml")
            if not os.path.exists(base_yaml_dest):
                subprocess.run(["cp", base_yaml_path, base_yaml_dest])

            # Create the test.yaml file
            job_name = f"{folder}_output"
            yaml_path = create_test_yaml(config_dir, job_name, fasta_path, smiles_string, test_rfaa_dir)

            # Generate the RFAA run command
            run_command = (
                f"PYTHONPATH={rfaa_software_dir} python -m rf2aa.run_inference "
                f"--config-name test.yaml --config-path {config_dir}"
            )

            # Execute commands on the GPU node
            print(f"Running RFAA for {folder}...")
            print(f'Command to run is: {run_command}')
            subprocess.run(run_command, shell=True)
            print(f"Completed RFAA for {folder}.")


if __name__ == "__main__":

    
    config = load_config()
    current_dir = os.getcwd()
    parent_dir = os.path.abspath(os.path.join(current_dir, ".."))
    main_directory = os.path.join(parent_dir, "input")
    print(f"Main directory: {main_directory}")
    rfaa_software_directory = config["RFAA_path"]
    base_yaml_file_path ="./base.yaml"
    csv_directory = config["dolphinbusters_path"]

    df = pd.read_csv(csv_directory)
    # Create a dictionary mapping BioDolphinID to SMILES strings
    dolphin_smiles = df.set_index('BioDolphinID')['lipid_Canonical_smiles'].to_dict()

    # Run
    prepare_and_run_rfaa(main_directory, rfaa_software_directory, base_yaml_file_path, dolphin_smiles)
