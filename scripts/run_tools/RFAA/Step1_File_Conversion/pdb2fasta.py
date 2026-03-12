import os, sys
import subprocess
from typing import Optional
import yaml
import pandas as pd
from pathlib import Path


# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *

def extract_fasta_from_pdb(pdb_path: str) -> Optional[str]:
    """
    Extracts FASTA sequence from a PDB file using Open Babel.
    
    This function uses the Open Babel command-line tool (obabel) to convert
    a PDB file to FASTA format, extracting the protein sequence information.
    
    Args:
        pdb_path (str): Path to the input PDB file to process.
    
    Returns:
        Optional[str]: The FASTA sequence including header if successful,
                      None if an error occurs during processing.
    
    """
    try:
        # Run obabel command to convert PDB to FASTA
        result = subprocess.run(
            ["obabel", pdb_path, "-o", "fasta"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        # Check for errors in running the command
        if result.returncode != 0:
            print(f"Error processing {pdb_path}: {result.stderr.strip()}")
            return None
        # Return the FASTA sequence (including FASTA header)
        return result.stdout.strip()
    except Exception as e:
        print(f"Exception while processing {pdb_path}: {e}")
        return None

def process_folders(main_path: str) -> None:
    """
    Processes all folders to extract FASTA sequences from 'protein.pdb' files.
    
    This function iterates through all subdirectories in the specified path,
    looks for 'protein.pdb' files, extracts their FASTA sequences using Open Babel,
    and saves the sequences as 'protein.fasta' files in the same directories.
    
    Args:
        main_path (str): Root directory path containing subdirectories with protein.pdb files.
    
    Returns:
        None: Function creates FASTA files in place and prints progress messages.
    """
    # Iterate over each folder in the pdb folder
    for foldername in os.listdir(main_path):
        print(f"Processing folder: {foldername}")
        folder_path = os.path.join(main_path, foldername)
        if os.path.isdir(folder_path):  # Check if it's a folder
            protein_pdb_path = os.path.join(folder_path, "protein.pdb")
            
            if os.path.isfile(protein_pdb_path):  # Check if 'protein.pdb' exists
                # Extract FASTA sequence
                fasta_sequence = extract_fasta_from_pdb(protein_pdb_path)
                
                if fasta_sequence:
                    # Save the FASTA sequence to a .fasta file in the same folder
                    fasta_path = os.path.join(folder_path, "protein.fasta")
                    with open(fasta_path, "w") as fasta_file:
                        fasta_file.write(fasta_sequence)
                    print(f"FASTA file created: {fasta_path}")
            else:
                print(f"'protein.pdb' not found in folder: {foldername}")pdb' not found in folder: {foldername}")


def create_prep_directory(LiPP_path: str, pdb_path: str) -> None:
    """
    Creates a preparation directory with selected BioDolphin datasets for RFAA processing.
    
    This function reads a CSV file containing BioDolphin IDs, creates an '../input' directory,
    and copies the corresponding subdirectories from the source PDB path to prepare
    the input structure needed for RFAA workflow.
    
    Args:
        LiPP_path (str): Path to CSV file containing BioDolphin IDs in 'BioDolphinID' column.
        pdb_path (str): Source directory path containing subdirectories named by BioDolphin IDs.
    
    Returns:
        None: Function creates directory structure and prints status messages.
    """
    BDIDs_df = pd.read_csv(LiPP_path)
    BDIDs_list = BDIDs_df['BioDolphinID'].tolist()
    # make a input directory in current directory
    Path("../input").mkdir(parents=True, exist_ok=True)


    for BD_id in BDIDs_list:
        # copy it directories from pdb_path
        src_dir = os.path.join(pdb_path, BD_id)
        dest_dir = os.path.join("../input", BD_id)
        if os.path.exists(src_dir) and not os.path.exists(dest_dir):
            subprocess.run(f"cp -r {src_dir} {dest_dir}", shell=True)
        else:
            print(f"Source directory {src_dir} does not exist or destination directory {dest_dir} already exists. Skipping.")




if __name__ == "__main__":
    config = load_config()
    pdb_path = config["pdb_path"]
    LiPP_path = config["LiPP_path"]
    print(f'Using dataset path: {pdb_path}')
    print(f'Using LiPP_path: {LiPP_path}')
    create_prep_directory(LiPP_path, pdb_path)
    process_folders('../input')
    print("Processing complete.")
