import os, sys
import subprocess
import yaml
import pandas as pd
from pathlib import Path


# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *

def split_alternate_conformations(file_path, split_alt_confs_script, conda_env_path):
    """Split alternate conformations from the PDB file."""
    try:
        out_path = file_path.replace(".pdb", "")
        split_command = f"{conda_env_path}/bin/python2 {split_alt_confs_script} -r {file_path} -o {out_path}"
        print(f"Running: {split_command}")
        subprocess.run(split_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error splitting alternate conformations for {file_path}: {e}")
        return False
    return True


def prepare_receptor(conf_path, pdbqt_file, prepare_receptor_script, conda_env_path):
    """Prepare receptor using MGLTools."""
    try:
        prepare_receptor_command = (
           f"{conda_env_path}/bin/python2 {prepare_receptor_script} -r {conf_path} -o {pdbqt_file} -A hydrogens"
        )

        print(f"Running: {prepare_receptor_command}")
        subprocess.run(prepare_receptor_command, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error preparing receptor for {conf_path}: {e}")
        return False


def process_protein_files(main_directory, split_alt_confs_script, prepare_receptor_script, conda_env_path):
    """Process protein PDB files and convert them to PDBQT format."""
    counter = 0
    
    # Loop through each folder inside the main directory
    for root, dirs, files in os.walk(main_directory):
        counter += 1
        print(counter)
        for file in files:
            file_path = os.path.join(root, file)

            if file == "protein.pdb":
                # Split alternate conformations if they exist
                print(f'file_path: {file_path}')
                if not split_alternate_conformations(file_path, split_alt_confs_script, conda_env_path):
                    continue

                # Check if conformation A exists after splitting
                conformation_A_found = False
                for conf_file in os.listdir(root):
                    print(f'conf_file: {conf_file}')
                    if conf_file.endswith("_A.pdb"):
                        # Define the path to conformation A
                        conf_A_path = os.path.join(root, conf_file)
                        print(f'conf_A_path: {conf_A_path}')
                        conformation_A_found = True
                        
                        # Define the output PDBQT file for conformation A
                        pdbqt_file = conf_A_path.replace("_A.pdb", ".pdbqt")
                        
                        # Prepare the receptor for conformation A using MGLTools
                        prepare_receptor(conf_A_path, pdbqt_file, prepare_receptor_script, conda_env_path)
                
                # If no conformation A was found, run on the original protein.pdb
                if not conformation_A_found:
                    pdbqt_file = file_path.replace(".pdb", ".pdbqt")
                    prepare_receptor(file_path, pdbqt_file, prepare_receptor_script, conda_env_path)

    print("Conversion completed!")


def create_prep_directory(LiPP_path: str, pdb_path: str) -> None:
    """Create a 'prep' directory for the converted inputs."""
    BDIDs_df = pd.read_csv(LiPP_path)
    BDIDs_list = BDIDs_df['BioDolphinID'].tolist()
    # make a prep directory in current directory
    Path("./prep").mkdir(parents=True, exist_ok=True)


    for BD_id in BDIDs_list:
        # copy it directories from pdb_path
        src_dir = os.path.join(pdb_path, BD_id)
        dest_dir = os.path.join("./prep", BD_id)
        if os.path.exists(src_dir) and not os.path.exists(dest_dir):
            subprocess.run(f"cp -r {src_dir} {dest_dir}", shell=True)
        else:
            print(f"Source directory {src_dir} does not exist or destination directory {dest_dir} already exists. Skipping.")




if __name__ == "__main__":

    config = load_config()
    conda_env_path = config["conda_path"]
    pdb_path = config["pdb_path"]
    LiPP_path = config["LiPP_path"]
    create_prep_directory(LiPP_path, pdb_path)
    
    main_directory = './prep'
    # Set up conda environment path:
    conda_env_path = f"{conda_env_path}/openbabel_env"
    
    # Path to the prepare_pdb_split_alt_confs.py script
    split_alt_confs_script = f"{conda_env_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py"
    prepare_receptor_script = f"{conda_env_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"


    process_protein_files(main_directory, split_alt_confs_script, prepare_receptor_script, conda_env_path)
