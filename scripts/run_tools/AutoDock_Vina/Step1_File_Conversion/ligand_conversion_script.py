import os, sys
import subprocess
# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *


def main():
    """Main function to convert ligand PDB files to PDBQT format."""
    #config = load_config()
    #main_directory = config["pdb_path"]
    main_directory = './prep'
    print(f'Using dataset path: {main_directory}')
    
    # Loop through each folder inside the main directory
    counter = 0
    for root, dirs, files in os.walk(main_directory):
        counter += 1
        print(counter)

        for file in files:
            # Only process files named 'lipid.pdb'
            if file == "lipid.pdb":
                pdb_file = os.path.join(root, file)
                pdbqt_file = pdb_file.replace(".pdb", ".pdbqt")

                # Use obabel to convert pdb to pdbqt
                obabel_command = f"obabel {pdb_file} -xh -O {pdbqt_file}"
                print(f"Running: {obabel_command}")

                # Run the obabel command using subprocess
                subprocess.run(obabel_command, shell=True, check=True)

    print("Conversion completed!")


if __name__ == "__main__":
    main()
   
