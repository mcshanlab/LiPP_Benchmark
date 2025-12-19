import cmd
import os, sys
from typing import List, Tuple, Optional, Union, Any
import openpyxl
import pymol2
import torch
import pandas as pd

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *

def extract_pae_inter_and_update_excel(input_dir: str, output_file: str) -> None:
    """
    Extracts PAE (Predicted Aligned Error) inter values from PyTorch files and updates an Excel file.
    
    This function walks through a directory structure to find .pt or .pth files, extracts the
    'pae_inter' values from them, calculates the mean value, and appends these values to an
    existing Excel file in the 'pae_inter' column.
    
    Args:
        input_dir (str): Directory path to search for PyTorch files (.pt or .pth).
        output_file (str): Path to the Excel file to be updated with PAE inter values.
    
    Returns:
        None: Function modifies the Excel file in place and prints status messages.
    
    """
    # Collect all PyTorch files in the directory
    pytorch_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".pt") or file.endswith(".pth"):
                pytorch_files.append(os.path.join(root, file))

    if not pytorch_files:
        print("No PyTorch files found in the input directory.")
        return

    # Load the existing Excel file
    try:
        excel_data = pd.read_excel(output_file)
    except FileNotFoundError:
        print(f"Error: The Excel file '{output_file}' does not exist.")
        return

    # Ensure the Excel file has a column to append pae_inter values
    if "pae_inter" not in excel_data.columns:
        excel_data["pae_inter"] = None

    # Extract pae_inter values from each PyTorch file and append to the Excel file
    for pytorch_file in pytorch_files:
        try:
            data = torch.load(pytorch_file, map_location="cpu")
            if "pae_inter" in data:
                pae_inter_value = data["pae_inter"]
                if isinstance(pae_inter_value, (list, torch.Tensor)):
                    mean_pae_inter = torch.tensor(pae_inter_value).float().mean().item()
                else:
                    mean_pae_inter = float(pae_inter_value)  # Directly use the float value
                
                print(f"Extracted mean pae_inter: {mean_pae_inter} from {pytorch_file}")

                # Add the mean pae_inter value to the next available row in the "pae_inter" column
                row_index = excel_data[excel_data["pae_inter"].isna()].index[0]
                excel_data.at[row_index, "pae_inter"] = mean_pae_inter
            else:
                print(f"Warning: 'pae_inter' key not found in {pytorch_file}.")
        except Exception as e:
            print(f"Error loading {pytorch_file}: {e}")

    # Save the updated Excel file
    excel_data.to_excel(output_file, index=False)
    print(f"Updated Excel file saved to {output_file}.")


def get_ligand_info(pdb_file: str) -> List[Tuple[str, str]]:
    """
    Extracts ligand information from a PDB file using PyMOL.
    
    This function loads a PDB file and identifies all non-polymer, non-solvent residues
    (ligands), returning their residue names and chain identifiers.
    
    Args:
        pdb_file (str): Path to the PDB file to analyze.
    
    Returns:
        List[Tuple[str, str]]: List of tuples containing (residue_name, chain_id)
                              for each unique ligand found in the structure.
                              Duplicates are removed.

    """
    ligand_info: List[Tuple[str, str]] = []
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb_file, "structure")

        # Iterate over all residues in the structure
        for model in pymol.cmd.get_object_list("structure"):
            pymol.cmd.iterate(f"{model} and not polymer and not solvent", 
                              "ligand_info.append((resn, chain))", 
                              space={'ligand_info': ligand_info})
    
    # Remove duplicates
    ligand_info = list(set(ligand_info))

    return ligand_info

def find_pdb_files(directory: str) -> List[str]:
    """
    Recursively finds all PDB files in the given directory and its subdirectories.
    
    Args:
        directory (str): Root directory path to search for PDB files.
    
    Returns:
        List[str]: List of absolute paths to all PDB files found.
                  Files are identified by having a '.pdb' extension.
    
    """
    pdb_files: List[str] = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(root, file))
    return pdb_files

def get_related_files(main_pdb: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Finds the corresponding lipid.pdb and protein.pdb files for a given main PDB file.
    
    This function constructs the expected paths for related files based on the main PDB
    filename and checks if they exist in the '../input' directory structure.
    
    Args:
        main_pdb (str): Path to the main PDB file (expected to end with '_output.pdb').
    
    Returns:
        Tuple[Optional[str], Optional[str]]: A tuple containing:
            - Path to lipid.pdb file if it exists, None otherwise
            - Path to protein.pdb file if it exists, None otherwise

    """
    base_dir = os.path.dirname(main_pdb)
    base_name = os.path.basename(main_pdb).replace("_output.pdb", "")
    #target_dir = os.path.join(base_dir, base_name)
    target_dir = os.path.join('../input', base_name)

    if os.path.exists(target_dir):
        lipid_pdb = os.path.join(target_dir, "lipid.pdb")
        protein_pdb = os.path.join(target_dir, "protein.pdb")
        if os.path.exists(lipid_pdb) and os.path.exists(protein_pdb):
            return lipid_pdb, protein_pdb
    return None, None

def process_pdb_files(input_dir: str, output_dir: str, output_excel: str, spyrmsd_directory: str) -> None:
    """
    Processes PDB files to calculate RMSD values using PyMOL and SpyRMSD, saving results to Excel.
    
    This function performs comprehensive RMSD analysis on protein-ligand structures by:
    1. Finding all '_output.pdb' files in the input directory
    2. Locating corresponding reference files (lipid.pdb and protein.pdb)
    3. Aligning structures using PyMOL
    4. Calculating RMSD using both PyMOL and SpyRMSD methods
    5. Saving aligned structures and results to files
    
    Args:
        input_dir (str): Directory containing the output PDB files to process.
        output_dir (str): Directory where aligned PDB files will be saved.
        output_excel (str): Path to Excel file where RMSD results will be stored.
        spyrmsd_directory (str): Directory where temporary files for SpyRMSD calculation will be saved.
    
    Returns:
        None: Function creates output files and updates the Excel file with RMSD values.
              Each row contains [Biodolphin ID, RMSD_Pymol, RMSD_SPY].

    """
    pdb_files = find_pdb_files(input_dir)
    print(f'pdb_files: {pdb_files}')

    # Initialize Excel workbook
    if os.path.exists(output_excel):
        workbook = openpyxl.load_workbook(output_excel)
    else:
        workbook = openpyxl.Workbook()
        workbook.active.append(["Biodolphin ID", "RMSD_Pymol", "RMSD_SPY"])
    sheet = workbook.active

    # PyMOL processing
    with pymol2.PyMOL() as pymol:
        pymol.cmd.reinitialize()

        for pdb_file in pdb_files:
            if "_output.pdb" not in pdb_file:
                continue

            lipid_pdb, protein_pdb = get_related_files(pdb_file)
            if not lipid_pdb or not protein_pdb:
                print(f"Missing related files for {pdb_file}")
                continue

            import spyrmsd
            from spyrmsd import io, rmsd

            biodolphin_id = os.path.basename(pdb_file).replace("_output.pdb", "")
            
            # Load PDB files
            pymol.cmd.load(pdb_file, "rosetta_output")
            pymol.cmd.load(lipid_pdb, "lipid")
            pymol.cmd.load(protein_pdb, "protein")

            # Create reference object
            pymol.cmd.select("reference_selection", "protein, lipid")
            pymol.cmd.create("reference", "reference_selection")
            pymol.cmd.delete("lipid")
            pymol.cmd.delete("protein")

            # Align reference and Rosetta output
            pymol.cmd.align("rosetta_output and polymer.protein", "reference and polymer.protein")
            
            # Calculate RMSD for ligands only
            print(lipid_pdb)
            parts = biodolphin_id.split('-')
            protein_chain = parts[1]  # Protein chain letter
            ligand_chain = parts[2]   # Ligand chain letter
            ligand_name = parts[3]    # Ligand name
            ligand_name = ligand_name[:3]

            ligands = get_ligand_info(pdb_file)
            name = ligands[0][0]
            chain = ligands[0][1]

            # Pymol RMSD
            pymol.cmd.select("reference_ligand", f"resn {ligand_name} and chain {ligand_chain}")
            pymol.cmd.select("rosetta_ligand", f"resn {name} and chain {chain}")
            rmsd = pymol.cmd.rms_cur("rosetta_ligand", "reference_ligand", matchmaker=4)

            # Save aligned
            pymol.cmd.save(os.path.join(output_dir, f"{biodolphin_id}_aligned.pdb"))
            pymol.cmd.save(filename=os.path.join(output_dir, f"{biodolphin_id}_aligned.pdb"), selection='rosetta_output')
            pymol.cmd.save(filename=os.path.join(output_dir, f"{biodolphin_id}_alignProtein.pdb"), selection='rosetta_output and polymer.protein')

            # SPYRMSD Prep
            pymol.cmd.create("reflig", "reference_ligand")
            pymol.cmd.create("rosettalig", "rosetta_ligand")
            ref_path = os.path.join(spyrmsd_directory, f"{biodolphin_id}_reflig.pdb")
            rosetta_path = os.path.join(spyrmsd_directory, f"{biodolphin_id}_rosettalig.pdb")
            pymol.cmd.save(ref_path, "reflig")
            pymol.cmd.save(rosetta_path, "rosettalig")

            # Re-add protein to the environment
            pymol.cmd.delete("rosetta_output")
            pymol.cmd.delete("reference")
            pymol.cmd.load(protein_pdb, "protein")

            # Save the session
            #pymol.cmd.save(os.path.join(output_dir, f"{biodolphin_id}_aligned.pse"))
            pymol.cmd.delete("all")

            # SPYRMSD Calculation
            diff_ref = io.loadmol(ref_path) # ligand reference file (can be in pdb or sdf)
            diff_out = io.loadmol(rosetta_path) # ligand output file from docking (can be in pdb or sdf)
            try:
                print(spyrmsd.rmsd.rmsdwrapper(diff_out, diff_ref))
                spyrmsd1 = spyrmsd.rmsd.rmsdwrapper(diff_out, diff_ref)
                spyrmsd1 = str(spyrmsd1)
            except Exception as e:
                print(f"Error: {e}")
                spyrmsd1 = "None: Error"

            # Append RMSD to Excel
            sheet.append([biodolphin_id, rmsd, spyrmsd1])
            
        # Save the workbook
        workbook.save(output_excel)


if __name__ == "__main__":
    
    config = load_config()
    input_directory = '../Step2_Run_RFAA'
    output_directory = './rfaa_aligned_output'
    output_excel_file = './rosetta.xlsx'
    spyrmsd_directory = './lipid_outputs'

    os.makedirs(output_directory, exist_ok=True)
    os.makedirs(spyrmsd_directory, exist_ok=True)

    process_pdb_files(input_directory, output_directory, output_excel_file, spyrmsd_directory)
    extract_pae_inter_and_update_excel(input_directory, output_excel_file)

