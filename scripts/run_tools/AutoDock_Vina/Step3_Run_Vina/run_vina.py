import os, sys
from typing import List
from Bio import PDB
import pandas as pd
import pymol
from pymol import cmd
from pymol.cgo import *
from random import randint

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)
from src.load import *

import spyrmsd
from spyrmsd import io, rmsd
import yaml


#############################################################################
#                                                                            
# drawBoundingBox.py -- Draws a box surrounding a selection 
#
#                                                                            
# AUTHOR: Jason Vertrees                                                   
# DATE  : 2/20/2009                                                          
# NOTES : See comments below.                                                
#                                                                            
#############################################################################
def drawBoundingBox(selection="(all)", padding=0.0, linewidth=2.0, r=1.0, g=1.0, b=1.0):     
        """                                                                  
        DESCRIPTION                                                          
                Given selection, draw the bounding box around it.          

        USAGE:
                drawBoundingBox [selection, [padding, [linewidth, [r, [g, b]]]]]

        PARAMETERS:
                selection,              the selection to enboxen.  :-)
                                        defaults to (all)

                padding,                defaults to 0

                linewidth,              width of box lines
                                        defaults to 2.0

                r,                      red color component, valid range is [0.0, 1.0]
                                        defaults to 1.0                               

                g,                      green color component, valid range is [0.0, 1.0]
                                        defaults to 1.0                                 

                b,                      blue color component, valid range is [0.0, 1.0]
                                        defaults to 1.0                                

        RETURNS
                string, the name of the CGO box

        NOTES
                * This function creates a randomly named CGO box that minimally spans the protein. The
                user can specify the width of the lines, the padding and also the color.                            
        """                                                                                                    

        ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

        print("Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))

        minX = minX - float(padding)
        minY = minY - float(padding)
        minZ = minZ - float(padding)
        maxX = maxX + float(padding)
        maxY = maxY + float(padding)
        maxZ = maxZ + float(padding)

        if padding != 0:
                 print("Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))

        boundingBox = [
                LINEWIDTH, float(linewidth),

                BEGIN, LINES,
                COLOR, float(r), float(g), float(b),

                VERTEX, minX, minY, minZ,       #1
                VERTEX, minX, minY, maxZ,       #2

                VERTEX, minX, maxY, minZ,       #3
                VERTEX, minX, maxY, maxZ,       #4

                VERTEX, maxX, minY, minZ,       #5
                VERTEX, maxX, minY, maxZ,       #6

                VERTEX, maxX, maxY, minZ,       #7
                VERTEX, maxX, maxY, maxZ,       #8


                VERTEX, minX, minY, minZ,       #1
                VERTEX, maxX, minY, minZ,       #5

                VERTEX, minX, maxY, minZ,       #3
                VERTEX, maxX, maxY, minZ,       #7

                VERTEX, minX, maxY, maxZ,       #4
                VERTEX, maxX, maxY, maxZ,       #8

                VERTEX, minX, minY, maxZ,       #2
                VERTEX, maxX, minY, maxZ,       #6


                VERTEX, minX, minY, minZ,       #1
                VERTEX, minX, maxY, minZ,       #3

                VERTEX, maxX, minY, minZ,       #5
                VERTEX, maxX, maxY, minZ,       #7

                VERTEX, minX, minY, maxZ,       #2
                VERTEX, minX, maxY, maxZ,       #4

                VERTEX, maxX, minY, maxZ,       #6
                VERTEX, maxX, maxY, maxZ,       #8

                END
        ]

        boxName = "box_" + str(randint(0,10000))
        while boxName in cmd.get_names():
                boxName = "box_" + str(randint(0,10000))

        cmd.load_cgo(boundingBox,boxName)
        return "Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)

cmd.extend("drawBoundingBox", drawBoundingBox)

def generate_grid(pdb_file: str, chain_id: str, residue_list: List, residue_numbers: List, failed_dockings=None):

    # Creating the PDB parser, please check https://biopython.org/docs/1.75/api/Bio.PDB
    parser = PDB.PDBParser(QUIET=True)

    # Parse the structure and creating empty list
    try:
        structure = parser.get_structure("protein", pdb_file)
    except FileNotFoundError:
        print(f"File not found: {pdb_file}")
        if failed_dockings is not None:
            failed_dockings.append(pdb_file)
        return None

    labeled_residue_coordinates = []
    not_labeled_residue_coordinates = []

    print("")
    print("")
    print("Current Protein and Chain:", pdb_file, chain_id)

    # Access the chain and extract coordinates for the specified residues
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    res_name = residue.get_resname()
                    res_id = residue.get_id()[1]
                    res_id = str(res_id)
                    res_id = "'" + res_id + "'"
                    if res_name in residue_list and res_id in residue_numbers:
                        print(f"Residue: {res_name} {res_id}")
                        for atom in residue:
                            print(f"Atom: {atom.get_name()}, Coordinates: {atom.coord}")
                            labeled_residue_coordinates.append((res_name, res_id, atom.get_name(), atom.coord))
                            not_labeled_residue_coordinates.append(atom.coord)
                    else:
                        print(f"No active residue found here: {res_id}")
            else:
                print(f"Chain {chain.id} not found in the structure.")

    # Split the list into 3 seperate coordinate lists
    x_values = [arr[0] for arr in not_labeled_residue_coordinates]
    y_values = [arr[1] for arr in not_labeled_residue_coordinates]
    z_values = [arr[2] for arr in not_labeled_residue_coordinates]

    # Calculate the centroid of the residues
    centroid_x = sum(x_values) / len(x_values)
    centroid_y = sum(y_values) / len(y_values)
    centroid_z = sum(z_values) / len(z_values)

    centroid_x = float(centroid_x)
    centroid_y = float(centroid_y)
    centroid_z = float(centroid_z)

    list_centers = [centroid_x, centroid_y, centroid_z]
    print("Autodock Center:", list_centers)
    return list_centers

def calculate_rmsd(ligand_file, output_file, output_directory, receptor_file):
    cmd.reinitialize()  # Reset PyMOL environment
    cmd.load(ligand_file, "ligand")
    cmd.load(output_file, "output")
    cmd.load(receptor_file, "receptor")
    print(cmd.get_object_list("all"))
    rmsd = cmd.rms_cur("ligand", "output")  # Calculate RMSD

    output_stem = os.path.splitext(os.path.basename(output_file))[0]
    aligned_path = os.path.join(output_directory, f"{output_stem}_aligned.pse")
    ref_path = os.path.join(output_directory, f"{output_stem}_reflig.pdb")
    vina_path = os.path.join(output_directory, f"{output_stem}_vinalig.pdb")

    cmd.save(aligned_path)
    pymol.cmd.save(ref_path, "ligand")
    pymol.cmd.save(vina_path, "output")

    cmd.delete("all")  # Reset again just in case
    return rmsd


def calculate_spy_rmsd(ligand_file, output_file, output_directory):
    output_stem = os.path.splitext(os.path.basename(output_file))[0]
    reflig_path = os.path.join(output_directory, f"{output_stem}_reflig.pdb")
    vinalig_path = os.path.join(output_directory, f"{output_stem}_vinalig.pdb")

    # SPYRMSD Calculation
    diff_ref = io.loadmol(reflig_path)  # ligand reference file (can be in pdb or sdf)
    diff_out = io.loadmol(vinalig_path)  # ligand output file from docking (can be in pdb or sdf)
    try:
        print(spyrmsd.rmsd.rmsdwrapper(diff_out, diff_ref))
        spyrmsd1 = spyrmsd.rmsd.rmsdwrapper(diff_out, diff_ref)
        spyrmsd1 = str(spyrmsd1)
    except Exception as e:
        print(f"Error: {e}")
        spyrmsd1 = "None: Error"
    return spyrmsd1

def extract_affinity(log_file):
    with open(log_file, 'r') as f:
        for line in f:
            if line.strip().startswith("1"):
                return float(line.split()[1])  # Extract the affinity value from the line
    return None  # Return None if no valid affinity is found


def build_vina_command(vina_binary, size_x, size_y, size_z, autodock_center, receptor_file, ligand_file, log_file, output_file):
    return (
        f"{vina_binary} --size_x {size_x} --size_y {size_y} --size_z {size_z} "
        f"--center_x {autodock_center[0]} --center_y {autodock_center[1]} --center_z {autodock_center[2]} "
        f"--receptor {receptor_file} --ligand {ligand_file} "
        f"--log {log_file} --out {output_file}"
    )


def append_results(output_excel, new_data):
    if os.path.exists(output_excel):
        with pd.ExcelWriter(output_excel, mode="a", engine="openpyxl", if_sheet_exists="overlay") as writer:
            new_data.to_excel(writer, index=False, header=False, startrow=writer.sheets["Sheet1"].max_row)
    else:
        new_data.to_excel(output_excel, index=False)


def process_rows(df, output_directory, vina_binary):
    failed_dockings = []
    output_excel = os.path.join(output_directory, "docking_results.xlsx")

    name = df['Name'].tolist()
    directory = df['Directory'].tolist()
    chain_ids = df['Protein Chain'].tolist()
    active_residue_names = df['Residue Names'].tolist()
    residue_ids = df['Residue Numbers'].tolist()
    x = df['X'].tolist()
    y = df['Y'].tolist()
    z = df['Z'].tolist()

    if len(directory) != len(chain_ids) != len(active_residue_names) != len(residue_ids):
        print("Corresponding lig and protein files have different number of lines.")

    for i in range(len(directory)):
        receptor_file = os.path.join(directory[i], "protein.pdbqt")
        ligand_file = os.path.join(directory[i], "lipid.pdbqt")
        log_file = os.path.join(output_directory, f"{os.path.basename(name[i])}_log.txt")
        output_file = os.path.join(output_directory, f"{os.path.basename(name[i])}_output.pdbqt")
        size_x = x[i]
        size_y = y[i]
        size_z = z[i]

        try:
            autodock_center = generate_grid(receptor_file, chain_ids[i], active_residue_names[i], residue_ids[i], failed_dockings)
            if autodock_center is None:
                print(f"Skipping {name[i]} due to missing autodock center.")
                continue
        except Exception as e:
            print(f"Error generating grid for {name[i]}: {e}")
            failed_dockings.append(name[i])
            continue

        try:
            command = build_vina_command(vina_binary, size_x, size_y, size_z, autodock_center, receptor_file, ligand_file, log_file, output_file)
        except Exception as e:
            print(f"Error building Vina command for {name[i]}: {e}")
            failed_dockings.append(name[i])
            continue

        print(f"Running Vina for: {directory[i]}")

        try:
            os.system(command)
        except Exception as e:
            print(f"Error: {e}")
            failed_dockings.append(name[i])
            continue

        try:
            rmsd = calculate_rmsd(ligand_file, output_file, output_directory, receptor_file)
            spy_rmsd = calculate_spy_rmsd(ligand_file, output_file, output_directory)
            affinity = extract_affinity(log_file)

            new_data = pd.DataFrame([[name[i], rmsd, spy_rmsd, affinity]], columns=["Name", "RMSD", "SPY_RMSD", "Affinity"])
            append_results(output_excel, new_data)
        except Exception as e:
            print(f"Error: {e}")
            continue

    print("\nFailed dockings:")
    for failed in failed_dockings:
        print(failed)
    print(f"\nDocking results saved to {output_excel}")


def parse_arguments(argv):
    if len(argv) < 4:
        print("Usage: python vina_script.py <input_excel> <task_id> <chunk_size>")
        sys.exit(1)

    excel_directory = argv[1]
    slurm_task_id = int(argv[2])
    chunk_size = int(argv[3])
    return excel_directory, slurm_task_id, chunk_size


def load_dataframe_slice(excel_directory, slurm_task_id=0, chunk_size=10):
    df_chunk = pd.read_excel(excel_directory, engine="openpyxl")
    #start_idx = slurm_task_id * chunk_size
    #end_idx = start_idx + chunk_size
    #return df_chunk.iloc[start_idx:end_idx]
    return df_chunk


def main():
    #excel_directory, slurm_task_id, chunk_size = parse_arguments(sys.argv)
    excel_directory = "../Step2_GridBox_Creation/output_residues.xlsx"
    output_directory = "./output"
    #vina_binary = "/storage/cedar/cedar0/cedarp-amcshan3-0/autodock_vina_1_1_2_linux_x86/bin/vina"
    config = load_config()
    vina_binary = config["AutoDockVina_path"]

    os.makedirs(output_directory, exist_ok=True)

    #df = load_dataframe_slice(excel_directory, slurm_task_id, chunk_size)
    df = load_dataframe_slice(excel_directory)
    process_rows(df, output_directory, vina_binary)


if __name__ == "__main__":
    main()
