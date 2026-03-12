import pandas as pd
import os, sys
import re
from openpyxl import Workbook
from pymol import cmd
from random import randint
from pymol.cgo import *

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *

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


# Function to extract the protein chain, ligand chain, and ligand name
def extract_chain_and_ligand(biodolphin_id):
    # Split the BioDolphinID by '-'
    parts = biodolphin_id.split('-')

    # Extract the protein chain, ligand chain, and ligand name
    protein_chain = parts[1]  # Protein chain letter
    ligand_chain = parts[2]   # Ligand chain letter
    ligand_name = parts[3]    # Ligand name
    ligand_name = ligand_name[:3]

    return protein_chain, ligand_chain, ligand_name

def generate_directory(biodolphin_id):
    return f"../Step1_File_Conversion/prep/{biodolphin_id}"

def process_pdb_file(row, ws, pdb_directory):
    """Process a single PDB file and extract residue information."""
    # Extract data from the CSV
    biodolphin_id = row['BioDolphinID']
    ligand_chain = row['Ligand_Chain']  
    protein_chain = row['Protein_Chain']
    ligand_resn = row['Ligand_Name']

    print(biodolphin_id, ligand_chain, protein_chain, ligand_resn)

    # Get the PDB file path
    pdb_file = os.path.join(pdb_directory, f"{biodolphin_id}.pdb")
    print(pdb_file)

    if os.path.exists(pdb_file):
        cmd.load(pdb_file)
        print(cmd.get_object_list())
        print("File Found")

        # PyMOL command to select residues within 15 Angstroms of ligand
        selectname = f"select_{biodolphin_id}"
        cmd.select(selectname, f"byres (chain {protein_chain} within 15 of chain {ligand_chain} and resn {ligand_resn})")
        cmd.select("ligand_selection", f"resn {ligand_resn}")
        
        # Iterate over selected residues and collect their information
        resnr = []
        resname = []

        def store_resnr(resi):
            resnr.append(resi)

        def store_resname(resn):
            resname.append(resn)

        cmd.iterate(f"{selectname} and name CA", "store_resnr(resi)", space=locals())
        cmd.iterate(f"{selectname} and name CA", "store_resname(resn)", space=locals())
        cmd.iterate(f"{selectname} and name CA", "print(resi,resn)")

        directory = generate_directory(biodolphin_id)
        
        # Ligand Size
        numbers = re.findall(r"[\d.]+", drawBoundingBox("ligand_selection"))
        dim1, dim2, dim3 = map(float, numbers)

        # Append data to Excel sheet
        ws.append([biodolphin_id, directory, protein_chain, str(resnr), str(resname), str(2 * dim1), str(2 * dim2), str(2 * dim3)])
        resnr.clear()
        resname.clear()
                
        # Clear the session for the next PDB file
        cmd.delete("all")
    else:
        print(f"PDB file for BioDolphin ID {biodolphin_id} not found.")


def main(csv_file, pdb_directory, output_excel):
    """Main function to process PDB files and generate pocket definitions."""
    # Load the CSV data
    df = pd.read_csv(csv_file)
    cmd.extend ("drawBoundingBox", drawBoundingBox)

    
    # Extract chain information
    df[['Protein_Chain', 'Ligand_Chain', 'Ligand_Name']] = df['BioDolphinID'].apply(lambda x: pd.Series(extract_chain_and_ligand(x)))
    df.to_csv('program_dataframe.csv', index=False)

    # Initialize Excel workbook
    wb = Workbook()
    ws = wb.active
    ws.title = "Residue Data"
    ws.append(["Name", "Directory", "Protein Chain", "Residue Numbers", "Residue Names", "X", "Y", "Z"])

    # Process each PDB file
    for index, row in df.iterrows():
        process_pdb_file(row, ws, pdb_directory)

    # Save the Excel file
    wb.save(output_excel)
    print(f"Residue data saved to {output_excel}")


if __name__ == "__main__":
    
    config = load_config()
    pdb_directory = config["pdb_complex_path"]
    csv_file = config["LiPP_path"]    
    output_excel = 'output_residues.xlsx'

    main(csv_file, pdb_directory, output_excel)
