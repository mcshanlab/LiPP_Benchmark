from typing import Tuple
import pymol
from pymol import cmd
import argparse
import os, sys
import pandas as pd
import json

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *



def WriteSheet(
    dataframe: pd.DataFrame,
    BDID: str,
    rmsd: float,
    protein_rmsd: float,
    score_ptm: float,
    score_iptm: float
) -> pd.DataFrame:
    """Add a new row to the results DataFrame with RMSD and confidence scores.
    
    Args:
        dataframe: Existing results DataFrame.
        BDID: BioDolphin identifier.
        rmsd: Ligand RMSD value.
        protein_rmsd: Protein RMSD value.
        score_ptm: pTM confidence score.
        score_iptm: ipTM confidence score.
        
    Returns:
        Updated DataFrame with new row appended.
    """
    df_newrow = pd.DataFrame(
        [[BDID, rmsd, protein_rmsd, score_ptm, score_iptm]],
        columns=['BioDolphinID', 'lipid_RMSD_pymol', 'protein_RMSD_pymol', 'ptm', 'iptm']
    )
    df_update = pd.concat([dataframe, df_newrow])
    return df_update


def ParseScore(score_file: str) -> Tuple[float, float]:
    """Parse confidence scores from AF3 summary JSON file.
    
    Args:
        score_file: Path to the AF3 summary confidences JSON file.
        
    Returns:
        Tuple of (pTM score, ipTM score).
    """
    with open(score_file) as json_file:
        json_dict = json.load(json_file)
        score_iptm = json_dict['iptm']
        score_ptm = json_dict['ptm']

    return score_ptm, score_iptm


def process_af_results(
    af_dir_path: str,
    output_dir_path: str,
    csv_path: str,
    pdb_path: str
) -> None:
    """Process AF3 results and calculate RMSD values using PyMOL.
    
    Iterates through AF3 output models, aligns them to reference structures,
    calculates protein and ligand RMSD values, and saves results to CSV and PDB files.
    
    Args:
        af_dir_path: Directory containing AF3 output results.
        output_dir_path: Directory where aligned structures will be saved.
        csv_path: Path to output CSV file for results.
        pdb_path: Base path to reference PDB structures.
        
    Returns:
        None. Results are saved to files.
    """
    df = pd.DataFrame(columns=['BioDolphinID', 'lipid_RMSD_pymol', 'protein_RMSD_pymol', 'ptm', 'iptm'])
    BDIDs = os.listdir(af_dir_path)

    for BDID in BDIDs:
        model = f'{af_dir_path}/{BDID}/{BDID.lower()}/{BDID.lower()}_model.cif'
        score_file = f'{af_dir_path}/{BDID}/{BDID.lower()}/{BDID.lower()}_summary_confidences.json'

        if os.path.exists(model):
            print(f'processing: {BDID}')
        else:
            print(f"The model for '{BDID}' does not exist.")
            continue

        score_ptm, score_iptm = ParseScore(score_file)

        real_path = f'{pdb_path}/{BDID}/'
        real_lipid = real_path + 'lipid.pdb'
        real_protein = real_path + 'protein.pdb'

        # Load the structures:
        cmd.load(real_lipid, "real_lipid")
        cmd.load(real_protein, "real_protein")
        cmd.load(model, "model_complex")

        # align the model to the real protein:
        cmd.align("model_complex and chain P", "real_protein")

        # create an object for the lipid of the model:
        cmd.select("model_lipid", "model_complex and chain L")
        cmd.create("model_lipid", "model_lipid")

        # create an object for the protein of the model:
        cmd.select("model_protein", "model_complex and chain P")
        cmd.create("model_protein", "model_protein")

        # calculate rmsd:
        # (1) lipid rmsd:
        lipid_rmsd = cmd.rms_cur("model_lipid", "real_lipid", matchmaker=4)
        # (2) protein rmsd:
        protein_rmsd = cmd.rms_cur("model_protein and name CA", "real_protein and name CA", matchmaker=4)

        # save the rmsd in the dataframe
        df = WriteSheet(df, BDID, lipid_rmsd, protein_rmsd, score_ptm, score_iptm)
        df.to_csv(csv_path, index=False)

        # save the aligned model pdb
        os.makedirs(output_dir_path, exist_ok=True)
        save_pdb = f'{output_dir_path}/{BDID}_alignComplex.pdb'
        save_lipid = f'{output_dir_path}/{BDID}_alignLipid.pdb'
        save_protein = f'{output_dir_path}/{BDID}_alignProtein.pdb'

        cmd.save(filename=save_pdb, selection='model_complex')
        cmd.save(filename=save_lipid, selection='model_lipid')
        cmd.save(filename=save_protein, selection='model_protein')
        cmd.delete('all')


if __name__ == "__main__":
    af_dir_path = './af_output'
    output_dir_path = './aligned_af_output'
    csv_path = './af_result_pymol.csv'
    config = load_config()
    pdb_path = config["pdb_path"]

    process_af_results(af_dir_path, output_dir_path, csv_path, pdb_path)
    
        
        