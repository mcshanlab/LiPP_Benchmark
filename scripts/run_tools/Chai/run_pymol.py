import pymol
from pymol import cmd
import argparse
import os
import sys
import pandas as pd
import json
import numpy as np
from typing import Tuple, Dict

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *



def WriteSheet(dataframe: pd.DataFrame, BDID: str, rmsd: float, score_ptm: float, score_iptm: float) -> pd.DataFrame:
    """
    Add a new row to the dataframe with RMSD and score information.
    
    Args:
        dataframe: The existing dataframe to append to
        BDID: BioDolphin ID
        rmsd: RMSD value
        score_ptm: PTM score
        score_iptm: iPTM score
        
    Returns:
        Updated dataframe with new row appended
    """
    df_newrow = pd.DataFrame([[BDID,rmsd, score_ptm, score_iptm]], columns=['BioDolphinID', 'lipid_RMSD_pymol', 'ptm', 'iptm'])
    df_update = pd.concat([dataframe, df_newrow])

    #df_update = dataframe.append({'BioDolphinID': BDID, 'RMSD_pymol': rmsd}, ignore_index=True)

    return df_update


def GetBestScore(BDID: str) -> Tuple[str, float, float]:
    """
    Get the best scoring model based on iPTM score.
    
    Args:
        BDID: BioDolphin ID
        
    Returns:
        Tuple containing (model_name, iptm_score, ptm_score)
    """
    # get all scores
    score_path = f'./chai_output/{BDID}'
    iptm_dict: Dict[str, float] = {}
    ptm_dict: Dict[str, float] = {}

    for i in range(5):
        model_name = f'pred.model_idx_{i}.cif'
        model_scorefile = f'{score_path}/scores.model_idx_{i}.npz'
        scores = np.load(model_scorefile)
        #print(f'scores: {scores}')
        iptm, ptm = scores['iptm'], scores['ptm']
        iptm_dict[model_name] = iptm[0]
        ptm_dict[model_name] = ptm[0]

    print(iptm_dict)
    print(ptm_dict)

    topmodel = max(iptm_dict, key=iptm_dict.get) # get the top model's name (based on maximum iptm score)
    topmodel_iptm, topmodel_ptm = iptm_dict[topmodel], ptm_dict[topmodel]
    #print(f'max iptm model: {model_name_max}')

    return topmodel, topmodel_iptm, topmodel_ptm


def align_and_calculate_rmsd(BDID: str, chai_dir_path: str, pdb_path: str, output_dir_path: str) -> Tuple[float, float, float]:
    """
    Align model to real protein and calculate RMSD for lipid.
    
    Args:
        BDID: BioDolphin ID
        chai_dir_path: Path to chai output directory
        pdb_path: Path to PDB structures
        output_dir_path: Path to save aligned outputs
        
    Returns:
        Tuple containing (rmsd, topmodel_ptm, topmodel_iptm)
    """
    print(f'processing: {BDID}')

    # get the scores for each model
    topmodel, topmodel_iptm, topmodel_ptm = GetBestScore(BDID)
    print(f'topmodel: {topmodel}')
    print(f'topmodel_iptm: {topmodel_iptm}')
    print(f'topmodel_ptm: {topmodel_ptm}')
    
    model = f'{chai_dir_path}/{BDID}/{topmodel}'

    # pymol alignment:
    real_path = f'{pdb_path}/{BDID}/'
    real_lipid = real_path+'lipid.pdb'
    real_protein = real_path+'protein.pdb'

    # Load the structures
    cmd.load(real_lipid, "real_lipid")
    cmd.load(real_protein, "real_protein")
    cmd.load(model, "model_complex")

    # align the model to the real protein
    print(cmd.align("model_complex and chain A", "real_protein"))

    # create an object for the lipid of the model
    cmd.select("model_lipid", "model_complex and chain B")
    cmd.create("model_lipid", "model_lipid")

    # create an object for the protein of the model
    cmd.select("model_protein", "model_complex and chain A")
    cmd.create("model_protein", "model_protein")

    # calculate rmsd
    rmsd = cmd.rms_cur("model_lipid", "real_lipid", matchmaker=4)
    print(rmsd)

    # save the aligned model pdb
    os.makedirs(output_dir_path, exist_ok=True)
    save_pdb = f'{output_dir_path}/{BDID}_alignComplex.pdb'
    save_lipid = f'{output_dir_path}/{BDID}_alignLipid.pdb'
    save_protein = f'{output_dir_path}/{BDID}_alignProtein.pdb'

    cmd.save(filename=save_pdb, selection='model_complex')
    cmd.save(filename=save_lipid, selection='model_lipid')
    cmd.save(filename=save_protein, selection='model_protein')
    cmd.delete('all')
    
    return rmsd, topmodel_ptm, topmodel_iptm


def process_all_models(chai_dir_path: str = './chai_output', 
                       output_dir_path: str = './aligned_chai_output',
                       csv_path: str = './chai_result_pymol.csv') -> None:
    """
    Process all models in the chai output directory.
    
    Args:
        chai_dir_path: Path to chai output directory
        output_dir_path: Path to save aligned outputs
        csv_path: Path to save results CSV
    """
    config = load_config()
    pdb_path = config["pdb_path"]

    df = pd.DataFrame(columns=['BioDolphinID', 'lipid_RMSD_pymol', 'ptm', 'iptm'])
    BDIDs = os.listdir(chai_dir_path)
    
    for BDID in BDIDs:
        rmsd, topmodel_ptm, topmodel_iptm = align_and_calculate_rmsd(
            BDID, chai_dir_path, pdb_path, output_dir_path
        )
        
        # save the rmsd in the dataframe
        df = WriteSheet(df, BDID, rmsd, topmodel_ptm, topmodel_iptm)
        df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    process_all_models()

        
        