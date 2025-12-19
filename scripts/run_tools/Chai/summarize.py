
import spyrmsd
from spyrmsd import io, rmsd
from spyrmsd.rmsd import rmsdwrapper
import argparse
import os, sys
import pandas as pd
import json
from typing import Tuple, Union


# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *



def GetRMSD_spy(ref_fn: str, out_fn: str) -> float:
    """
    Calculate RMSD using spyrmsd.
    
    Args:
        ref_fn: Path to ligand reference file (can be in pdb or sdf)
        out_fn: Path to ligand output file from docking (can be in pdb or sdf)
        
    Returns:
        RMSD value
    """
    ref = io.loadmol(ref_fn) # ligand reference file (can be in pdb or sdf)
    out = io.loadmol(out_fn) # ligand output file from docking (can be in pdb or sdf)
    rmsd_spy = rmsdwrapper(ref, out)[0].item() 
    return rmsd_spy

def AppendSheet(df: pd.DataFrame, BDID: str, rmsd_spy: Union[float, Exception], 
                lddt: Union[float, Exception], tm_score: Union[float, Exception], 
                rmsd_protein: Union[float, Exception]) -> pd.DataFrame:
    """
    Append a new row to the dataframe with calculated metrics.
    
    Args:
        df: The existing dataframe to append to
        BDID: BioDolphin ID
        rmsd_spy: RMSD value from spyrmsd (or Exception if failed)
        lddt: LDDT score (or Exception if failed)
        tm_score: TM score (or Exception if failed)
        rmsd_protein: Protein RMSD (or Exception if failed)
        
    Returns:
        Updated dataframe with new row appended
    """
    df_newrow = pd.DataFrame([[BDID,rmsd_spy, lddt, tm_score, rmsd_protein]], columns=['BioDolphinID', 'lipid_RMSD_spy', 'protein_lddt_openstc', 'protein_tm_openstc', 'protein_RMSD_openstc'])
    df_update = pd.concat([df, df_newrow])

    return df_update

def GetProteinScore(BDID: str) -> Tuple[float, float, float]:
    """
    Get the protein score from openstructure output.
    
    Args:
        BDID: BioDolphin ID
        
    Returns:
        Tuple containing (lddt, tm_score, rmsd)
            - lddt: LDDT score
            - tm_score: TM score
            - rmsd: C-alpha RMSD of the protein backbone
    """
    openstructure_json = f'./open_structure_output/{BDID}.json'

    with open(openstructure_json, "r") as f:
        data = json.load(f)
        lddt, tm_score, rmsd_protein = data.get("lddt"), data.get("tm_score"), data.get("rmsd")

    return lddt, tm_score, rmsd_protein


def process_single_bdid(BDID: str, pdb_path: str) -> Tuple[str, Union[float, Exception], Union[float, Exception], Union[float, Exception], Union[float, Exception]]:
    """
    Process a single BioDolphin ID to calculate metrics.
    
    Args:
        BDID: BioDolphin ID
        pdb_path: Path to PDB structures
        
    Returns:
        Tuple containing (BDID, rmsd_spy, lddt, tm_score, rmsd_protein)
    """
    print(f'processing BDID: {BDID}')
    
    # define paths
    real_lipid = f'{pdb_path}/{BDID}/lipid.pdb'
    chai_lipid = f'./aligned_chai_output/{BDID}_alignLipid.pdb'
    
    # calculate rmsd with spyrmsd
    try:
        rmsd_spy = GetRMSD_spy(real_lipid, chai_lipid)
        print(f'spy_rmsd: {rmsd_spy}')
    except Exception as e:
        print(e)
        print('cannot calculate spy_rmsd')
        rmsd_spy = e
    
    # Get protein score from openstructure
    try:
        lddt, tm_score, rmsd_protein = GetProteinScore(BDID)
        print(f'protein score from openstructure: lddt={lddt}, tm_score={tm_score}, rmsd_protein={rmsd_protein}')
    except Exception as e:
        print(e)
        print('cannot get protein score from openstructure')
        lddt, tm_score, rmsd_protein = e, e, e
    
    return BDID, rmsd_spy, lddt, tm_score, rmsd_protein


def summarize_results(chai_dir_path: str = './chai_output',
                     result_path: str = './chai_result_pymol.csv',
                     output_path: str = 'chai_result_all.csv') -> None:
    """
    Summarize all Chai results including RMSD and protein scores.
    
    Args:
        chai_dir_path: Path to chai output directory
        result_path: Path to the result file from pymol
        output_path: Path to save the final combined results
    """
    config = load_config()
    pdb_path = config["pdb_path"]
    
    print(f'reading outputs from this path: {chai_dir_path}')
    print(f'path to the output directory with the result from pymol: {result_path}')
    
    # get a list of all biodolphin IDs
    BDIDs = os.listdir(chai_dir_path)

    result_df = pd.read_csv(result_path)
    df = pd.DataFrame(columns=['BioDolphinID', 'RMSD_spy'])
    
    # loop over the output files
    for BDID in BDIDs:
        BDID, rmsd_spy, lddt, tm_score, rmsd_protein = process_single_bdid(BDID, pdb_path)
        df = AppendSheet(df, BDID, rmsd_spy, lddt, tm_score, rmsd_protein)

    # combine it with the existing result file
    final_df = pd.merge(result_df, df, on='BioDolphinID', how='left')
    final_df = final_df[['BioDolphinID', 'lipid_RMSD_pymol', 'lipid_RMSD_spy', 'protein_lddt_openstc', 'protein_tm_openstc', 'protein_RMSD_openstc','iptm', 'ptm']]
    final_df = final_df.sort_values('BioDolphinID')
    final_df.to_csv(output_path, index=False)
    print(f'saved results to {output_path}')


if __name__ == "__main__":
    summarize_results()
