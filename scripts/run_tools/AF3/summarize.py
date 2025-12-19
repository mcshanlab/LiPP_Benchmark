import spyrmsd
from spyrmsd import io, rmsd
from spyrmsd.rmsd import rmsdwrapper
import argparse
import os
import sys
import pandas as pd
import json

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *
from typing import Tuple, Any, Optional

def GetRMSD_spy(ref_fn: str, out_fn: str) -> float:
    """Calculate RMSD between reference and output ligand using spyrmsd.
    
    Args:
        ref_fn: Path to reference ligand file (pdb or sdf)
        out_fn: Path to output ligand file from docking (pdb or sdf)
    
    Returns:
        RMSD value as a float
    """
    ref = io.loadmol(ref_fn) # ligand reference file (can be in pdb or sdf)
    out = io.loadmol(out_fn) # ligand output file from docking (can be in pdb or sdf)
    rmsd_spy = rmsdwrapper(ref, out)[0].item() 
    return rmsd_spy

def AppendSheet(df: pd.DataFrame, BDID: str, rmsd_spy: Any, lddt: Any, tm_score: Any, rmsd_protein: Any) -> pd.DataFrame:
    """Append a new row with RMSD and protein scores to the dataframe.
    
    Args:
        df: DataFrame to append to
        BDID: BioDolphin ID
        rmsd_spy: RMSD value from spyrmsd
        lddt: lDDT score from openstructure
        tm_score: TM-score from openstructure
        rmsd_protein: RMSD of protein backbone from openstructure
    
    Returns:
        Updated DataFrame with new row appended
    """
    df_newrow = pd.DataFrame([[BDID, rmsd_spy, lddt, tm_score, rmsd_protein]], columns=['BioDolphinID', 'lipid_RMSD_spy', 'protein_lddt_openstc', 'protein_tm_openstc', 'protein_RMSD_openstc'])
    df_update = pd.concat([df, df_newrow])

    return df_update


def GetProteinScore(BDID: str) -> Tuple[Optional[Any], Optional[Any], Optional[Any]]:
    """Get protein scores from openstructure output JSON file.
    
    Args:
        BDID: BioDolphin ID
    
    Returns:
        Tuple containing (lddt, tm_score, rmsd) from openstructure output.
        Returns None for each value if file is not found or error occurs.
        - lddt: Local Distance Difference Test score
        - tm_score: Template Modeling score
        - rmsd: C-alpha RMSD of the protein backbone
    """
    openstructure_json = f'./open_structure_output/{BDID}.json'

    try:
        with open(openstructure_json, "r") as f:
            data = json.load(f)

    except Exception as e:
        print(f"Skipped json file for {BDID} due to error: {e}")

    return data.get("lddt"), data.get("tm_score"), data.get("rmsd")


def main() -> None:
    """Main function to process AF3 output and calculate RMSD scores.
    
    Processes all AlphaFold3 output structures, calculates lipid RMSD using spyrmsd,
    retrieves protein scores from openstructure output, and generates a combined
    results CSV file.
    """
    # define arguments
    af_dir_path = './af_output'
    result_path = './af_result_pymol.csv'
    
    config = load_config()
    pdb_path = config["pdb_path"]
        
    
    # get a list of all biodolphin 
    BDIDs = os.listdir(af_dir_path)


    result_df = pd.read_csv(result_path)
    df = pd.DataFrame(columns=['BioDolphinID', 'lipid_RMSD_spy', 'protein_lddt_openstc', 'protein_tm_openstc', 'protein_RMSD_openstc'])
    
    # loop over the output files:
    for BDID in BDIDs:
        
        print(f'processing BDID: {BDID}')
        
        # define paths
        real_lipid = f'{pdb_path}/{BDID}/lipid.pdb'
        AF_lipid = f'./aligned_af_output/{BDID}_alignLipid.pdb'
        
        # calculate rmsd with spyrmsd
        try:
            rmsd_spy = GetRMSD_spy(real_lipid,AF_lipid)
            print(f'spy_rmsd: {rmsd_spy}')
        except Exception as e:
            print(e)
            print('cannot calculate spy_rmsd')
            rmsd_spy = e

        # Get protein score from openstructure
        try:
            lddt, tm_score, rmsd_protein = GetProteinScore(BDID)
        except Exception as e:
            print(e)
            print('cannot get protein score from openstructure')
            lddt, tm_score, rmsd_protein = e, e, e
        df = AppendSheet(df, BDID, rmsd_spy, lddt, tm_score, rmsd_protein)
            

    # combine it with the existing result file
    final_df = pd.merge(result_df, df, on='BioDolphinID', how='left')
    final_df = final_df[['BioDolphinID', 'lipid_RMSD_pymol', 'lipid_RMSD_spy', 'protein_RMSD_pymol', 'protein_lddt_openstc', 'protein_tm_openstc', 'protein_RMSD_openstc', 'iptm', 'ptm']]
    final_df = final_df.sort_values('BioDolphinID')
    final_df.to_csv('af_result_all.csv', index=False)
    print('saved results to af_result_all.csv')


if __name__ == "__main__":
    main()


        
        
