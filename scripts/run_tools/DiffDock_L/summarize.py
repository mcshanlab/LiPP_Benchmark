

import spyrmsd
from spyrmsd import io, rmsd
from spyrmsd.rmsd import rmsdwrapper
import argparse
import os, sys
from typing import Union, Any
import pandas as pd
# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)
from src.load import *


def GetRMSD_spy(ref_fn: str, out_fn: str) -> float:
    """
    Calculates RMSD between reference and predicted ligand structures using SpyRMSD.
    
    SpyRMSD is specifically designed for small molecules and provides more accurate
    RMSD calculations compared to traditional algorithms by considering molecular
    symmetry and optimal atom matching.
    
    Args:
        ref_fn (str): Path to the reference ligand file (PDB or SDF format).
        out_fn (str): Path to the predicted ligand file from docking (PDB or SDF format).
    
    Returns:
        float: RMSD value in Angstroms calculated using SpyRMSD algorithm.
    """
    ref = io.loadmol(ref_fn)  # ligand reference file (can be in pdb or sdf)
    out = io.loadmol(out_fn) # ligand output file from docking (can be in pdb or sdf)
    rmsd_spy = rmsdwrapper(ref, out)[0].item() 
    return rmsd_spy

def AppendSheet(df: pd.DataFrame, BDID: str, rmsd_spy: Union[float, Exception]) -> pd.DataFrame:
    """
    Appends SpyRMSD results to the DataFrame.
    
    This function creates a new row containing BioDolphin ID and SpyRMSD value,
    then concatenates it with the existing DataFrame to build the SpyRMSD results table.
    
    Args:
        df (pd.DataFrame): Existing DataFrame containing previous SpyRMSD results.
        BDID (str): BioDolphin ID identifier for the current protein-ligand complex.
        rmsd_spy (Union[float, Exception]): SpyRMSD value or Exception object if calculation failed.
    
    Returns:
        pd.DataFrame: Updated DataFrame with the new row appended.
    
    """
    df_newrow = pd.DataFrame([[BDID, rmsd_spy]], columns=['BioDolphinID', 'RMSD_spy'])
    df_update = pd.concat([df, df_newrow])

    return df_update


def calculate_spyrmsd_for_system(BDID: str, dataset_path: str, diff_dir_path: str) -> Union[float, Exception]:
    """
    Calculates SpyRMSD for a single BioDolphin system with error handling.
    
    This function handles the complete SpyRMSD calculation workflow for one system,
    including file path construction, existence checking, and error handling.
    
    Args:
        BDID (str): BioDolphin ID identifier for the protein-ligand complex.
        dataset_path (str): Directory path containing reference structures.
        diff_dir_path (str): Directory path containing DiffDock-L prediction outputs.
    
    Returns:
        Union[float, Exception]: SpyRMSD value if successful, Exception object if failed.
    """
    print(f'Processing BDID: {BDID}')
    
    # Define file paths
    real_lipid = f'{dataset_path}/{BDID}/lipid.pdb'
    diff_lipid = f'{diff_dir_path}/{BDID}/rank1.sdf'
    
    # Calculate RMSD with SpyRMSD
    try:
        rmsd_spy = GetRMSD_spy(real_lipid, diff_lipid)
        print(f'SpyRMSD: {rmsd_spy:.3f}')
        return rmsd_spy
    except Exception as e:
        print(f'Error calculating SpyRMSD: {e}')
        print('Cannot calculate SpyRMSD for this system')
        return e


def process_spyrmsd_calculations(diff_dir_path: str, result_path: str, dataset_path: str) -> pd.DataFrame:
    """
    Processes SpyRMSD calculations for all systems and combines with PyMOL results.
    
    This function is the main processing pipeline that:
    1. Reads existing PyMOL results to get the list of systems to process
    2. Calculates SpyRMSD for each system
    3. Combines PyMOL and SpyRMSD results into a comprehensive evaluation
    4. Formats and sorts the final results
    
    Args:
        diff_dir_path (str): Directory containing DiffDock-L output subdirectories.
        result_path (str): Path to CSV file containing PyMOL results.
        dataset_path (str): Directory containing reference structures.
    
    Returns:
        pd.DataFrame: Combined DataFrame with PyMOL RMSD, SpyRMSD, and confidence scores.
    """
    print(f'Reading outputs from: {diff_dir_path}')
    print(f'Reading PyMOL results from: {result_path}')
    
    # Load existing PyMOL results
    result_df = pd.read_csv(result_path)
    BDIDs = result_df['BioDolphinID'].tolist()
    
    print(f'Processing {len(BDIDs)} systems for SpyRMSD calculation...')
    
    # Initialize SpyRMSD results DataFrame
    df = pd.DataFrame(columns=['BioDolphinID', 'RMSD_spy'])
    
    # Calculate SpyRMSD for each system
    for BDID in BDIDs:
        rmsd_spy = calculate_spyrmsd_for_system(BDID, dataset_path, diff_dir_path)
        df = AppendSheet(df, BDID, rmsd_spy)
    
    # Combine PyMOL and SpyRMSD results
    final_df = pd.merge(result_df, df, on='BioDolphinID', how='left')
    
    # Reorder and rename columns for clarity
    final_df = final_df[['BioDolphinID', 'RMSD_pymol', 'RMSD_spy', 'Confidence_Score']]
    final_df = final_df.rename(columns={
        'RMSD_pymol': 'lipid_RMSD_pymol', 
        'RMSD_spy': 'lipid_RMSD_spy'
    })
    
    # Sort by BioDolphin ID for consistent ordering
    final_df = final_df.sort_values('BioDolphinID')
    
    return final_df




if __name__ == "__main__":
    """
    Main execution block for DiffDock-L results summarization.
    
    Configuration:
        - diff_dir_path: Directory containing DiffDock-L outputs
        - result_path: CSV file with PyMOL results
        - dataset_path: Reference structures directory (from config)
    """
    # Define input and output paths
    diff_dir_path = './diff_outputs'
    result_path = './diff_result_pymol.csv'
    output_path = './diff_result_all.csv'
    
    # Load configuration
    config = load_config()
    dataset_path = config["pdb_path"]
    
    print("Starting DiffDock-L results summarization with SpyRMSD...")
    print(f"DiffDock outputs: {diff_dir_path}")
    print(f"PyMOL results: {result_path}")
    print(f"Reference dataset: {dataset_path}")
    print(f"Output file: {output_path}")
    
    # Process all SpyRMSD calculations and combine results
    final_df = process_spyrmsd_calculations(diff_dir_path, result_path, dataset_path)
    
    # Save final results
    final_df.to_csv(output_path, index=False)
    print(f'\nSaved comprehensive results to {output_path}')
    print(f'Total systems processed: {len(final_df)}')
    

