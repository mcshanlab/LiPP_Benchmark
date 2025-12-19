
import pymol
from pymol import cmd
import argparse
import os, sys
from typing import Tuple, List, Optional
import pandas as pd
import json
import numpy as np
# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)
from src.load import *





def WriteSheet(dataframe: pd.DataFrame, BDID: str, rmsd: float, score: str) -> pd.DataFrame:
    """
    Appends a new row of evaluation results to the DataFrame.
    
    This function creates a new row containing BioDolphin ID, RMSD value, and confidence score,
    then concatenates it with the existing DataFrame to build the complete results table.
    
    Args:
        dataframe (pd.DataFrame): Existing DataFrame containing previous results.
        BDID (str): BioDolphin ID identifier for the current protein-ligand complex.
        rmsd (float): Root Mean Square Deviation value calculated by PyMOL.
        score (str): Confidence score extracted from DiffDock-L output filename.
    
    Returns:
        pd.DataFrame: Updated DataFrame with the new row appended.
    
    """
    df_newrow = pd.DataFrame([[BDID, rmsd, score]], columns=['BioDolphinID', 'RMSD_pymol', 'Confidence_Score'])
    df_update = pd.concat([dataframe, df_newrow])
    return df_update


def GetConfidence(BD_dir_path: str) -> Tuple[str, bool]:
    """
    Extracts confidence score from DiffDock-L output filenames.
    
    This function searches for files with the 'rank1_confidence' prefix in the specified directory
    and extracts the confidence score from the filename. It also validates that exactly one
    confidence file exists to avoid ambiguity.
    
    Args:
        BD_dir_path (str): Directory path containing DiffDock-L output files for one BioDolphin ID.
    
    Returns:
        Tuple[str, bool]: A tuple containing:
            - confidence_score (str): Extracted confidence value from filename
            - skip (bool): True if processing should be skipped (duplicate or missing files),
                          False if processing can continue

    """
    skip = False
    prefix = "rank1_confidence"
    matching_files = []
    for root, dirs, files in os.walk(BD_dir_path):
        for filename in files:
            if filename.startswith(prefix):
                matching_files.append(filename)
    if len(matching_files) != 1:
        skip = True
    confidence_score = matching_files[0].replace("rank1_confidence", "").replace(".sdf", "")
    return confidence_score, skip


def calculate_rmsd_for_system(BDID: str, diff_dir_path: str, dataset_path: str) -> Tuple[Optional[float], Optional[str], bool]:
    """
    Calculates RMSD between predicted and reference ligand structures for a single system.
    
    This function processes one BioDolphin ID by loading both the predicted ligand structure
    from DiffDock-L and the reference structure, then calculating RMSD using PyMOL.
    
    Args:
        BDID (str): BioDolphin ID identifier for the protein-ligand complex.
        diff_dir_path (str): Directory path containing DiffDock-L prediction outputs.
        dataset_path (str): Directory path containing reference structures.
    
    Returns:
        Tuple[Optional[float], Optional[str], bool]: A tuple containing:
            - rmsd (Optional[float]): Calculated RMSD value, None if calculation failed
            - confidence_score (Optional[str]): Extracted confidence score, None if extraction failed
            - skip (bool): True if system should be skipped, False if processing succeeded
    """
    model_lipid = f'{diff_dir_path}/{BDID}/rank1.sdf'
    
    if not os.path.exists(model_lipid):
        print(f'Skipping BDID: {BDID} - no output files found')
        return None, None, True
    
    confidence_score, skip = GetConfidence(f'{diff_dir_path}/{BDID}/')
    if skip:
        print(f'Skipping BDID: {BDID} - duplicate or missing confidence files')
        return None, None, True
    
    # Setup paths for reference structure
    real_path = f'{dataset_path}/{BDID}/'
    real_lipid = real_path + 'lipid.pdb'
    
    if not os.path.exists(real_lipid):
        print(f'Skipping BDID: {BDID} - reference lipid.pdb not found')
        return None, None, True
    
    try:
        # Load structures and calculate RMSD
        cmd.load(real_lipid, "real_lipid")
        cmd.load(model_lipid, "model_lipid")
        
        # Calculate RMSD using PyMOL
        rmsd = cmd.rms_cur("model_lipid", "real_lipid", matchmaker=4)
        print(f'BDID: {BDID}, RMSD: {rmsd:.3f}, Confidence: {confidence_score}')
        
        # Clean up PyMOL session
        cmd.delete('all')
        
        return rmsd, confidence_score, False
        
    except Exception as e:
        print(f'Error processing BDID: {BDID} - {str(e)}')
        cmd.delete('all')  # Clean up on error
        return None, None, True


def process_diffdock_results(diff_dir_path: str, csv_path: str) -> None:
    """
    Processes all DiffDock-L results and generates evaluation CSV.
    
    This function is the main processing pipeline that:
    1. Iterates through all BioDolphin IDs in the DiffDock output directory
    2. Calculates RMSD for each valid system
    3. Compiles results into a pandas DataFrame
    4. Saves the complete evaluation results to a CSV file
    
    Args:
        diff_dir_path (str): Directory containing DiffDock-L output subdirectories.
        csv_path (str): Output path for the evaluation results CSV file.
    
    Returns:
        None: Function creates CSV file and prints progress messages.
    """
    # Initialize results DataFrame
    df = pd.DataFrame(columns=['BioDolphinID', 'RMSD_pymol', 'Confidence_Score'])
    BDIDs = os.listdir(diff_dir_path)
    
    # Load configuration for dataset path
    config = load_config()
    dataset_path = config["pdb_path"]
    
    processed_count = 0
    skipped_count = 0
    
    print(f"Processing {len(BDIDs)} BioDolphin systems...")
    
    for BDID in BDIDs:
        rmsd, confidence_score, skip = calculate_rmsd_for_system(BDID, diff_dir_path, dataset_path)
        
        if not skip:
            df = WriteSheet(df, BDID, rmsd, confidence_score)
            processed_count += 1
        else:
            skipped_count += 1
    
    # Save results and print summary
    df.to_csv(csv_path, index=False)
    print(f"\nProcessing complete!")
    print(f"Processed: {processed_count} systems")
    print(f"Skipped: {skipped_count} systems")
    print(f"Results saved to: {csv_path}")


if __name__ == "__main__":
    """    
    Configuration:
        - diff_dir_path: Directory containing DiffDock-L outputs
        - csv_path: Output file for evaluation results
    """
    # Define input and output paths
    diff_dir_path = './diff_outputs'
    csv_path = './diff_result_pymol.csv'
    
    print("Starting DiffDock-L evaluation with PyMOL...")
    print(f"Input directory: {diff_dir_path}")
    print(f"Output CSV: {csv_path}")
    
    # Process all results
    process_diffdock_results(diff_dir_path, csv_path)


        