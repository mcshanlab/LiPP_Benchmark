from pathlib import Path
import logging

import numpy as np
import torch
import sys
import os
import re
import shutil

from chai_lab.chai1 import run_inference
import argparse
import pandas as pd
from typing import List, Optional, Any

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *

def RunSingleFasta(fasta_path: Path, output_dir: Path) -> Optional[Any]:
    """Run inference for a single FASTA file and write results to `output_dir`.

    Returns whatever `run_inference` returns (often candidates) or None on failure.
    """
    try:
        candidates = run_inference(
            fasta_file=fasta_path,
            output_dir=output_dir,
            # 'default' setup
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device=torch.device("cuda:0"),
            use_esm_embeddings=True,
        )
        return candidates
    except Exception as exc:
        logging.exception(f"Failed to run inference for {fasta_path}: {exc}")
        return None


def parse_args() -> argparse.Namespace:
    """Parse command line arguments and return the namespace."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', "--start", type=int, default=0)
    parser.add_argument('-e', "--end", type=int, default=None)
    parser.add_argument('--items', nargs='+', help='A list of BDIDs to run on', default=None)
    return parser.parse_args()


def get_ids_from_range(dataset_path: str, start: int, end: int) -> List[str]:
    """Read the CSV at `dataset_path` and return a list of `BioDolphinID` strings
    using the slice [start:end] (end adjusted if out of range).
    """
    df = pd.read_csv(dataset_path)

    if start >= len(df):
        logging.error("start index is out of range")
        return []

    if end is None or end > len(df):
        df = df.iloc[start:]
    else:
        df = df.iloc[start:end]

    return df['BioDolphinID'].astype(str).tolist()


def process_ids(ids: List[str]) -> None:
    """Process a list of BioDolphin IDs by running inference on each FASTA."""
    for BD_id in ids:
        print(f'processing BioDolphin ID: {BD_id}')

        output_dir = Path(f'./chai_output/{BD_id}')
        if output_dir.exists():
            logging.warning(f"Removing old output directory: {output_dir}")
            shutil.rmtree(output_dir, ignore_errors=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        fasta_path = Path(f'./chai_input_fasta/{BD_id}.fasta')
        RunSingleFasta(fasta_path, output_dir)


def main(args: Optional[argparse.Namespace] = None) -> None:
    """Main entry point. Either use provided args or parse them from CLI."""
    if args is None:
        args = parse_args()

    IDs_to_run = args.items

    if IDs_to_run is None:
        print('preparing to run on a range of dataset using the index')
        config = load_config()
        dataset_path = config["LiPP_path"]
        startIDX, endIDX = args.start, args.end

        ids = get_ids_from_range(dataset_path, startIDX, endIDX)
        if not ids:
            print("No IDs to run (start index out of range or no rows). Ending process.")
            return

        process_ids(ids)
    else:
        print(f'preparing to run on: {IDs_to_run}')
        process_ids(IDs_to_run)
    



if __name__ == "__main__":
    main()
            
        
