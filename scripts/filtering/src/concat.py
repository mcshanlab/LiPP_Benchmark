import os
from typing import List

import pandas as pd
from src.filter import *

def getFinal(repo_PATH: str, intermediate_PATH: str, final_fn: str = 'DolphinBuster_vr1.1.csv', res: float = 2) -> pd.DataFrame:
    """Concatenate distance-filtered subsets, apply resolution filter, and save.

    Args:
        repo_PATH: Root repository path where the final CSV is written.
        intermediate_PATH: Directory under repo_PATH containing df_filter_distance*.csv files.
        final_fn: Filename for the combined output CSV.
        res: Resolution cutoff to re-filter after concatenation.

    Returns:
        pd.DataFrame: Combined and resolution-filtered dataframe.
    """

    df_PATH = os.path.join(repo_PATH, intermediate_PATH)
    files: List[str] = [f for f in os.listdir(df_PATH) if f.startswith('df_filter_distance') and f.endswith('.csv')]
    print(f'concatenating files: {files}')
    dfs = [pd.read_csv(os.path.join(df_PATH, f)) for f in files]
    df_all = pd.concat(dfs, ignore_index=True)
    df_all = df_all.loc[:, ~df_all.columns.str.contains('^Unnamed')]
    df_all = Filter_Resolution(df_all, cutoff=res)
    df_all = df_all.sort_values(by='BioDolphinID')
    df_all.to_csv(os.path.join(repo_PATH, final_fn), index=False)

    return df_all

