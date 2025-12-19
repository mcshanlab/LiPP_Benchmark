import pandas as pd
import numpy as np
import argparse
import random
import itertools
from typing import List, Tuple, Union, Optional
from plot import *
import os, sys
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# Add the scripts directory to Python path to allow importing src module
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, '..', '..')
sys.path.insert(0, scripts_dir)

from src.load import *



def AddPrefix(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """
    Add a prefix to column names except 'BioDolphinID'.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with RMSD data from a prediction method.
    prefix : str
        Prefix to add to column names (e.g., 'AF', 'CHAI', 'VINA', 'RS', 'DD').
    
    Returns
    -------
    pd.DataFrame
        DataFrame with prefixed column names. 'BioDolphinID' column remains unchanged.
    """
    df.columns = [f'{prefix}_{col}' if col != 'BioDolphinID' else col for col in df.columns]
    return df



def TagClasses(df_rmsds: pd.DataFrame, df_dolbuster: pd.DataFrame, add_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Merge RMSD results with DolphinBusters metadata tags (protein families, lipid classes, etc.).
    
    Parameters
    ----------
    df_rmsds : pd.DataFrame
        Input dataframe containing RMSD results with 'BioDolphinID' column.
    df_dolbuster : pd.DataFrame
        DolphinBusters dataset containing metadata annotations.
    add_cols : Optional[List[str]]
        Columns to merge from DolphinBusters. Defaults to biological metadata columns.
    
    Returns
    -------
    pd.DataFrame
        Merged dataframe with RMSD results plus metadata columns joined on 'BioDolphinID'.
    """
    if add_cols is None:
        add_cols = ['BioDolphinID', 'protein_Pfam', 'protein_Pfam_ID', 'lipid_Lipidmaps_categories', 'lipid_Lipidmaps_terms', 'lipid_Molecular_weight', 'PDB_release_date']
    df_dolbuster_labels = df_dolbuster[add_cols]
    tagged_df = pd.merge(df_rmsds, df_dolbuster_labels, on='BioDolphinID', how='left')
    return tagged_df




def AvgRMSD(df: pd.DataFrame, rmsd_method: str = 'spy', prefix: str = 'AF') -> float:
    """
    Calculate average RMSD for a prediction method.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with RMSD values.
    rmsd_method : str, default='spy'
        RMSD calculation method ('spy' or 'pymol').
    prefix : str, default='AF'
        Method prefix in column names (e.g., 'AF', 'CHAI', 'VINA', 'RS', 'DD').
    
    Returns
    -------
    float
        Average RMSD value across all predictions in Angstroms.
    """
    if rmsd_method == 'pymol':
        col_name = f'{prefix}_RMSD_pymol'
    if rmsd_method == 'spy':
        col_name = f'{prefix}_RMSD_spy'
    avg_rmsd = df[col_name].mean()
    return avg_rmsd


def SucessRates(df: pd.DataFrame, cutoffs: Optional[List[float]] = None, rmsd_method: str = 'spy', prefix: str = 'AF') -> List[float]:
    """
    Calculate success rates at multiple RMSD cutoff thresholds.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with RMSD values.
    cutoffs : Optional[List[float]], default=[2, 2.5, 3]
        RMSD cutoff values in Angstroms. Defaults to [2.0, 2.5, 3.0] Angstroms.
    rmsd_method : str, default='spy'
        RMSD calculation method ('spy' or 'pymol').
    prefix : str, default='AF'
        Method prefix in column names.
    
    Returns
    -------
    List[float]
        Success rates (0-1) for each cutoff value. Returns empty list if no data.
    """
    if cutoffs is None:
        cutoffs = [2, 2.5, 3]
    
    if rmsd_method == 'pymol':
        col_name = f'{prefix}_lipid_RMSD_pymol'
        rmsd_list = df[col_name].tolist()
    if rmsd_method == 'spy':
        col_name = f'{prefix}_lipid_RMSD_spy'
        rmsd_list = df[col_name].tolist()

    sucessrates = []
    for cutoff in cutoffs:
        NumSucessed = sum(rmsd < cutoff for rmsd in rmsd_list)
        if len(rmsd_list) == 0:
            sucessrate = 0
        else:
            sucessrate = NumSucessed/len(rmsd_list)
        sucessrates.append(sucessrate)
    return sucessrates
    

def Culmulate_SucessRate(df: pd.DataFrame, cutoff: Optional[np.ndarray] = None, rmsd_method: str = 'spy', prefix: str = 'AF') -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate cumulative success rate across a range of RMSD cutoff values.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with RMSD values.
    cutoff : Optional[np.ndarray], default=np.arange(0, 11, 0.1)
        Array of RMSD cutoff values in Angstroms. Defaults to 0.0 to 10.9 in 0.1 Angstrom steps.
    rmsd_method : str, default='spy'
        RMSD calculation method ('spy' or 'pymol').
    prefix : str, default='AF'
        Method prefix in column names.
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (cumulative_freq, cutoff) where cumulative_freq is success rate per cutoff.
    """
    if cutoff is None:
        cutoff = np.arange(0, 11, 0.1)
    
    if rmsd_method == 'pymol':
        col_name = f'{prefix}_lipid_RMSD_pymol'
        rmsd_list = df[col_name].tolist()
    if rmsd_method == 'spy':
        col_name = f'{prefix}_lipid_RMSD_spy'
        rmsd_list = df[col_name].tolist()
    
    culmulate_freq = []

    for i, c in enumerate(cutoff):
        NumSucessed = sum(rmsd < c for rmsd in rmsd_list)
        sucessrate = NumSucessed/len(rmsd_list)
        culmulate_freq.append(sucessrate)
    
    return np.array(culmulate_freq), cutoff


def GetConfidenceInterval(df_subset: pd.DataFrame, successrate: float) -> List[float]:
    """
    Calculate 95% confidence interval for a success rate using Wilson score method.
    
    Parameters
    ----------
    df_subset : pd.DataFrame
        Subset of data used to compute the success rate.
    successrate : float
        Success rate value between 0 and 1.
    
    Returns
    -------
    List[float]
        [lower_bound, upper_bound] of 95% confidence interval using beta method.
    """
    successcount = int(successrate * len(df_subset))
    lower, upper = proportion_confint(count=successcount, nobs=len(df_subset), alpha=0.05, method='beta')
    return [lower, upper]



def GetSubset_LipidClass(df: pd.DataFrame, subset: str = 'Glycerophospholipid') -> pd.DataFrame:
    """
    Filter dataframe to include only a specific lipid classification.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with 'lipid_Lipidmaps_categories' column.
    subset : str, default='Glycerophospholipid'
        Lipid class: 'Glycerophospholipid' (GP), 'Fatty_Acyl' (FA), 'Glycerolipid' (GL),
        'Polyketide' (PK), 'Prenol_Lipid' (PR), 'Saccharolipid' (SL), 'Sphingolipid' (SP),
        or 'Sterol_Lipid' (ST).
    
    Returns
    -------
    pd.DataFrame
        Filtered dataframe containing only entries matching the specified lipid class.
    """
    if subset == 'Glycerophospholipid':
        keyword = 'GP'
    elif subset == 'Fatty_Acyl':
        keyword = 'FA'
    elif subset == 'Glycerolipid':
        keyword = 'GL'
    elif subset == 'Polyketide':
        keyword = 'PK'
    elif subset == 'Prenol_Lipid':
        keyword = 'PR'
    elif subset == 'Saccharolipid':
        keyword = 'SL'
    elif subset == 'Sphingolipid':
        keyword = 'SP'
    elif subset == 'Sterol_Lipid':
        keyword = 'ST'
    
    subset_df = df[df['lipid_Lipidmaps_categories'].str.contains(keyword, na=False)]

    return subset_df


def GetSubset_LipidMW(df: pd.DataFrame, range_start: float = 0, range_end: float = 100) -> pd.DataFrame:
    """
    Filter dataframe to include only lipids within a molecular weight range.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with 'lipid_Molecular_weight' column in Daltons.
    range_start : float, default=0
        Minimum molecular weight (Da, inclusive).
    range_end : float, default=100
        Maximum molecular weight (Da, exclusive).
    
    Returns
    -------
    pd.DataFrame
        Filtered dataframe with lipids in specified MW range [range_start, range_end).
    """
    subset_df = df[(df['lipid_Molecular_weight'] >= range_start) & (df['lipid_Molecular_weight'] < range_end)]
    return subset_df



def GetSubset_newreleased(df: pd.DataFrame, cutoff: Union[str, pd.Timestamp]) -> pd.DataFrame:
    """
    Filter dataframe to include only structures released after a specified date.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with 'PDB_release_date' column.
    cutoff : Union[str, pd.Timestamp]
        Cutoff date in string format (e.g., '2021-10-01') or as pd.Timestamp.
        Structures released AFTER this date are included.
    
    Returns
    -------
    pd.DataFrame
        Filtered dataframe containing only newly released structures (date > cutoff).
    """
    df['PDB_release_date'] = pd.to_datetime(df['PDB_release_date'])
    filtered_df = df[df['PDB_release_date'] > cutoff]
    return filtered_df


def LipidCounts(df: pd.DataFrame) -> None:
    """
    Extract lipid CCD codes from BioDolphinID and print their value counts.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with 'BioDolphinID' column (format: PDB-CHAIN-CHAIN-LIPID_CCD).
    
    Returns
    -------
    None
        Prints lipid CCD codes and their frequency counts to console.
    """
    df['lipid_CCD'] = df['BioDolphinID'].apply(lambda x: x.split('-')[-1][:3])
    lipid_counts = df['lipid_CCD'].value_counts()
    print(lipid_counts.to_string())



def MapMetadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract PDB IDs from BioDolphinID and merge with PDB metadata.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with 'BioDolphinID' column.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with extracted 'pdb_id' column and merged PDB metadata from
        '../PDB_Classifications/pdb_metadata.csv' using left join on 'pdb_id'.
    """
    df = df.copy()
    df['pdb_id'] = df['BioDolphinID'].apply(lambda x: x.split('-')[0][-4:])
    path_metadata = '../PDB_Classifications/pdb_metadata.csv'
    df_metadata = pd.read_csv(path_metadata)
    df_merged = pd.merge(df, df_metadata, on='pdb_id', how='left')
    return df_merged


def load_dataframes(chai_csv_path: str, af_csv_path: str, vina_cvs_path: str, rs_csv_path: str, dd_csv_path: str) -> pd.DataFrame:
    """
    Load RMSD result CSV files for all prediction methods and merge them.
    
    Parameters
    ----------
    chai_csv_path : str
        Path to Chai-1 RMSD results CSV file.
    af_csv_path : str
        Path to AlphaFold 3 RMSD results CSV file.
    vina_cvs_path : str
        Path to AutoDock Vina RMSD results CSV file.
    rs_csv_path : str
        Path to RoseTTAFoldAA RMSD results CSV file.
    dd_csv_path : str
        Path to DiffDock-L RMSD results CSV file.
    
    Returns
    -------
    pd.DataFrame
        Merged dataframe with prefixed columns from all methods.
        Only entries present in ALL files are included (inner join on 'BioDolphinID').
    """
    print('Start adding prefix to the dataframes')
    chai_df = pd.read_csv(chai_csv_path)
    chai_df = AddPrefix(chai_df, 'CHAI')

    af_df = pd.read_csv(af_csv_path)
    af_df = AddPrefix(af_df, 'AF')

    vina_df = pd.read_csv(vina_cvs_path)
    vina_df = AddPrefix(vina_df, 'VINA')

    rs_df = pd.read_csv(rs_csv_path)
    rs_df = AddPrefix(rs_df, 'RS')

    dd_df = pd.read_csv(dd_csv_path)
    dd_df = AddPrefix(dd_df, 'DD')

    dfs = [af_df, chai_df, vina_df, rs_df, dd_df]
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='BioDolphinID', how='inner')
    return merged_df


def load_unique_data(dolphinbusters_path: str, test_data_cutoff: Union[str, pd.Timestamp], unique_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load DolphinBusters metadata and prepare datasets for analysis.
    
    Parameters
    ----------
    dolphinbusters_path : str
        Path to BioDolphinBusters CSV file containing metadata.
    test_data_cutoff : Union[str, pd.Timestamp]
        Cutoff date for separating General vs Untrained datasets.
    unique_df : pd.DataFrame
        Input dataframe with RMSD results and 'BioDolphinID'.
    
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (unique_df, new_released_df) where:
        - unique_df: Full dataset with metadata and PDB classifications
        - new_released_df: Subset released after cutoff date, saved to 'unique_new_rmsd.csv'
    """
    df_dolbuster = pd.read_csv(dolphinbusters_path)
    num_rows = df_dolbuster.shape[0]
    print(f'number of data in the DolphinBusters Dataset evaluated: {num_rows}')
    print('-----------------------------------------------------------')

    unique_df = TagClasses(unique_df, df_dolbuster, add_cols=['BioDolphinID','protein_Sequence'])
    unique_df = MapMetadata(unique_df)
    print(f'Number of  entries: {unique_df.shape[0]}')

    # Newly Released Dataset
    new_released_df = GetSubset_newreleased(unique_df, cutoff=test_data_cutoff)
    print(f'Number of  entries released after {test_data_cutoff}: {new_released_df.shape[0]}')
    new_released_df.to_csv('./unique_new_rmsd.csv', index=False)

    print('-----------------------------------------------------------')
    
    return unique_df, new_released_df


def analyze_success_rates(df_subset: pd.DataFrame, data_subset: str, plot_path: str, f, colors: List[str]) -> None:
    """
    Analyze and visualize success rates across all prediction methods.
    
    Parameters
    ----------
    df_subset : pd.DataFrame
        Dataset to analyze with RMSD columns from all methods.
    data_subset : str
        Label for dataset type ('General' or 'Untrained').
    plot_path : str
        Directory path where plots will be saved.
    f : file object
        Open file handle for writing results summary.
    colors : List[str]
        RGB hex color codes for method visualization (AF, CHAI, VINA, RS, DD).
    
    Returns
    -------
    None
        Generates PNG plots and writes results to file:
        - Bar.png, Bar_new.png, Cumulative.png plots
        - Success rates and 95% CIs written to file handle
    """
    # (1) Success Rate for all data:
    sucessrates_chai,  sucessrates_af, sucessrates_vina, sucessrates_rs,  sucessrates_dd = SucessRates(df_subset, prefix='CHAI'), SucessRates(df_subset, prefix='AF'), SucessRates(df_subset, prefix='VINA'), SucessRates(df_subset, prefix='RS'), SucessRates(df_subset, prefix='DD')
    culsucessrate_chai, culsucessrate_af, culsucessrate_vina, culsucessrate_rs, culsucessrate_dd = Culmulate_SucessRate(df_subset, prefix='CHAI'), Culmulate_SucessRate(df_subset, prefix='AF'), Culmulate_SucessRate(df_subset, prefix='VINA'), Culmulate_SucessRate(df_subset, prefix='RS'), Culmulate_SucessRate(df_subset, prefix='DD')
    plotCumulative(x=culsucessrate_chai[1], y_list=[culsucessrate_af[0], culsucessrate_chai[0], culsucessrate_vina[0], culsucessrate_rs[0], culsucessrate_dd[0]], label_list=['AlphaFold3', 'Chai-1', 'AutoDock Vina', 'RoseTTAFoldAA', 'DiffDock-L'], color_list=colors, filename=f'{plot_path}/Cumulative.png')
    plot_bar_new(categories=['AlphaFold 3', 'Chai-1', 'AutoDock Vina', 'RoseTTAFoldAA', 'DiffDock-L'], values=[sucessrates_af, sucessrates_chai, sucessrates_vina, sucessrates_rs, sucessrates_dd], colors=colors, figname=f'{plot_path}/Bar_new.png')
    
    confi_intervals_chai, confi_intervals_af, confi_intervals_vina, confi_intervals_rs, confi_intervals_dd = GetConfidenceInterval(df_subset, sucessrates_chai[0]), GetConfidenceInterval(df_subset, sucessrates_af[0]), GetConfidenceInterval(df_subset, sucessrates_vina[0]), GetConfidenceInterval(df_subset, sucessrates_rs[0]), GetConfidenceInterval(df_subset, sucessrates_dd[0])
    ci_intervals = [confi_intervals_af, confi_intervals_chai, confi_intervals_vina, confi_intervals_rs, confi_intervals_dd]
    plot_bar(categories=['AlphaFold 3', 'Chai-1', 'AutoDock Vina', 'RoseTTAFoldAA', 'DiffDock-L'], values=[sucessrates_af[0], sucessrates_chai[0], sucessrates_vina[0], sucessrates_rs[0], sucessrates_dd[0]], colors=colors, figname=f'{plot_path}/Bar.png', yerrors=ci_intervals, lipidclass='All')
    f.write(f"# All data: \n")
    f.write(f'Number of data points: {df_subset.shape[0]}\n')
    f.write(f'Chai-1 Success Rate: {sucessrates_chai[0]}, Confidence Interval: {confi_intervals_chai}\n')
    f.write(f'AlphaFold 3 Success Rate: {sucessrates_af[0]}, Confidence Interval: {confi_intervals_af}\n')
    f.write(f'Autodock Vina Success Rate: {sucessrates_vina[0]}, Confidence Interval: {confi_intervals_vina}\n')
    f.write(f'RoseTTAFoldAA Success Rate: {sucessrates_rs[0]}, Confidence Interval: {confi_intervals_rs}\n')
    f.write(f'DiffDock-L Success Rate: {sucessrates_dd[0]}, Confidence Interval: {confi_intervals_dd}\n\n')


def analyze_lipid_classes(df_subset: pd.DataFrame, plot_path: str, f, colors: List[str]) -> None:
    """
    Analyze success rates stratified by lipid classification.
    
    Parameters
    ----------
    df_subset : pd.DataFrame
        Dataset to analyze with lipid classification columns.
    plot_path : str
        Directory path where plots will be saved.
    f : file object
        Open file handle for writing results summary.
    colors : List[str]
        RGB hex color codes for method visualization.
    
    Returns
    -------
    None
        Generates PNG plots for each of 8 lipid classes:
        - Bar_LipidClass.png: Comprehensive comparison
        - Bar_*.png: Individual plots for each lipid class
    """
    lipid_classes = ['Glycerophospholipid', 'Fatty_Acyl', 'Glycerolipid',  'Polyketide', 'Prenol_Lipid', 'Saccharolipid', 'Sphingolipid', 'Sterol_Lipid']
    sucessrates_lipid_dict = {}
    ci_intervals_lipid_dict = {}
    num_datapoints_list = []

    for lipid_class in lipid_classes:
        df_subset_lipid = GetSubset_LipidClass(df_subset, subset=lipid_class)
        sucessrate_chai,  sucessrate_af, sucessrate_vina, sucessrate_rs,  sucessrate_dd = SucessRates(df_subset_lipid, prefix='CHAI')[0], SucessRates(df_subset_lipid, prefix='AF')[0], SucessRates(df_subset_lipid, prefix='VINA')[0], SucessRates(df_subset_lipid, prefix='RS')[0], SucessRates(df_subset_lipid, prefix='DD')[0]
        confi_intervals_chai, confi_intervals_af, confi_intervals_vina, confi_intervals_rs, confi_intervals_dd = GetConfidenceInterval(df_subset_lipid, sucessrate_chai), GetConfidenceInterval(df_subset_lipid, sucessrate_af), GetConfidenceInterval(df_subset_lipid, sucessrate_vina), GetConfidenceInterval(df_subset_lipid, sucessrate_rs), GetConfidenceInterval(df_subset_lipid, sucessrate_dd)
        ci_intervals = [confi_intervals_af, confi_intervals_chai, confi_intervals_vina, confi_intervals_rs, confi_intervals_dd]
        plot_bar(categories=['AlphaFold 3', 'Chai-1', 'AutoDock Vina', 'RoseTTAFoldAA', 'DiffDock-L'], values=[sucessrate_af, sucessrate_chai, sucessrate_vina, sucessrate_rs, sucessrate_dd], colors=colors, figname=f'{plot_path}/Bar_{lipid_class}.png', lipidclass=lipid_class, yerrors=ci_intervals)
        sucessrates_lipid_dict[lipid_class] = [sucessrate_af, sucessrate_chai, sucessrate_vina, sucessrate_rs, sucessrate_dd]
        ci_intervals_lipid_dict[lipid_class] = ci_intervals
        num_data = df_subset_lipid.shape[0]
        num_datapoints_list.append(num_data)
    plot_bar_lipid(sucessrates_lipid_dict, ci_intervals_lipid_dict, num_datapoints_list, colors=colors, figname=f'{plot_path}/Bar_LipidClass.png')


def analyze_molecular_weights(df_subset: pd.DataFrame, plot_path: str, f, colors: List[str]) -> None:
    """
    Analyze success rates stratified by molecular weight ranges.
    
    Parameters
    ----------
    df_subset : pd.DataFrame
        Dataset to analyze with 'lipid_Molecular_weight' column.
    plot_path : str
        Directory path where plots will be saved.
    f : file object
        Open file handle for writing results summary.
    colors : List[str]
        RGB hex color codes for method visualization.
    
    Returns
    -------
    None
        Generates PNG plot:
        - Bar_MW.png: Success rates across 20 MW bins (0-2000 Da in 100 Da steps)
    """
    sucessrates_MW_dict = {}
    ci_intervals_MW_dict = {}
    num_datapoints_list = []
    for start in range(0, 2001, 100):
        end = start + 100
        unique_df_subset_MW = GetSubset_LipidMW(df_subset, range_start=start, range_end=end)
        num_data = unique_df_subset_MW.shape[0]
        if not num_data == 0:
            sucessrate_chai,  sucessrate_af, sucessrate_vina, sucessrate_rs,  sucessrate_dd = SucessRates(unique_df_subset_MW, prefix='CHAI')[0], SucessRates(unique_df_subset_MW, prefix='AF')[0], SucessRates(unique_df_subset_MW, prefix='VINA')[0], SucessRates(unique_df_subset_MW, prefix='RS')[0], SucessRates(unique_df_subset_MW, prefix='DD')[0]
            sucessrates_MW_dict[f'{start}-{end}'] = [sucessrate_af, sucessrate_chai, sucessrate_vina, sucessrate_rs, sucessrate_dd]
            confi_intervals_chai, confi_intervals_af, confi_intervals_vina, confi_intervals_rs, confi_intervals_dd = GetConfidenceInterval(unique_df_subset_MW, sucessrate_chai), GetConfidenceInterval(unique_df_subset_MW, sucessrate_af), GetConfidenceInterval(unique_df_subset_MW, sucessrate_vina), GetConfidenceInterval(unique_df_subset_MW, sucessrate_rs), GetConfidenceInterval(unique_df_subset_MW, sucessrate_dd)
            ci_intervals = [confi_intervals_af, confi_intervals_chai, confi_intervals_vina, confi_intervals_rs, confi_intervals_dd]
            ci_intervals_MW_dict[f'{start}-{end}'] = ci_intervals

            num_datapoints_list.append(num_data)
    # Set without intervals for now
    ci_intervals_MW_dict = None
    plot_bar_MW(sucessrates_MW_dict, ci_intervals_MW_dict, num_datapoints_list, colors=colors, figname=f'{plot_path}/Bar_MW.png')


def analyze_protein_errors_and_scoring(df_subset: pd.DataFrame, plot_path: str, f, colors: List[str]) -> None:
    """
    Analyze protein structure prediction errors and method scoring power.
    
    Parameters
    ----------
    df_subset : pd.DataFrame
        Dataset to analyze with protein error and lipid RMSD columns.
    plot_path : str
        Directory path where plots will be saved.
    f : file object
        Open file handle for writing results summary.
    colors : List[str]
        RGB hex color codes for method visualization.
    
    Returns
    -------
    None
        Generates plots analyzing protein errors (AF, CHAI, RS) and scoring correlations:
        - Scatter plots of protein errors vs lipid RMSD
        - Bar plots comparing protein prediction errors
        - Scoring power analysis plots
    """
    # (4) Effects on protein structure predicion errors
    plot_scatter_protein_errors_seaborn(df_subset, prefix_list=['AF', 'CHAI', 'RS'], fig_path=plot_path)
    plot_bar_protein_errors_seaborn(df_subset, prefix_list=['AF', 'CHAI', 'RS'], fig_path=plot_path)

    # (5) Scoring power of all methods 
    plot_scoring_power(df_subset, fig_path=plot_path)


def run_analysis(unique_df: pd.DataFrame, new_released_df: pd.DataFrame, test_data_cutoff: Union[str, pd.Timestamp]) -> None:
    """
    Execute complete comparative analysis pipeline for both General and Untrained datasets.
    
    Parameters
    ----------
    unique_df : pd.DataFrame
        Full dataset with all methods' results and metadata.
    new_released_df : pd.DataFrame
        Subset of structures released after test_data_cutoff.
    test_data_cutoff : Union[str, pd.Timestamp]
        Date cutoff for separating datasets (displayed in results).
    
    Returns
    -------
    None
        Generates comprehensive analysis outputs:
        - result_successrate.txt: Summary of metrics, success rates, and CIs
        - Results_General/: Plots for full dataset
        - Results_Untrained/: Plots for newly released structures
        - Comparison plots for 5 methods across multiple analysis categories
    """
    print('Start analyzing and ploting')
    colors = ['#a1c9f4', '#ffb482', '#8de5a1', '#ff9f9b', '#c9a0dc']

    with open("result_successrate.txt", "w", encoding="utf-8") as f:

        for data_subset in ['General', 'Untrained']:
            if data_subset == 'General':
                print('Analyzing General Data')
                f.write("---------General Dataset:----------\n")
                df_subset = unique_df
                plot_path = './Results_General'
                os.makedirs(plot_path, exist_ok=True)
                
            elif data_subset == 'Untrained':
                print('Analyzing Untrained Data')
                f.write(f"\n---------Untrained Dataset (with {test_data_cutoff} as cutoff ): ----------\n")
                df_subset = new_released_df
                plot_path = './Results_Untrained'
                os.makedirs(plot_path, exist_ok=True)

            # Run analysis functions
            analyze_success_rates(df_subset, data_subset, plot_path, f, colors)
            analyze_lipid_classes(df_subset, plot_path, f, colors)
            analyze_molecular_weights(df_subset, plot_path, f, colors)
            analyze_protein_errors_and_scoring(df_subset, plot_path, f, colors)

            print(f'Finished plotting for {data_subset} Data')


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()

    config = load_config()
    dolphinbusters_path = config["dolphinbusters_path"]
    print(f'Using dolphinbusters_path: {dolphinbusters_path}')

    TEST_DATA_CUTOFF = '2021-10-01' 

    print('Start running the program:')

    print('Using the provided unique_rmsd.csv file')
    unique_df = pd.read_csv('./unique_rmsd.csv')
    unique_df, new_released_df = load_unique_data(dolphinbusters_path, TEST_DATA_CUTOFF, unique_df)

    # Run analysis
    run_analysis(unique_df, new_released_df, TEST_DATA_CUTOFF)
            
        





