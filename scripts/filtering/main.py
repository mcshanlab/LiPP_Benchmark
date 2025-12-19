import argparse
from typing import Optional
import os

import pandas as pd
#import sys
#sys.path.append('./src')
from src.filter import *
from src.parsePDB import *
from src.concat import *
from src.stats import *

import warnings
warnings.filterwarnings('ignore')


def ReadDataset(dataset_path: str) -> pd.DataFrame:
    """Read BioDolphin dataset from txt/csv.

    Args:
        dataset_path: Path to .txt (tab-separated) or .csv file.

    Returns:
        pd.DataFrame: Loaded dataset contents.
    """
    if dataset_path.endswith(".txt"):
        df = pd.read_csv(dataset_path, sep='\t')
    elif dataset_path.endswith(".csv"):
        df = pd.read_csv(dataset_path)
    else:
        raise Exception("invalid biodolphin dataset format")
    return df
    
def CalculateSize(dataframe: pd.DataFrame) -> int:
    """Return number of rows in dataframe.

    Args:
        dataframe: Input dataframe whose row count is needed.

    Returns:
        int: Number of rows.
    """
    return dataframe.shape[0]


def GetPDB_info(df_result: pd.DataFrame, savePATH: str) -> pd.DataFrame:
    """Augment dataset with PDB metadata and lipid quality info.

    Args:
        df_result: Dataframe of complexes to annotate.
        savePATH: Directory to write intermediate CSVs.

    Returns:
        pd.DataFrame: Input dataframe with added PDB info columns.
    """
    PDB_list = list(set(df_result['complex_PDB_ID'].tolist()))
    df_date, df_protein, df_lipid = GetPDBdata(PDB_list, savePATH)

    # tag the df_result with dates
    print(f'shape of the dataset: {df_result.shape}')

    df_date = df_date.rename(columns={'PDBID': 'complex_PDB_ID'})
    df_date = df_date.drop_duplicates(subset='complex_PDB_ID')
    df_result = pd.merge(df_result, df_date, on='complex_PDB_ID', how='left')
    print(f'shape of the dataset: {df_result.shape}')

    # tag the df_result with lipid quality info
    df_lipid = df_lipid.rename(columns={'PDBID': 'complex_PDB_ID', 'Ligand_resNum':'complex_Residue_number_of_the_ligand', 'Ligand_Chain':'complex_Ligand_Chain', 'Ligand_CompID':'lipid_Ligand_ID_CCD'})
    
    df_lipid['complex_Residue_number_of_the_ligand']=df_lipid['complex_Residue_number_of_the_ligand'].astype(str)
    df_result['complex_Residue_number_of_the_ligand']=df_result['complex_Residue_number_of_the_ligand'].astype(str)

    df_result = pd.merge(df_result, df_lipid, on=['complex_PDB_ID', 'complex_Residue_number_of_the_ligand', 'complex_Ligand_Chain', 'lipid_Ligand_ID_CCD'], how='left')
    #print(f'shape of the dataset: {df_result.shape}')
    df_result.to_csv(f'{savePATH}/df_PDBinfo.csv')

    return df_result


def process_first_filters(BD_dataset: str, PATH_INTERMEDIATE: str, res_cutoff: float, biodolphin_pdb_path: str) -> pd.DataFrame:
    """Run the initial filtering pipeline and persist intermediates.

    Args:
        BD_dataset: Path to the BioDolphin dataset file.
        PATH_INTERMEDIATE: Directory to store intermediate CSVs and logs.
        res_cutoff: Resolution cutoff to filter structures.
        biodolphin_pdb_path: Path to the directory containing BioDolphin PDB files.

    Returns:
        pd.DataFrame: Filtered dataframe after early-stage checks.
    """
    print('Processing the first filters')
    with open(f"{PATH_INTERMEDIATE}/outfile.txt", "a") as f:
        print(f'Reading the dataset: {BD_dataset}')
        df = ReadDataset(BD_dataset)

        num_data = CalculateSize(df)
        print(f'original number of data: {num_data}')
        f.write(f'original number of data: {num_data} \n')

        df_result = df

        df_result = Filter_Unk_Residue(df_result)
        print(f'# Number of data after filtering unknown protein residue: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering unknown protein residue: {CalculateSize(df_result)} \n')

        df_result = Filter_Mod_Residue(df_result)
        print(f'# Number of data after filtering modified residue: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering modified residue: {CalculateSize(df_result)} \n')

        df_result = Filter_Resolution(df_result, cutoff=res_cutoff)
        print(f'# Number of data after filtering resolution: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering resolution: {CalculateSize(df_result)} \n')

        #df_result = Filter_RDKIT(df_result, biodolphin_pdb_path)
        print(f'# Number of data after filtering rdkit check: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering rdkit check: {CalculateSize(df_result)} \n')
        df_result.to_csv(f"{PATH_INTERMEDIATE}/df_filter_rdkit.csv")

        df_result = GetPDB_info(ReadDataset(f"{PATH_INTERMEDIATE}/df_filter_rdkit.csv"), PATH_INTERMEDIATE)

        df_result = Filter_Lipid_Covalent(df_result)
        print(f'# Number of data after filtering lipid covalent bond: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering lipid covalent bond: {CalculateSize(df_result)} \n')
        df_result.to_csv(f"{PATH_INTERMEDIATE}/df_filter_covalent.csv")

        df_result = Filter_Lipid_Clash(df_result)
        print(f'# Number of data after filtering lipid clashes: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering lipid clashes: {CalculateSize(df_result)} \n')
        df_result.to_csv(f"{PATH_INTERMEDIATE}/df_filter_lipidclash.csv")

        df_result = Filter_Lipid_Stereo(df_result)
        print(f'# Number of data after filtering lipid stereochemical errors: {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering lipid stereochemical errors: {CalculateSize(df_result)} \n')
        df_result.to_csv(f"{PATH_INTERMEDIATE}/df_filter_lipidstereo.csv")

        df_result = Filter_Lipid_Complete(df_result)
        print(f'# Number of data after filtering lipid completness : {CalculateSize(df_result)}')
        f.write(f'# Number of data after filtering lipid completness : {CalculateSize(df_result)} \n')
        df_result.to_csv(f"{PATH_INTERMEDIATE}/df_filtered_lipidcomplete.csv")

    return df_result


def process_last_filter(PATH_INTERMEDIATE: str, subset_IDX: int, subindex: Optional[int] = None, original_pdb_path: str = "") -> pd.DataFrame:
    """Apply final distance filter to a subset and save the result.

    Args:
        PATH_INTERMEDIATE: Directory containing intermediate CSVs.
        subset_IDX: 1-based subset index determining the 100-row slice.
        subindex: Optional extra index (unused placeholder for compatibility).
        original_pdb_path: Path to the directory containing original PDB files.

    Returns:
        pd.DataFrame: Distance-filtered subset dataframe.
    """
    print('Processing the last filter')

    df_result = ReadDataset(f"{PATH_INTERMEDIATE}/df_filtered_lipidcomplete.csv")
    print(f'Read input dataset size: {CalculateSize(df_result)}')

    START_IDX = (subset_IDX-1) * 100 
    END_IDX = START_IDX + 100
    SAVE_FILE = f'df_filter_distance_{subset_IDX}.csv'

    if END_IDX > (CalculateSize(df_result)):
        df_result_subset = df_result.iloc[START_IDX:]
        print(f'process subset as from {START_IDX} to end of the dataframe')
    else:
        df_result_subset = df_result.iloc[START_IDX:END_IDX]
        print(f'START_IDX: {START_IDX}')
        print(f'END_IDX: {END_IDX}')

    print(f'Subset size: {CalculateSize(df_result_subset)}')

    print(df_result_subset.columns)

    df_result_subset = Filter_PL_distance(df_result_subset, pdb_original_path=original_pdb_path)
    print(f'# Number of data in subset after filtering distance: {CalculateSize(df_result_subset)}')

    df_result_subset = df_result_subset.loc[:, ~df_result_subset.columns.str.contains('^Unnamed')]
    df_result_subset.to_csv(f"{PATH_INTERMEDIATE}/{SAVE_FILE}")

    return df_result_subset


def combine_dataframes(PATH_INTERMEDIATE: str, resolution: int = 2) -> pd.DataFrame:
    """Combine filtered subsets and compute final stats.

    Args:
        PATH_INTERMEDIATE: Directory containing filtered subset CSVs.
        resolution: Resolution cutoff to apply when combining.

    Returns:
        pd.DataFrame: Final combined dataframe after all filters.
    """
    print('combining the dataframes')
    with open(f"{PATH_INTERMEDIATE}/outfile.txt", "a") as f:
        final_df = getFinal(repo_PATH='./', intermediate_PATH=PATH_INTERMEDIATE, final_fn='DolphinBuster_res2.csv', res=resolution )
        print(f'# Number of data after filtering distance and with resoltion < {resolution}: {CalculateSize(final_df)}')
        f.write(f'# Number of data after filtering distance and with resoltion < {resolution}: {CalculateSize(final_df)} \n')

        PrintStats(final_df, f)

    return final_df


def main() -> None:
    """CLI entrypoint for running filtering stages or combining results.

    Parses arguments and dispatches to the requested pipeline stage.

    Returns:
        None
    """
    print('Starting the program')
    parser = argparse.ArgumentParser(description="A simple script to process a file.")
    parser.add_argument('-d','--dataset', type=str, help="The full path of the dataset file.", default='../../data/BioDolphin_vr1.1.csv')
    parser.add_argument('--biodolphin_pdb_path', type=str, help="The full path of the pdb files of BioDolphin.", default='/storage/cedar/cedar0/cedarp-amcshan3-0/plesk_data/entry_selected/pdbs_selected')
    parser.add_argument('--original_pdb_path', type=str, help="The full path of the original pdb files in BioDolphin.", default='/storage/cedar/cedar0/cedarp-amcshan3-0/plesk_data/pdbs_original')

    parser.add_argument('--last_filter', action='store_true', help="Process the last filter")
    parser.add_argument('-s', '--subset',type=int, help="The index of the subset to process")
    parser.add_argument('-i', '--subindex',type=int, default=None, help="The subindex of the subset to process")
    parser.add_argument('--combine', action='store_true', help="Combine the dataframes")

    args = parser.parse_args()
    BD_dataset = args.dataset
    biodolphin_pdb_path = args.biodolphin_pdb_path
    original_pdb_path = args.original_pdb_path

    PATH_INTERMEDIATE = './intermediate_files'
    os.makedirs(PATH_INTERMEDIATE, exist_ok=True)
    res_cutoff = 2
    print(f'PATH_INTERMEDIATE is set to {PATH_INTERMEDIATE}')
    print(f'res_cutoff is set to {res_cutoff}')

    if (not args.last_filter) and (not args.combine):
        process_first_filters(BD_dataset, PATH_INTERMEDIATE, res_cutoff, biodolphin_pdb_path)

    if args.last_filter:
        process_last_filter(PATH_INTERMEDIATE, args.subset, args.subindex, original_pdb_path)

    if args.combine:
        combine_dataframes(PATH_INTERMEDIATE, resolution=2)


if __name__ == "__main__":
    main()

