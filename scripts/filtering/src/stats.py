from typing import TextIO

import pandas as pd


def PrintStats(df: pd.DataFrame, writefile: TextIO) -> None:
    """Compute and print/write dataset statistics to console and file.

    Generates summary statistics including unique CCD identifiers, lipid
    categories, PDB structures, and UniProt IDs.

    Args:
        df: Input dataframe with columns 'lipid_Ligand_ID_CCD',
            'lipid_Lipidmaps_categories', 'complex_PDB_ID', 'protein_UniProt_ID'.
        writefile: File object (opened in write/append mode) to write statistics.

    Returns:
        None: Prints to console and writes to file as side effects.
    """
    writefile.write(f'--------------- final statistics ------------------- \n')

    # unique CCD
    unique_CCD_list = df['lipid_Ligand_ID_CCD'].unique().tolist()
    print(f'unique number of CCD: {len(unique_CCD_list)}')
    writefile.write(f'# unique number of CCD: {len(unique_CCD_list)} \n')

    # lipid categories
    print(df['lipid_Lipidmaps_categories'].value_counts())
    writefile.write(f'{df['lipid_Lipidmaps_categories'].value_counts()} \n')

    # unique PDB IDs
    unique_PDB_list = df['complex_PDB_ID'].unique().tolist()
    print(f'unique number of PDB IDs: {len(unique_PDB_list)}')
    writefile.write(f'unique number of PDB IDs: {len(unique_PDB_list)} \n')

    # unique uniprotIDs
    unique_UniProt_list = df['protein_UniProt_ID'].unique().tolist()
    print(f'unique number of UniProtIDs: {len(unique_UniProt_list)}')
    writefile.write(f'unique number of UniProtIDs: {len(unique_UniProt_list)} \n')

    return None


