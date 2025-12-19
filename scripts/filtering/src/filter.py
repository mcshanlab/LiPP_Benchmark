import os
from pathlib import Path
from typing import List, Tuple, Any

import pandas as pd
from Bio import PDB
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem

from src.pairwisedistances import *



def Filter_Resolution(df: pd.DataFrame, cutoff: float = 2) -> pd.DataFrame:
    """Filter structures by X-ray resolution threshold.

    Args:
        df: Input dataframe with 'complex_Resolution' column.
        cutoff: Maximum resolution (Angstroms) to keep; default 2.

    Returns:
        pd.DataFrame: Dataframe with entries between 0 and cutoff resolution.
    """
    new_df = df[df['complex_Resolution'] < cutoff]
    new_df = new_df[new_df['complex_Resolution'] > 0]
    return new_df

def Filter_Unk_Residue(df: pd.DataFrame) -> pd.DataFrame:
    """Remove structures with unknown amino acid residues in protein.

    Args:
        df: Input dataframe with 'protein_Sequence' column.

    Returns:
        pd.DataFrame: Dataframe without entries containing 'UNK' in protein sequence.
    """
    new_df = df[~df['protein_Sequence'].str.contains('UNK')]
    return new_df

def Filter_Mod_Residue(df: pd.DataFrame) -> pd.DataFrame:
    """Remove structures with modified amino acid residues in protein.

    Args:
        df: Input dataframe with 'protein_Sequence' column.

    Returns:
        pd.DataFrame: Dataframe without entries containing '(' (mod indicator) in protein sequence.
    """
    new_df = df[~df['protein_Sequence'].str.contains('(', regex=False)]
    return new_df

def Filter_Unk_Atom(df: pd.DataFrame, pdb_path: str) -> pd.DataFrame:
    """Remove structures with unknown atoms detected via PDB inspection.

    Args:
        df: Input dataframe with 'BioDolphinID' column.
        pdb_path: Path to the directory containing PDB files.

    Returns:
        pd.DataFrame: Dataframe without entries containing unknown atoms.
    """
    df['HAS_UNKNOWN_ATOM'] = df['BioDolphinID'].apply(lambda bdid: Identify_UnkAtom(bdid, pdb_path))
    new_df = df[df['HAS_UNKNOWN_ATOM'] == False]
    return new_df
        

def Identify_UnkAtom(BDID: str, pdb_path: str) -> bool:
    """Check if a PDB structure contains unknown atoms (marked as 'X').

    Args:
        BDID: BioDolphin ID used to construct the PDB file path.

    Returns:
        bool: True if unknown atoms detected, False otherwise.
    """
    #pdb_path = "/storage/cedar/cedar0/cedarp-amcshan3-0/plesk_data/entry_selected/pdbs_selected"
    pdb_file = f'{pdb_path}/{BDID}.pdb'
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("curr_struct", pdb_file)
    HAS_UNKNOWN = False
    for atom in structure.get_atoms():
        if atom.get_name() == "X": #X1
            HAS_UNKNOWN = True
            print(f'detect unknown atom from {BDID}.pdb')
            return HAS_UNKNOWN
    return HAS_UNKNOWN


def Filter_PL_distance(df: pd.DataFrame, pdb_original_path: str) -> pd.DataFrame:
    """Filter by protein-lipid distance constraints using PyMOL.

    Args:
        df: Input dataframe with 'complex_PDB_ID', 'complex_Receptor_Chain',
            'complex_Ligand_Chain', 'complex_Residue_number_of_the_ligand'.
        pdb_original_path: Path to the directory containing original PDB files.

    Returns:
        pd.DataFrame: Dataframe with valid protein-lipid distances (0.2-5 Å).
    """
    df['DIST_IN_RANGE'] = df.apply(lambda row: Calc_dist_range(row, pdb_original_path), axis=1)
    new_df = df[df['DIST_IN_RANGE'] == True]
    return new_df

def Calc_dist_range(row: pd.Series, pdb_original_path: str) -> bool:
    """Extract distance parameters from dataframe row and test distance.

    Args:
        row: Dataframe row with PDB ID, chain, and residue info.
        pdb_original_path: Path to the directory containing original PDB files.

    Returns:
        bool: True if distance constraints are satisfied.
    """
    pdbid, proteinChain, lipidChain, lipidResNum = row['complex_PDB_ID'], row['complex_Receptor_Chain'], row['complex_Ligand_Chain'], row['complex_Residue_number_of_the_ligand']
    return test_distance(pdbid, proteinChain, lipidChain, lipidResNum, pdb_original_path)
        
        
def test_distance(pdbid: str, proteinChain: str, lipidChain: str, lipidResNum: str, pdb_original_path: str) -> bool:
    """Test if protein-lipid distance constraints are satisfied using PyMOL.

    Checks that:
    - Target protein is within 5 Å of target lipid.
    - Min distance between them is >= 0.2 Å.
    - Other protein chains are >5 Å from target lipid.

    Args:
        pdbid: PDB ID.
        proteinChain: Chain ID of target protein.
        lipidChain: Chain ID of target lipid.
        lipidResNum: Residue number of lipid in its chain.
        pdb_original_path: Path to the directory containing original PDB files.

    Returns:
        bool: True if all constraints satisfied.
    """
    print('------------------------------------------------')
    cmd.delete("all")
    pdbpath = pdb_original_path
    pdb = f'{pdbpath}/pdb{pdbid}.ent'
    if os.path.exists(pdb):
        print(f'pdb path: {pdb}')
    else:
        pdb = f'{pdbpath}/{pdbid}.cif'
        print(f'cif path: {pdb}')
    cmd.load(pdb)
    KEEP = True

    # Step1: Check the target protein and lipid distance first
    # select protein
    selectname, selestr = 'protein_sele', f'chain {proteinChain} and polymer'
    print(f'select name: {selectname}')
    print(f'select string: {selestr}')
    cmd.select(selectname, selestr)
    # select lipid
    selectname, selestr = 'lipid_sele', f'chain {lipidChain} and resi {lipidResNum}'
    cmd.select(selectname, selestr)
    print(f'select name: {selectname}')
    print(f'select string: {selestr}')
    # calculate distance
    cutoff='5'
    dist_list = pairwise_dist('protein_sele','lipid_sele',cutoff,output='N',show='N')
    # check the distance
    if len(dist_list) == 0: # if no protein atoms is within 5 A from the lipid, then reject this pair
        print('### not keeping because target protein is too far away from (>5 A) from target lipid')
        KEEP = False
        cmd.deselect()
        return KEEP
    min_dist = min(dist_list)
    if min_dist < 0.2: # if there is a distance less than 0.2 between the lipid and the protein, then reject this pair
        print('### not keeping because target lipid and target protein is too close (<0.2 A)')
        KEEP = False
        cmd.deselect()
        return KEEP
    
    # Step2: Check the other proteins and lipid distance next -> reject if any proteins in the pdb are too close to the lipid
    chains_list = cmd.get_chains()
    chains_list.remove(proteinChain)
    if len(chains_list) ==0:
        print("no other chains detected")
        cmd.deselect()
        print(f'Keeping= {KEEP}')
        return KEEP
    chain_selection = " or ".join([f"chain {chain}" for chain in chains_list])
    selection_string = f"({chain_selection}) and polymer"
    cmd.select('other_proteins', selection_string)
    dist_list = pairwise_dist('other_proteins','lipid_sele',cutoff,output='N',show='N')

    if (len(dist_list) != 0):
        print('### not keeping because there are other protein partners that are close to the target lipid')
        KEEP = False
    else:
        print(f'Keeping= {KEEP}')
    cmd.deselect()
    return KEEP

def Filter_RDKIT(df: pd.DataFrame, biodolphin_pdb_path: str) -> pd.DataFrame:
    """Filter lipids by RDKIT sanitization and completeness check.

    Args:
        df: Input dataframe with 'BioDolphinID' and 'lipid_Canonical_smiles'.

    Returns:
        pd.DataFrame: Dataframe with only complete, valid lipid structures.
    """
    df['PASS_RDKIT'] = df.apply(lambda row: CheckRDKIT(row, biodolphin_pdb_path), axis=1)
    new_df = df[df['PASS_RDKIT'] == True]
    return new_df
    
def CheckRDKIT(row: pd.Series, biodolphin_pdb_path: str) -> bool:
    """Extract ID and SMILES from row and perform completeness check.

    Args:
        row: Dataframe row with 'BioDolphinID' and 'lipid_Canonical_smiles'.

    Returns:
        bool: True if lipid is complete and valid.
    """
    BDID, lipid_SMILES = row['BioDolphinID'], row['lipid_Canonical_smiles']
    return CheckComplete(BDID, lipid_SMILES, biodolphin_pdb_path)

def CheckComplete(BDID: str, lipid_SMILES: str, biodolphin_pdb_path: str) -> bool:
    """Validate lipid structure completeness via SMILES and PDB alignment.

    Args:
        BDID: BioDolphin ID to construct lipid PDB file path.
        lipid_SMILES: Canonical SMILES string for the lipid.

    Returns:
        bool: True if PDB structure matches SMILES template.
    """
    PASS_COMPLETE = True
    pdb_lipid_file = f'{biodolphin_pdb_path}/{BDID}/lipid.pdb'
    #pdb_lipid_file = f'/storage/cedar/cedar0/cedarp-amcshan3-0/docking_inputs/full_pdb_inputs/{BDID}/lipid.pdb'
    print(f'checking: {pdb_lipid_file}')


    # Read the pdb file of the target lipid
    try:
        # Define complete ligand from CCD data using SMILES
        std_mol = Chem.MolFromSmiles(lipid_SMILES)
        pdb_mol = Chem.MolFromPDBFile(pdb_lipid_file, proximityBonding=True, sanitize=True, removeHs=False)
        pdb_mol_fix = AllChem.AssignBondOrdersFromTemplate(std_mol, pdb_mol)
    except Exception as e:
        PASS_COMPLETE = False
        print(f'Did not pass. {e}')
        return PASS_COMPLETE

    if (pdb_mol_fix.HasSubstructMatch(std_mol)) and (std_mol.HasSubstructMatch(pdb_mol_fix)):
        PASS_COMPLETE = True
        print('Pass RDKIT!!!')
    else:
        PASS_COMPLETE = False
        print('Not the same with original ccd smiles')

    return PASS_COMPLETE
    

# Filter: lipid that is complete, and has no sterochemical errors or clashes

def Filter_Lipid_Covalent(df: pd.DataFrame) -> pd.DataFrame:
    """Remove lipids with covalent bonds to protein.

    Args:
        df: Input dataframe with 'Ligand_Covalent' boolean column.

    Returns:
        pd.DataFrame: Dataframe with non-covalent lipids only.
    """
    new_df = df.loc[df['Ligand_Covalent'] != True]
    return new_df

def Filter_Lipid_Clash(df: pd.DataFrame) -> pd.DataFrame:
    """Remove lipids with atomic clashes.

    Args:
        df: Input dataframe with 'Ligand_Num_Clash' integer column.

    Returns:
        pd.DataFrame: Dataframe with clash-free lipids (count == 0).
    """
    new_df = df.loc[df['Ligand_Num_Clash'] == 0]
    return new_df

def Filter_Lipid_Stereo(df: pd.DataFrame) -> pd.DataFrame:
    """Remove lipids with stereochemical errors.

    Args:
        df: Input dataframe with 'Ligand_Stereo' integer column.

    Returns:
        pd.DataFrame: Dataframe with no stereo errors (count == 0).
    """
    new_df = df.loc[df['Ligand_Stereo'] == 0]
    return new_df


def Filter_Lipid_Complete(df: pd.DataFrame) -> pd.DataFrame:
    """Remove lipids with missing atoms.

    Args:
        df: Input dataframe with 'Ligand_Completeness' binary column.

    Returns:
        pd.DataFrame: Dataframe with complete lipids (completeness == 1).
    """
    new_df = df.loc[df['Ligand_Completeness'] == 1]
    return new_df









