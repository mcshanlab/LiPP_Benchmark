"""
Run PoseBusters validation on protein-ligand structures.

Uses the PoseBusters suite to validate the chemical validity of docked ligands
against reference structures.
"""

import argparse
from pathlib import Path

import pandas as pd
from posebusters import PoseBusters



def run_buster_table(csv_path):
    """
    Validate all molecules in a CSV table using PoseBusters.

    Args:
        csv_path (str): Path to CSV file containing molecules to validate.

    Returns:
        pd.DataFrame: Full validation report with all checks.
    """
    buster = PoseBusters(config="redock")
    df = pd.read_csv(csv_path)
    df_result = buster.bust_table(mol_table=df, full_report=True)
    return df_result


def run_buster(pred_sdf, true_sdf, cond_pdb):
    """
    Validate a single set of predicted and reference structures.

    Args:
        pred_sdf (str): Path to predicted ligand structure (SDF format).
        true_sdf (str): Path to reference/true ligand structure (SDF format).
        cond_pdb (str): Path to protein structure (PDB format).

    Returns:
        pd.DataFrame: Validation report for the structure pair.
    """
    pred_file = Path(pred_sdf)
    true_file = Path(true_sdf)
    cond_file = Path(cond_pdb)

    buster = PoseBusters(config="redock")
    df = buster.bust([pred_file], true_file, cond_file)
    return df






def main():
    """Run PoseBusters validation and generate reports."""
    parser = argparse.ArgumentParser(
        description="Validate protein-ligand docking predictions using PoseBusters"
    )
    parser.add_argument(
        "-c",
        "--csv-file",
        type=str,
        help="Path to input CSV file containing molecules to validate",
        required=True,
    )

    args = parser.parse_args()
    csv_file = args.csv_file

    # Define columns for condensed report
    short_report_cols = [
        "mol_pred_loaded",
        "mol_true_loaded",
        "mol_cond_loaded",
        "sanitization",
        "inchi_convertible",
        "all_atoms_connected",
        "molecular_formula",
        "molecular_bonds",
        "double_bond_stereochemistry",
        "tetrahedral_chirality",
        "bond_lengths",
        "bond_angles",
        "internal_steric_clash",
        "aromatic_ring_flatness",
        "non-aromatic_ring_non-flatness",
        "double_bond_flatness",
        "internal_energy",
        "protein-ligand_maximum_distance",
        "minimum_distance_to_protein",
        "minimum_distance_to_organic_cofactors",
        "minimum_distance_to_inorganic_cofactors",
        "minimum_distance_to_waters",
        "volume_overlap_with_protein",
        "volume_overlap_with_organic_cofactors",
        "volume_overlap_with_inorganic_cofactors",
        "volume_overlap_with_waters",
        "rmsd_≤_2å",
    ]

    # Run PoseBusters validation
    df_result_long = run_buster_table(csv_file)

    # Save condensed report
    df_result_short = df_result_long[short_report_cols]
    # Generate output filename based on input CSV filename
    input_path = Path(csv_file)
    output_file = f"{input_path.stem}_results_short.csv"
    df_result_short.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
    
    

    








