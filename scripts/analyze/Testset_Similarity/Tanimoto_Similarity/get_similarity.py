import pandas as pd
from functools import lru_cache
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def read_data():
    all_df = pd.read_csv("../ProteinCartography/demo/cluster_benchmark/BioDolphin_vr1.1.csv")
    test_df = pd.read_csv("../ProteinCartography/demo/cluster_benchmark/test_set.csv")

    CCD_all = (
        all_df["lipid_Ligand_ID_CCD"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s != ""]
        .drop_duplicates()
        .tolist()
    )
    CCD_test = (
        test_df["lipid_Ligand_ID_CCD"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s != ""]
        .drop_duplicates()
        .tolist()
    )

    return CCD_all, CCD_test


@lru_cache(maxsize=None)
def mol_from_ccd(ccd_id: str) -> Chem.Mol:
    """
    Fetch a PDB CCD ligand by its 3-letter (sometimes longer) CCD ID and return an RDKit Mol.
    Uses the 'ideal' SDF download from RCSB.
    """
    ccd_id = ccd_id.strip().upper()
    url = f"https://files.rcsb.org/ligands/download/{ccd_id}_ideal.sdf"

    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        raise RuntimeError(
            f"Failed to fetch {ccd_id} from {url} (HTTP {r.status_code}). "
            "If this CCD exists but SDF is unavailable, try the CIF download: "
            f"https://files.rcsb.org/ligands/download/{ccd_id}.cif"
        )

    suppl = Chem.SDMolSupplier()
    suppl.SetData(r.text, sanitize=False, removeHs=False, strictParsing=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise RuntimeError(f"Could not parse SDF for CCD {ccd_id}.")

    mol = mols[0]

    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        print(f"[WARN] Sanitization failed for CCD {ccd_id}: {exc}. Continuing with relaxed molecule.")
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol) 

    return mol

@lru_cache(maxsize=None)
def _fp_from_ccd(ccd_id: str, radius: int, nbits: int, mode="morgan"):
    mol = mol_from_ccd(ccd_id)
    if mode == "morgan":
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nbits)
    elif mode == "atom_pair":
        return AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=nbits)

def tanimoto_ccd(ccd_a: str, ccd_b: str, radius: int = 3, nbits: int = 2048, mode: str = "morgan") -> float:
    """
    Compute Tanimoto similarity between two CCD ligands.
    """
    fp_a = _fp_from_ccd(ccd_a.strip().upper(), radius=radius, nbits=nbits, mode=mode)
    fp_b = _fp_from_ccd(ccd_b.strip().upper(), radius=radius, nbits=nbits, mode=mode)

    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))

def Get_SimilarityTable(CCD_reference_list, CCD_test_list):
    CCD_reference_list = [str(ccd).strip().upper() for ccd in CCD_reference_list if str(ccd).strip()]
    CCD_test_list = [str(ccd).strip().upper() for ccd in CCD_test_list if str(ccd).strip()]

    similarity_table = pd.DataFrame(index=CCD_test_list, columns=CCD_reference_list, dtype=float)

    for test_ccd in CCD_test_list:
        for ref_ccd in CCD_reference_list:
            try:
                similarity_table.loc[test_ccd, ref_ccd] = tanimoto_ccd(test_ccd, ref_ccd, mode="atom_pair")
            except Exception as exc:
                print(f"[WARN] Failed similarity for {test_ccd} vs {ref_ccd}: {exc}")
                similarity_table.loc[test_ccd, ref_ccd] = float("nan")

    similarity_table.index.name = "CCD_test"
    similarity_table.to_csv("CCD_similarity_atompair.csv")
    return similarity_table
    
    

if __name__ == "__main__":
    CCD_all, CCD_test = read_data()
    Get_SimilarityTable(CCD_all, CCD_test)
    
