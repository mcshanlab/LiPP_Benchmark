import pandas as pd
import requests


def get_pdb_metadata(pdb_ids):
    url = "https://data.rcsb.org/rest/v1/core/entry/"
    results = {}

    for pdb_id in pdb_ids:
        response = requests.get(url + pdb_id)
        if response.status_code == 200:
            data = response.json()
            classification = data.get("struct_keywords", {}).get("pdbx_keywords")
            method = data.get("exptl", [{}])[0].get("method")
            results[pdb_id] = {
                "classification": classification,
                "method": method
            }
        else:
            results[pdb_id] = {"error": f"Failed to fetch {pdb_id}"}
    return results


def metadata_to_dataframe(metadata):
    df = pd.DataFrame.from_dict(metadata, orient='index')
    df.index.name = 'pdb_id'
    return df


def get_pdb_query(path_csv):
    df = pd.read_csv(path_csv)
    bd_ids = df['BioDolphinID'].tolist()
    pdb_ids = [bd_id.split('-')[0][-4:] for bd_id in bd_ids]
    return list(set(pdb_ids))


def map(csv_path, output_path):
    pdb_list = get_pdb_query(csv_path)
    metadata = get_pdb_metadata(pdb_list)
    df = metadata_to_dataframe(metadata)
    df.to_csv(output_path, index=True)


if __name__ == "__main__":
    map(csv_path="../Evaluation/unique_rmsd.csv", output_path="./pdb_metadata.csv")
