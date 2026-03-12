import pandas as pd
import re
import ast

# File paths
uniprot_file = "uniprot_features.tsv"
biodolphin_file = "../BioDolphin_vr1.1.csv"
output_file = "uniprot_features_updated.tsv"

# Load data
df_uniprot = pd.read_csv(uniprot_file, sep="\t")
df_bd = pd.read_csv(biodolphin_file)

# Fill defaults for missing columns
default_values = {
    "Entry": "",
    "Entry Name": "",
    "Protein names": "",
    "Gene Names (primary)": "",
    "Annotation": "",
    "Organism": "",
    "Taxonomic lineage": "[]",
    "Length": "",
    "Fragment": "",
    "Sequence": "",
    "Reviewed": "reviewed",
    "Gene Names": "",
    "Protein existence": "Evidence at protein level",
    "Sequence version": "1",
    "RefSeq": "",
    "GeneID": "",
    "EMBL": "",
    "AlphaFoldDB": "",
    "PDB": "",
    "Pfam": "",
    "InterPro": "",
    "Lineage": "[]"
}

# Map BioDolphin columns to Uniprot columns
mapping = {
    "protein_UniProt_ID": "Entry",
    "protein_Name": "Entry Name",
    "protein_Synonyms": "Protein names",
    "protein_Gene": "Gene Names (primary)",
    "protein_Organism": "Organism",
    "protein_Pfam": "Pfam",
    "protein_InterPro": "InterPro",
    "protein_Sequence": "Sequence",
}

# Prepare BioDolphin dictionary for fast lookup
bd_dict = df_bd.set_index("BioDolphinID").to_dict(orient="index")

# Function to extract PF IDs from string or list
def extract_pf_ids(pfam_value):
    try:
        # Convert string representation of list to Python list
        if isinstance(pfam_value, str):
            pfam_list = ast.literal_eval(pfam_value)
        elif isinstance(pfam_value, list):
            pfam_list = pfam_value
        else:
            return ""
        ids = []
        for item in pfam_list:
            match = re.search(r"\[(PF\d+)\]", str(item))
            if match:
                ids.append(match.group(1))
        return ";".join(ids) + ";" if ids else ""
    except:
        return ""

# Update uniprot dataframe
for i, row in df_uniprot.iterrows():
    protid = row["protid"]  # keep protid as-is
    # Fill defaults first
    for col, val in default_values.items():
        df_uniprot.at[i, col] = val
    if protid in bd_dict:
        bd_row = bd_dict[protid]
        for bd_col, uni_col in mapping.items():
            if bd_col in bd_row:
                value = bd_row[bd_col]
                if uni_col == "Pfam":
                    value = extract_pf_ids(value)
                elif isinstance(value, list):
                    value = ";".join([str(v) for v in value]) + ";"
                df_uniprot.at[i, uni_col] = value


# Add dataset labels
# Load the test set BioDolphinIDs
test_ids = set(
    pd.read_csv("../test_set.csv")["BioDolphinID"]
)
print(f'test_ids: {test_ids}')

df_uniprot["Dataset"] = (
    df_uniprot["protid"]
    .isin(test_ids)
    .map({True: "Test", False: "General"})
)



# Save updated dataframe
df_uniprot.to_csv(output_file, sep="\t", index=False)

print(f"Updated uniprot_features.tsv saved to {output_file}")
