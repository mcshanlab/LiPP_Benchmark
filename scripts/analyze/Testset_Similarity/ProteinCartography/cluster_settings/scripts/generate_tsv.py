import os
import pandas as pd

# Folder containing .pdb files
folder = "../input"  # current directory

# Output TSV file
output_file = "uniprot_features.tsv"

# Define all columns
columns = [
    "protid", "Entry", "Entry Name", "Protein names", "Gene Names (primary)",
    "Annotation", "Organism", "Taxonomic lineage", "Length", "Fragment",
    "Sequence", "Reviewed", "Gene Names", "Protein existence", "Sequence version",
    "RefSeq", "GeneID", "EMBL", "AlphaFoldDB", "PDB", "Pfam", "InterPro", "Lineage"
]

# Collect rows
rows = []

for filename in os.listdir(folder):
    if filename.endswith(".pdb"):
        protid = os.path.splitext(filename)[0]  # remove .pdb extension
        row = {col: "" for col in columns}
        row["protid"] = protid
        rows.append(row)

# Create DataFrame
df = pd.DataFrame(rows, columns=columns)

# Save as TSV
df.to_csv(output_file, sep="\t", index=False)

print(f"TSV file created: {output_file} with {len(rows)} rows.")