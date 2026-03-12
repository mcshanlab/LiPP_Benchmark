## Instructions to run the scipts to analyze and visualize the evalaution results

### Overview
The `compare.py` script provides visualizations  of different protein-ligand docking and structure prediction methods on the LiPP benchmark dataset. It calculates success rates, generates visualizations, and produces a results report comparing:
- AlphaFold 3 (AF)
- Chai-1 (CHAI)
- AutoDock Vina (VINA)
- RoseTTAFoldAA (RS)
- DiffDock-L (DD)

### Prerequisites

### Required Files
Before running the script, ensure the following input CSV files are present in the same directory:
- `unique_rmsd.csv` - Contains concatenated results of lipid and protein errors for all methods across the benchmark dataset. 

### To reproduce the results from the publication, run the following:

`python compare.py`






