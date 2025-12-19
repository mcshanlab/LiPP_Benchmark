Filter protein-lipid complexes from BioDolphin dataset

## Environments

install the environments with `environment.yml` \

`conda activate filter_env` 

## Steps to run the code

### 1. Place the BioDolphin dataset files in the `filtering` directory (ex: `BioDolphin_vr1.1.csv`) 

### 2. Run all filters (except distance filters) to the dataset: 

- `python main.py` 

### 3. Run the last filter for distance calculation (longer runtime, best to submit on a cluster node): 

- `python main.py --last_filter`

### 4. Combine the files from step3 and get final statistics: 

- `python main.py --combine` 

