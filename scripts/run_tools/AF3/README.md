Steps to run AF3 on the LiPP dataset


### 1. Install conda environments in the environments/ directory 
- `conda env create -f af3_env.yml`
- `conda env create -f af3_pymol_env.yml`
- `conda env create -f af3_spy_env.yml`


### 2. Run preparation scripts to prepare the inputs for AF3
- `conda activate af3`
- `python preprocess.py`

### 3. Submit shell script to run AF3
- Run `./run.sh` with a customized slurm script or equivalent on your cluster

### 4. Superimpose protein sturctures
- `conda activate af3_pymol`
- `python run_pymol.py`

### 5. Run Openstructure to calculate protein prediction scores
- `./run_opensturcture.sh`

### 6. Summarize the results
- `conda activate af3_spy`
- `python summarize.py`
