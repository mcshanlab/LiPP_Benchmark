Steps to run Chai-l on the LiPP dataset


### 1. Install conda environments in the environments/ directory 
- `conda env create -f chai_env.yml`
- `conda env create -f chai_pymol_env.yml`
- `conda env create -f chai_spy_env.yml`

### 2. Run preparation scripts to prepare the inputs for Chai-l
- `conda activate chai`
- `python preprocess.py`

### 3. Run predictions with Chai-l using the environment with Chai-l installation
- run `python predict.py`, (on a cluster with gpu)

### 4. Superimpose protein sturctures
- `conda activate chai_pymol`
- `python run_pymol.py`

### 5. Run Openstructure to calculate protein prediction scores
- `./run_opensturcture.sh`

### 6. Summarize the results
- `conda activate chai_spy`
- `python summarize.py`
