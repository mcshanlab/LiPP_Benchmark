Steps to run DiffDock-L on the DolphinBusters dataset

1. Install conda environments from environments/
- `conda env create -f diff_env.yml`
- `conda env create -f diff_pymol.yml`
- `conda env create -f diff_spy.yml`

2. Run preparation scripts to prepare the inputs for DiffDock-L
- `conda activate diff_env`
- `python preprocess.py`, which will output the file `all_input.csv`

3.  Run predictions with DiffDock-L's container with `all_input.csv` as inputs
- `bash run.sh` submit it to a slurm script or equivalent on your cluster with gpus

5. Evaluate the results
- `conda activate diff_pymol`
- `python run_pymol.py`
- `conda activate diff_spy`
- `python summarize.py`