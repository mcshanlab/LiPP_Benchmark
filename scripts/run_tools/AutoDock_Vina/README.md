Steps to run AutoDock Vina

### 1. Install the conda environment from yml files in environments

- `conda env create -f ligand_env.yml`
- `conda env create -f protein_env.yml`  
- `conda env create -f vina_new_env.yml`



### 2. Run preparation scripts to prepare the inputs for AutoDock Vina

- `cd Step1:File_Conversion`

##### Process protein files:

- `conda activate openbabel_env`
- `python protein_conversionV4.py`

##### Process lipid files:

- `conda activate new_openbabel_env`
- `python ligand_conversion_script.py`

##### Define Grid Box:

- `cd Step2:GridBox_Creation`
- `conda activate vina_env_new`
- `python pocket_definitionV5.py`



### 3. Submit shell script to run AutoDock Vina

- `conda activate vina_env_new`
- `python run_vina.py`

