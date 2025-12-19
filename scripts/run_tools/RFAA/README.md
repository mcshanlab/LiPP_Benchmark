Steps to run RFAA

### 1. Install the conda environment from yml files in environments

- `conda env create -f pdb2fasta.yml`
- `conda env create -f rfaa.yml`  
- `conda env create -f eval.yml`


### 2. Run Step1 to prepare inputs for RFAA
- `cd Step1_File_conversion`
- `conda activate openbabel_env`
- `python pdb2fasta.py`
- This will create ../input with the processed inputs

### 3. Run Step2 to run RFAA
- `cd ../Step2_Run_RFAA`
- Edit base.yaml to the correct path of your RFAA installation
- `conda activate RFAA`
- `python run_rosetta.py`

### 4. Run Step3 to evalaute the results
- `cd ../Step3_Evaluation`
- `conda activate rfaa_eval`
- `python rmsd_rosetta.py`


### 5. Run Step4 to calculate protein errors for the results
- `cd ../Step4_Summarize`
- `./run_openstructure.sh`


