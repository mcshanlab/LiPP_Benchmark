## Scripts to run PoseBusters Test Suite

### Step1: Prepare csv files as inputs for PoseBusters

#### Example to analyze outputs of AF3:

```bash
python prepare_input_csv.py
  --output_name AF3
  --lipid_path "../BondOrder/fixed_sdfs/AF/{BDID}_alignLipid.sdf"
  --protein_path "../../../run_tools/AF3/aligned_af_output/{BDID}_alignProtein.pdb"
```

- lipid_path: The path to the lipid sdf files with assigned bond orders
- protein_path: The path to the protein pdb files

### Step2: Run PoseBusters

```bash
python run_buster.py --csv-file path/to/input.csv
```
This will generate an output csv file  in the current directory containing the validation results.

### Step3: Plot the results

```bash
python plot_buster.py
```
This will plot the test results. (Change the prefix in the script if you named the result csv files in step1 and step2 differently)