# How to do Large-Scale Predictions of Protein-Metabolite Interactions Using our Boltz-2 Pipeline

This guide explains the main steps needed to run the project.

## Step 0: Overall Setup and Storage
- This repository is meant to be cloned into your home directory
- You supply your own inputs in csv format and supply their relative or absolute filepaths as arguments to run_pipeline.sh (it doesn't really matter where they are)
- Boltz outputs are directed to /scratch/rai/vast1/stewartp/boltz_results/<boltz_results_your_job_name>
- Running parse_outputs_to_csv.py puts the outputs wherever you have configured in environment_variables.sh

In summary, everything should be configurable to your environment, locally or on the HPC.

---

## Step 1: Set up the Python environment and install `requirements.txt`

This project was tested using a Conda environment with Python 3.11.

Create the environment:

```bash
conda create -n boltz-pmi python=3.11
```

Activate the environment:

```bash
conda activate boltz-pmi
```

Install the required packages:

```bash
pip install -r requirements.txt
```

Your `requirements.txt` includes the packages needed by the scripts, such as:

```txt
pandas
numpy
matplotlib
seaborn
scipy
gseapy
mygene
pyrolite
```

---

## Step 2: Set up the path variables in `environment_variables.sh`

Before running the scripts, update the environment file so the paths match your machine.

Open:

```bash
environment_variables.sh
```

Make sure the variables point to the correct folders.

Example local setup:

```bash
### PATHS
# Modify this to point to the local project directory
export PROJECT_ROOT="/path/to/project/directory" # Where you cloned the github repo

# SCRATCH_ROOT is the location for outputs. If not
# using the University of Utah HPC, this can be set
# to the PROJECT_ROOT or another local directory
export SCRATCH_ROOT="/path/to/scratch/space"

export SLURM_SCRIPT_DIR="$PROJECT_ROOT/slurm_scripts"
export SCRIPT_DIR="$PROJECT_ROOT/scripts"
export LOG_DIR="$PROJECT_ROOT/logs"
export YAML_LIST_DIR="$PROJECT_ROOT/slurm_scripts/yaml_lists"
export PROCESSED_DIR="$PROJECT_ROOT/processed_outputs"
export PLOT_DIR="$PROJECT_ROOT/plots"


export OUTPUT_DIR="$SCRATCH_ROOT/boltz_results"
export YAML_DIR="$SCRATCH_ROOT/yaml"
export MSA_BASE_DIR="$SCRATCH_ROOT/msa_per_protein/proteome" # Use this if you have downloaded local MSA files
```

Once you have established your environment variables, make them available by running:
```bash
source environment_variables.sh
```

### Important variables

- `PROJECT_ROOT` = root of the repository 
- `SCRATCH_ROOT` = scratch space to store large boltz outputs
- `INPUT_DIR` = folder containing input CSV files
- `PROCESSED_DIR` = processed Boltz output CSV files
- `PLOT_DIR` = where figures will be saved

---

## Step 3: Customize Slurm Submission

Open: 
```bash
dynamic_array_worker.sh
```

Replace "# YOUR EMAIL HERE" with your email to receive job notifications

Open: 
```bash
submit_dynamic_array.sh
```

Replace lines 76-93 with job submission details specific to your account. Replace "# YOUR EMAIL HERE". The default job submission is set up for a user with access to the RAI cluster.

## Step 4: Run the pipeline

# Run the pipeline:
```bash
bash run_pipeline.sh --protein_csv <protein file> --ligand_csv <ligand file> [--proteome] [--dry-run]
```

- `--proteome`: download proteome if no protein file is given
- `--dry-run`: run without submitting jobs

The protein file may be in CSV or FASTA format. The ligand file must be in CSV format. 

## Step 5: Compile Boltz results into CSV

The `parse_output_to_csv.py` script has two modes: process a specific output folder or process all output folders

### Process a specific output folder
Run this script with Python:

```bash
python parse_output_to_csv.py --output_folder job_name
```

Multiple output folders can be supplied:

```bash
python parse_output_to_csv.py --output_folder job_1 job_2 job_3
```

### Process all output folders

```bash
python parse_output_to_csv.py --process_all
```

This processes every subdirectory found directly under `OUTPUT_DIR`.

### Specify the output CSV name
When processing a single folder, use --job_name to specify the name used in the output CSV:

```bash
python parse_output_to_csv.py \
    --output_folder job_name \
    --job_name my_analysis
```

This creates: `PROCESSED_DIR/confidence_predictions_my_analysis.csv`

If multiple folders are processed, the script ignores --job_name and uses each folder's name instead.

### Enable verbose output
Use --verbose to report missing confidence or affinity files:

```bash
python parse_output_to_csv.py \
    --output_folder job_name \
    --verbose
```

### Output
For each successfully processed output folder, the script creates:
`confidence_predictions_<job_name>.csv`
The CSV contains the original confidence and affinity fields along with:

- protein
- ligand
- micromolar_affinity_pred_value
- molar_affinity_pred_value
- neg_log_molar_affinity_pred_value

The derived affinity values are calculated from affinity_pred_value, whcih is assumed to be expressed as log10(micromolar affinity):

```
micromolar_affinity = 10 ** affinity_pred_value
molar_affinity      = micromolar_affinity / 1e6
negative_log_molar  = -log10(molar_affinity)
```

### Missing Files
A prediction is included in the CSV only when both required files are present:
```
confidence_<stem>_model_0.json
affinity_<stem>.json
```
Missing files are skipped. With --verbose, the script reports the number of missing confidence and affinity files after processing each output folder.

## Notes
- Full Boltz job submission requires the HPC/Slurm environment.
- Local testing is mainly for validating paths, generating YAMLs, and creating figures.
- Folder names are expected to follow Boltz naming convention `boltz_results_{protein}__{ligand}`.
- Folders that do not contain both required JSON files are skipped.