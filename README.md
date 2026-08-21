# How to do Large-Scale Predictions of Protein-Metabolite Interactions Using our Boltz-2 Pipeline

This guide explains the main steps needed to run the project.

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

## Notes
- Full Boltz job submission requires the HPC/Slurm environment.
- Local testing is mainly for validating paths, generating YAMLs, and creating figures.