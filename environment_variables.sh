### INSTRUCTIONS
# This script defines the path variables to simplify
# other scripts and allow for quick changes. Run this
# as follows to make the variables available

# source /path/to/environment_variables.sh

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