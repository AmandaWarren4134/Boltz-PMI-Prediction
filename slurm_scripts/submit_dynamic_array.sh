#!/bin/bash

module load boltz2

JOB_NAME="${1}"
YAML_PER_JOB="${2:-1}"

# Ensure environment variables are loaded
# Modify the path below to target 
# the environment file
ENV_FILE="/uufs/chpc.utah.edu/common/home/u6073680/Capstone/BIO-465-Capstone-Projectcapstone_path_env.sh"
if [ -f "$ENV_FILE" ]; then
    . "$ENV_FILE"
fi

# Check for environment variables
for var in PROJECT_ROOT OUTPUT_DIR YAML_DIR LOG_DIR YAML_LIST_DIR SCRIPT_DIR; do
    if [ -z "${!var}" ]; then
        echo "Error: $var is not set. Environment variables are not loaded." >&2
        exit 1
    fi
done

# Paths
OUTPUT_LOC="$OUTPUT_DIR/$JOB_NAME"
YAML_LOC="$YAML_DIR/$JOB_NAME"
LOG_LOC="$LOG_DIR/boltz_array_$JOB_NAME"

# Make directories
mkdir -p "$OUTPUT_LOC" "$LOG_LOC" "$YAML_LIST_DIR"

# Generate YAML list file
# YAML_LIST="$YAML_LIST_DIR/yaml_list_$JOB_NAME.txt"
# ls "$YAML_LOC"/*.yaml > "$YAML_LIST"
# NUM_YAML=$(wc -l < "$YAML_LIST")

# Produce YAML list file
python3 "$SCRIPT_DIR/prep_yaml_lists.py" --output_folder "$JOB_NAME" --create_list
YAML_LIST="$YAML_LIST_DIR/yaml_list_$JOB_NAME.txt"
NUM_YAML=$(wc -l < "$YAML_LIST")

# Ensure there are .YAML files in YAML_LOC
if [ "$NUM_YAML" -eq 0 ]; then
    echo "Error: No YAML files found in $YAML_LOC"
    exit 1
fi

# Use batching if there are more than 1000 YAML files
if [ "$NUM_YAML" -gt 1000 ] && [ "$YAML_PER_JOB" -eq 1 ]; then
    echo "Warning: Found $NUM_YAML YAML files in $JOB_NAME. Batching will be used."

    # Use smallest YAML_PER_JOB which 
    # keeps the array size at/below 1000
    YAML_PER_JOB=$(( (NUM_YAML + 999) / 1000 ))
    echo "Setting YAML_PER_JOB to $YAML_PER_JOB."
fi

echo "Found $NUM_YAML YAML files in $YAML_LOC"
echo "YAML list saved to $YAML_LIST"

# Calculate array indices
# Number of jobs is ceiling of (NUM_YAML / YAML_PER_JOB)
if [ "$YAML_PER_JOB" -eq 1 ]; then
    MAX_INDEX=$((NUM_YAML - 1))
else
    TASKS_NEEDED=$(( (NUM_YAML + YAML_PER_JOB - 1) / YAML_PER_JOB ))
    MAX_INDEX=$((TASKS_NEEDED - 1))
fi

sbatch \
    --job-name=$JOB_NAME \
    --time=00:30:00 \
    --cpus-per-task=1 \
    --gres=gpu:h200:1 \
    --mem=64G \
    --cluster=granite \
    --partition=rai-gpu-grn \
    --qos=rai-gpu-grn \
    --account=rai \
    --mail-user=aw998@byu.edu \
    --mail-type=BEGIN,END,FAIL \
    --array=0-$MAX_INDEX%10 \
    --export=ALL,YAML_LIST="$YAML_LIST",OUTPUT_DIR="$OUTPUT_LOC",YAML_PER_JOB="$YAML_PER_JOB",NUM_YAML="$NUM_YAML" \
    --output="$LOG_LOC/task_%A_%a.log" \
    --error="$LOG_LOC/task_%A_%a.log" \
    "$SLURM_SCRIPT_DIR/dynamic_array_worker.sh"

echo "Submitted a total of $((MAX_INDEX + 1)) jobs with $NUM_YAML YAML files, $YAML_PER_JOB per job, for $JOB_NAME."
