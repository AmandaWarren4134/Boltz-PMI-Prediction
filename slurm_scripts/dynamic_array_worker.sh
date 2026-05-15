#!/bin/bash
#SBATCH --job-name=boltz_worker
#SBATCH --mail-user=aw998@byu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
module load boltz2
module load mmseqs2

job_id="$SLURM_ARRAY_TASK_ID"

if [ "$YAML_PER_JOB" -gt 1 ]; then
    echo "YAML_PER_JOB is set to $YAML_PER_JOB. Each job will process multiple YAML files."

    yaml_index_start=$((job_id * YAML_PER_JOB))
    yaml_index_end=$((yaml_index_start + YAML_PER_JOB - 1))

    # Don't go past the total number of YAMLs
    if [ "$yaml_index_end" -ge "$NUM_YAML" ]; then
        yaml_index_end=$((NUM_YAML - 1))
    fi

    echo "----------------------------------------"
    echo "Running task $SLURM_ARRAY_TASK_ID on $(hostname)"
    echo "Processing the following YAML files: $((yaml_index_start + 1)) through $((yaml_index_end + 1))"
    echo "Progress: $((yaml_index_end * job_id + 1)) / $NUM_YAML"
    echo "----------------------------------------"

    # Loop through YAML files
    for ((i = yaml_index_start; i <= yaml_index_end; i++)); do
        YAML_FILE=$(sed -n "$((i + 1))p" "$YAML_LIST")
        echo "Processing YAML: $YAML_FILE"

        boltz predict "$YAML_FILE" \
            --use_msa_server \
            --msa_server_url=http://colabfold02.int.chpc.utah.edu:8088 \
            --num_workers 4 \
            --out_dir "$OUTPUT_DIR"

        echo "Completed YAML: $YAML_FILE"
        echo "----------------------------------------"
    done

    echo ""
    echo "----------------------------------------"
    echo "Task $SLURM_ARRAY_TASK_ID completed YAML indices: $yaml_index_start to $yaml_index_end"
    echo "----------------------------------------"

else
    echo "YAML_PER_JOB is set to 1. Each job will process a single YAML file."
    
    # Pick YAML for this task
    YAML_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$YAML_LIST")
    JOB_NAME=${YAML_FILE:-"boltz_job"}

    echo "----------------------------------------"
    echo "Running task $SLURM_ARRAY_TASK_ID on $(hostname)"
    echo "Processing YAML: $YAML_FILE"
    echo "Progress: $((SLURM_ARRAY_TASK_ID + 1)) / $NUM_YAML"
    echo "----------------------------------------"

    # Run Boltz
    boltz predict "$YAML_FILE" \
        --use_msa_server \
        --msa_server_url=http://colabfold02.int.chpc.utah.edu:8088 \
        --num_workers 4 \
        --out_dir "$OUTPUT_DIR"

    echo ""
    echo "----------------------------------------"
    echo "Task $SLURM_ARRAY_TASK_ID completed YAML: $YAML_FILE"
    echo "----------------------------------------"
fi
