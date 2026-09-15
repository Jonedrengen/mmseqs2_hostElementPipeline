#!/usr/bin/env bash
#SBATCH --job-name=mmseqs_worker
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=mmseqs_worker_%A_%a.out
#SBATCH --error=mmseqs_worker_%A_%a.err
#SBATCH --time=04:00:00

# Author: Jon Slotved
# Description: 
#   Worker script for running SLURM jobs in the host element pipeline
#   will run 1 row of the manifest file corresponding to the current SLURM array task

#duplication.. But necessary for standalone worker script
write_log() {
    local log_message="${1:-No log message provided}"
    local log_type="${2:-INFO}"
    local log_file="$3"
    local time=""
    time=$(date +"%Y-%m-%d %H:%M:%S")
    #stdout
    echo "[$log_type] [$time] $log_message"
    #optional to log file
    if [[ ! -z "$log_file" ]]; then
        echo "[$log_type] [$time] $log_message" >> "$log_file"
    fi
}

#inputs for the SLURM worker script
while getopts "p:e:m:c:" opt; do
    case $opt in
        p) pipeline_dir="$OPTARG" ;;
        e) conda_env_prefix="$OPTARG" ;;
        m) manifest_file="$OPTARG" ;;
        c) current_chunk="$OPTARG" ;;
        *) echo "you should not be passing anything here" ;;
    esac
done
source "$pipeline_dir/subscripts/mmseqs_functionality.sh"

data_row_for_worker=$(awk -F',' -v chunk="$current_chunk" -v task_id="$SLURM_ARRAY_TASK_ID" '$1 == chunk && $2 == task_id {print}' "$manifest_file")
IFS=',' read -r _ _ coverage_mode max_sequence_length query_database_prefix reference_database_prefix results_directory temporary_directory <<< "$data_row_for_worker"



#run search script
run_mmseqs_search_and_convert "$conda_env_prefix" \
                              "$query_database_prefix" \
                              "$reference_database_prefix" \
                              "$results_directory" \
                              "$temporary_directory" \
                              "$coverage_mode" \
                              "$max_sequence_length"


