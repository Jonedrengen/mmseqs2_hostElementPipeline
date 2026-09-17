#!/usr/bin/env bash
#SBATCH --job-name=mmseqs_worker
#SBATCH --output=mmseqs_worker_%A_%a.out
#SBATCH --error=mmseqs_worker_%A_%a.err

# Author: Jon Slotved
# Description: 
#   Worker script for running SLURM jobs in the host element pipeline
#   will run 1 row of the SLURM metadata file corresponding to the current SLURM array task

#duplication.. But necessary for standalone worker script
write_log() {
    local log_message="${1:-No log message provided}"
    local log_type="${2:-INFO}"
    local log_file="${3:-}"
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
while getopts "p:e:m:c:r:" opt; do
    case $opt in
        p) pipeline_dir="$OPTARG" ;;
        e) conda_env_prefix="$OPTARG" ;;
        m) slurm_meta_info_file="$OPTARG" ;;
        c) current_chunk="$OPTARG" ;;
        r) reference_fasta_file="$OPTARG" ;;
        *) echo "you should not be passing anything here" ;;
    esac
done

source "$pipeline_dir/subscripts/mmseqs_functionality.sh"

#find row, based on current chunk and SLURM array task ID
#datarow: current_chunk,array_task_id,sample_name,trimmed_fasta,query_database_prefix,reference_database_prefix,results_directory,temporary_directory,coverage_modes,max_sequence_lengths,mmseqs_min_seq_id,mmseqs_coverage,mmseqs_sensitivity,mmseqs_max_seqs
data_row_for_worker=$(grep "^${current_chunk},${SLURM_ARRAY_TASK_ID}," "$slurm_meta_info_file")
#read row: chunk, task, sample, trimmed fasta, query db, reference db, results, temporary, coverage modes, max lengths, MMseqs settings
IFS=',' read -r _ _ sample_name trimmed_fasta query_database_prefix reference_database_prefix results_directory temporary_directory coverage_modes max_sequence_lengths mmseqs_min_seq_id mmseqs_coverage mmseqs_sensitivity mmseqs_max_seqs <<< "$data_row_for_worker"
write_log "chunk=$current_chunk task=$SLURM_ARRAY_TASK_ID " "INFO"
write_log "metadata row: $data_row_for_worker" "INFO"

write_log "Starting MMseqs worker for $sample_name" "INFO"
write_nucl_query_db "$conda_env_prefix" \
                    "$trimmed_fasta" \
                    "$query_database_prefix"

#run search script
run_mmseqs_search_and_convert "$conda_env_prefix" \
                              "$query_database_prefix" \
                              "$reference_database_prefix" \
                              "$results_directory" \
                              "$temporary_directory" \
                              "$coverage_modes" \
                              "$max_sequence_lengths" \
                              "$mmseqs_min_seq_id" \
                              "$mmseqs_coverage" \
                              "$mmseqs_sensitivity" \
                              "$mmseqs_max_seqs" \
                              "$results_directory/logs/run.log"

combine_mmseqs_results_per_isolate "$conda_env_prefix" \
                                   "$results_directory" \
                                   "$reference_fasta_file" \
                                   "$pipeline_dir/subscripts/mmseq2_results_replicate_combine.py" \
                                   "$results_directory/logs/run.log"

#moving slurm stuff
mv "$SLURM_SUBMIT_DIR/mmseqs_worker_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" "$(dirname "$results_directory")"
mv "$SLURM_SUBMIT_DIR/mmseqs_worker_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" "$(dirname "$results_directory")"

write_log "Finished MMseqs worker for $sample_name" "INFO"


