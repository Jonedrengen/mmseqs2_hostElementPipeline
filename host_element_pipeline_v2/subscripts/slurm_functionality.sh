#!/usr/bin/env bash

# Author: Jon Slotved
# Description: Functionality script for handling SLURM job submissions in the host element pipeline

# this is the same as Edwards implementation with chunk in col 1 and task in col 2
# Structure of the SLURM metadata file:
# chunk_id,task_id,sample_name,trimmed_fasta,query_database_prefix,reference_database_prefix,results_directory,temporary_directory,cov_modes,max_seq_lengths,mmseqs_min_seq_id,mmseqs_coverage,mmseqs_sensitivity,mmseqs_max_seqs
write_slurm_meta_info_file() {
    # Arguments:
    #   1: processing_files_dir       <path> e.g. /analysis/processing_files
    #   2: trimmed_fasta_dir           <path> e.g. /analysis/500_bpTrimmed_fastas
    #   3: sample_id_list_file         <path> e.g. /analysis/sample_ID_list.txt
    #   4: reference_database_prefix  <path> e.g. /analysis/tmp/reference_db/reference_nucl_db_type_2
    #   5: cov_modes_array             <space-separated INTs>
    #   6: max_seq_lengths_array       <space-separated INTs>
    #   7: mmseqs_min_seq_id           <number>
    #   8: mmseqs_coverage             <number>
    #   9: mmseqs_sensitivity          <number>
    #   10: mmseqs_max_seqs            <integer>
    #   11: max_jobs_per_array         <integer>
    #   12: slurm_meta_info_file_name  <path> e.g. /analysis/slurm_meta_info.csv
    #   13: log_file                   <path> e.g. /analysis/logs/run.log (optional; last argument)

    local processing_files_dir="$1"
    local trimmed_fasta_dir="$2"
    local sample_id_list_file="$3"
    local reference_database_prefix="$4"
    local cov_modes_array="$5"
    local max_seq_lengths_array="$6"
    local mmseqs_min_seq_id="$7"
    local mmseqs_coverage="$8"
    local mmseqs_sensitivity="$9"
    local mmseqs_max_seqs="${10}"
    local max_jobs_per_array="${11}"
    local slurm_meta_info_file_name="${12}"
    local log_file="${13:-}"

    local array_task_id_counter=0
    local current_chunk=1

    # Write one row per isolate; the worker owns all searches for that isolate.
    : > "$slurm_meta_info_file_name"
    while read -r sample_filename; do
        local sample_name
        sample_name="${sample_filename%.*}"
        local sample_directory="$processing_files_dir/$sample_name"
        local trimmed_fasta="$trimmed_fasta_dir/$sample_filename"
        local query_database_prefix="$sample_directory/query_nucl_db/${sample_name}_nucl_db_type_2"
        local results_directory="$sample_directory/results_db"

        local temporary_directory="$sample_directory/tmp"

         printf "%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
             "$current_chunk" "$array_task_id_counter" "$sample_name" "$trimmed_fasta" "$query_database_prefix" "$reference_database_prefix" "$results_directory" "$temporary_directory" "$cov_modes_array" "$max_seq_lengths_array" "$mmseqs_min_seq_id" "$mmseqs_coverage" "$mmseqs_sensitivity" "$mmseqs_max_seqs" \
               >> "$slurm_meta_info_file_name"

        ((array_task_id_counter++))
        if [[ $array_task_id_counter -ge $max_jobs_per_array ]]; then
            ((current_chunk++))
            array_task_id_counter=0
        fi
    done < "$sample_id_list_file"

    if [[ -f "$slurm_meta_info_file_name" ]]; then
        write_log "SLURM metadata file written successfully: $slurm_meta_info_file_name" "INFO" "$log_file"
    else
        write_log "Failed to write SLURM metadata file: $slurm_meta_info_file_name" "ERROR" "$log_file"
    fi
}

start_slurm_runners() {
    local conda_env_prefix="$1"
    local pipeline_dir="$2"
    local slurm_meta_info_file_name="$3"

    local max_jobs_per_array="$4"
    local max_parallel_jobs_per_array="$5"
    local slurm_cpus_per_job="$6"
    local slurm_memory_per_job="$7"
    local slurm_partition="$8"

    local slurm_worker_script="$9"
    local log_file="${10:-}"

    write_log "Starting SLURM runners with SLURM metadata file: $slurm_meta_info_file_name" "INFO" "$log_file"

    local num_jobs
    local slurm_chunks
    num_jobs=$(wc -l < "$slurm_meta_info_file_name")

    #used Edwards implementation - Jon Slotved
    #find slurm array size based on the number of jobs and maximum parallel jobs
    if (( num_jobs % max_jobs_per_array == 0 )); then
        slurm_chunks=$((num_jobs / max_jobs_per_array))
    else
        slurm_chunks=$((num_jobs / max_jobs_per_array + 1)) # This is a ceiling int calculation
    fi
    write_log "running $num_jobs jobs across $slurm_chunks arrays with a maximum of $max_parallel_jobs_per_array parallel jobs per array" "INFO" "$log_file"

    local current_chunk
    local array_start
    local array_end
    for ((current_chunk=1; current_chunk<=slurm_chunks; current_chunk++))
    do
        array_start=$(grep "^${current_chunk}," "$slurm_meta_info_file_name" | head -n 1 | cut -d',' -f2)
        array_end=$(grep "^${current_chunk}," "$slurm_meta_info_file_name" | tail -n 1 | cut -d',' -f2)
        echo "Array start for chunk $current_chunk: $array_start"
        echo "Array end for chunk $current_chunk: $array_end"

        sbatch --array="$array_start-$array_end%$max_parallel_jobs_per_array" \
               --cpus-per-task="$slurm_cpus_per_job" \
               --mem="$slurm_memory_per_job" \
               --partition="$slurm_partition" \
               --time=04:00:00 \
               --job-name="mmseqs_worker_gogogogo" \
               "$slurm_worker_script" -p "$pipeline_dir" -e "$conda_env_prefix" -m "$slurm_meta_info_file_name" -c "$current_chunk"
    done

}

start_slurm_compiler() {
    local conda_env_prefix="$1"
    local pipeline_dir="$2"
    local output_dir="$3"
    local host_file="$4"
    local reference_fasta_file="$5"
    local slurm_memory_per_job="$6"
    local slurm_partition="$7"
    local slurm_compiler_script="$8"
    local log_file="${9:-}"

    write_log "SLURM compiler" "INFO" "$log_file"

    sbatch --dependency=singleton \
           --cpus-per-task=1 \
           --mem="$slurm_memory_per_job" \
           --partition="$slurm_partition" \
           --time=04:00:00 \
           --job-name="mmseqs_worker_gogogogo" \
           "$slurm_compiler_script" -p "$pipeline_dir" \
                                    -e "$conda_env_prefix" \
                                    -o "$output_dir" \
                                    -h "$host_file" \
                                    -r "$reference_fasta_file" \
                                    -l "$log_file"
}

