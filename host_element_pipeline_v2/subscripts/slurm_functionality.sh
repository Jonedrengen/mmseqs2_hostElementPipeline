#!/usr/bin/env bash

# Author: Jon Slotved
# Description: Functionality script for handling SLURM job submissions in the host element pipeline

# this is the same as Edwards implementation with chunk in col 1 and task in col 2
# Structure of the SLURM metadata file:
# chunk_id,task_id,sample_name,trimmed_fasta,query_database_prefix,reference_database_prefix,results_directory,temporary_directory,cov_modes,max_seq_lengths
write_slurm_meta_info_file() {
    local processing_files_dir="$1"
    local trimmed_fasta_dir="$2"
    local sample_id_list_file="$3"
    local reference_database_prefix="$4"
    local cov_modes_array="$5"
    local max_seq_lengths_array="$6"
    local slurm_meta_info_file_name="$7"
    local log_file="${8:-}"

    #task and job counters for SLURM array jobs (1 job=1000 tasks)
    local max_array_size=100
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

        printf "%d,%d,%s,%s,%s,%s,%s,%s,%s,%s\n" \
               "$current_chunk" "$array_task_id_counter" "$sample_name" "$trimmed_fasta" "$query_database_prefix" "$reference_database_prefix" "$results_directory" "$temporary_directory" "$cov_modes_array" "$max_seq_lengths_array" \
               >> "$slurm_meta_info_file_name"

        ((array_task_id_counter++))
        if [[ $array_task_id_counter -ge $max_array_size ]]; then
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
    local slurm_cpus_per_job="$5"
    local slurm_memory_per_job="$6"
    local slurm_partition="$7"

    local slurm_worker_script="$8"
    local log_file="${9:-}"

    write_log "Starting SLURM runners with SLURM metadata file: $slurm_meta_info_file_name" "INFO" "$log_file"

    local num_jobs
    local slurm_chunks
    local slurm_calc_run_parallel
    num_jobs=$(wc -l < "$slurm_meta_info_file_name")

    #used Edwards implementation - Jon Slotved
    #find slurm array size based on the number of jobs and maximum simultaneous jobs
    if (( num_jobs % max_jobs_per_array == 0 )); then
        slurm_chunks=$((num_jobs / max_jobs_per_array))
    else
        slurm_chunks=$((num_jobs / max_jobs_per_array + 1)) # This is a ceiling int calculation
    fi
    # Most used nodes 12, minimum is 1, it will only use 1 if its submits more than 6 batches
    # Please adjust these numbers accordingly to your specifications or HPC needs
    # Example: if samplelist contains 1000 files, it will submit 1 SlurmArray job that will run use 12 compute nodes at a time.
    # Example: if samplelist contains 10000 files, it will submit 10 SlurmArray jobs that will run only 1 compute node per SlurmArray job at a time.
    if [[ $slurm_chunks == 1 ]]
    then
    slurm_calc_run_parallel=12

    elif [[ $slurm_chunks == 2 ]]
    then
    slurm_calc_run_parallel=6

    elif [[ $slurm_chunks == 3 ]]
    then
    slurm_calc_run_parallel=4

    elif [[ $slurm_chunks == 4 ]]
    then
    slurm_calc_run_parallel=3

    elif [[ $slurm_chunks == 5 ]]
    then
    slurm_calc_run_parallel=2

    elif [[ $slurm_chunks == 6 ]]
    then
    slurm_calc_run_parallel=2

    else
    slurm_calc_run_parallel=1
    fi

    write_log "running $(( slurm_calc_run_parallel * slurm_chunks )) jobs: $slurm_calc_run_parallel jobs across $slurm_chunks chunks in parallel" "INFO" "$log_file"

    local current_chunk
    local array_start
    local array_end
    for ((current_chunk=1; current_chunk<=slurm_chunks; current_chunk++))
    do
        array_start=$(grep "^${current_chunk}," "$slurm_meta_info_file_name" | head -n 1 | cut -d',' -f2)
        array_end=$(grep "^${current_chunk}," "$slurm_meta_info_file_name" | tail -n 1 | cut -d',' -f2)
        echo "Array start for chunk $current_chunk: $array_start"
        echo "Array end for chunk $current_chunk: $array_end"

        sbatch --array="$array_start-$array_end%$slurm_calc_run_parallel" \
               --cpus-per-task="$slurm_cpus_per_job" \
               --mem="$slurm_memory_per_job" \
               --partition="$slurm_partition" \
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

    write_log "Starting SLURM compiler" "INFO" "$log_file"

    sbatch --dependency=singleton \
           --cpus-per-task=1 \
           --mem="$slurm_memory_per_job" \
           --partition="$slurm_partition" \
           --job-name="mmseqs_worker_gogogogo" \
           "$slurm_compiler_script" -p "$pipeline_dir" \
                                    -e "$conda_env_prefix" \
                                    -o "$output_dir" \
                                    -h "$host_file" \
                                    -r "$reference_fasta_file" \
                                    -l "$log_file"
}

