#!/usr/bin/env bash

# Author: Jon Slotved
# Description: Functionality script for handling SLURM job submissions in the host element pipeline

# Structure of the manifest file:
# chunk_id,task_id,sample_name,trimmed_fasta,query_database_prefix,reference_database_prefix,results_directory,temporary_directory,cov_modes,max_seq_lengths
write_manifest_file() {
    local processing_files_dir="$1"
    local trimmed_fasta_dir="$2"
    local sample_id_list_file="$3"
    local reference_database_prefix="$4"
    local cov_modes_array="$5"
    local max_seq_lengths_array="$6"
    local manifest_file="$7"

    #task and job counters for SLURM array jobs (1 job=1000 tasks)
    local max_array_size=100
    local array_task_id_counter=0
    local current_chunk=1

    # Write one row per isolate; the worker owns all searches for that isolate.
    : > "$manifest_file"
    while read -r sample_filename; do
        local sample_name
        sample_name="${sample_filename%.*}"
        local sample_directory="$processing_files_dir/$sample_name"
        local trimmed_fasta="$trimmed_fasta_dir/$sample_filename"
        local query_database_prefix="$sample_directory/query_nucl_db/${sample_name}_nucl_db_type_2"
        local results_directory="$sample_directory/results_db"

        local temporary_directory="$sample_directory/tmp"

        printf "%d,%d,%s,%s,%s,%s,%s,%s,%s,%s\n" \
               "$current_chunk" "$array_task_id_counter" "$sample_name" \
               "$trimmed_fasta" "$query_database_prefix" \
               "$reference_database_prefix" "$results_directory" \
               "$temporary_directory" "$cov_modes_array" "$max_seq_lengths_array" \
               >> "$manifest_file"

        ((array_task_id_counter++))
        if [[ $array_task_id_counter -ge $max_array_size ]]; then
            ((current_chunk++))
            array_task_id_counter=0
        fi
    done < "$sample_id_list_file"

    if [[ -f "$manifest_file" ]]; then
        write_log "Manifest file written successfully: $manifest_file" "INFO" 
    else
        write_log "Failed to write manifest file: $manifest_file" "ERROR"
    fi
}

start_slurm_runners() {
    local manifest_file="$1"

    local max_jobs_per_array="$2"
    local slurm_cpus_per_job="$3"
    local slurm_memory_per_job="$4"
    local slurm_partition="$5"

    local slurm_worker_script="$6"

    write_log "Starting SLURM runners with manifest file: $manifest_file" "INFO"

    numJobs=$(wc -l < "$manifest_file")

    #used Edwards implementation - Jon Slotved
    #find slurm array size based on the number of jobs and maximum simultaneous jobs
    if (( $numJobs % $max_jobs_per_array == 0 )); then
        Slurm_chunks=$(($numJobs / $max_jobs_per_array))
    else
        Slurm_chunks=$(($numJobs / $max_jobs_per_array + 1)) # This is a ceiling int calculation
    fi
    # Most used nodes 12, minimum is 1, it will only use 1 if its submits more than 6 batches
    # Please adjust these numbers accordingly to your specifications or HPC needs
    # Example: if samplelist contains 1000 files, it will submit 1 SlurmArray job that will run use 12 compute nodes at a time.
    # Example: if samplelist contains 10000 files, it will submit 10 SlurmArray jobs that will run only 1 compute node per SlurmArray job at a time.
    if [ $Slurm_chunks == 1 ]
    then
    Slurm_CalcRunParallel=12

    elif [ $Slurm_chunks == 2 ]
    then
    Slurm_CalcRunParallel=6

    elif [ $Slurm_chunks == 3 ]
    then
    Slurm_CalcRunParallel=4

    elif [ $Slurm_chunks == 4 ]
    then
    Slurm_CalcRunParallel=3

    elif [ $Slurm_chunks == 5 ]
    then
    Slurm_CalcRunParallel=2

    elif [ $Slurm_chunks == 6 ]
    then
    Slurm_CalcRunParallel=2

    else
    Slurm_CalcRunParallel=1
    fi

    write_log "running $(( Slurm_CalcRunParallel * Slurm_chunks )) jobs with $Slurm_CalcRunParallel jobs across $Slurm_chunks chunks in parallel" "INFO"

    for ((current_chunk=1; current_chunk<=Slurm_chunks; current_chunk++))
    do
        array_start=$(cat "$manifest_file" | grep "^$current_chunk" | head -n 1 | awk -F',' '{print $2}')
        array_end=$(cat "$manifest_file" | grep "^$current_chunk" | tail -n 1 | awk -F',' '{print $2}')
        echo "Array start for chunk $current_chunk: $array_start"
        echo "Array end for chunk $current_chunk: $array_end"

        sbatch --array="$array_start-$array_end%$Slurm_CalcRunParallel" \
               --cpus-per-task="$slurm_cpus_per_job" \
               --mem="$slurm_memory_per_job" \
               --partition="$slurm_partition" \
               --job-name="mmseqs_worker_${current_chunk}_${array_start}-${array_end}" \
               "$slurm_worker_script" -p "$pipeline_dir" -e "$conda_env_prefix" -m "$manifest_file" -c "$current_chunk"
               
    done
    
}

