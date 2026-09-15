#!/usr/bin/env bash

# Author: Jon Slotved
# Description: Functionality script for handling SLURM job submissions in the host element pipeline

#structure of the manifest file: job_counter,array_counter,cov_mode,max_seq_length,query_database_prefix,reference_database_prefix,sample_directory,temporary_directory
write_manifest_file() {
    local processing_files_dir="$1"
    local reference_database_prefix="$2"
    local cov_modes_array="$3"
    local max_seq_lengths_array="$4"
    local manifest_file="$5"

    #task and job counters for SLURM array jobs (1 job=1000 tasks)
    local max_array_size=5
    local array_task_id_counter=0
    local current_chunk=1

    #write manifest file
    touch "$manifest_file"
    for sample_directory in "$processing_files_dir"/*; do
        local sample_name
        sample_name=$(basename "$sample_directory")
        local query_database_prefix="$sample_directory/query_nucl_db/${sample_name}_nucl_db_type_2"
        local results_directory="$sample_directory/results_db"
        for cov_mode in $cov_modes_array; do
            for max_seq_length in $max_seq_lengths_array; do
                local temporary_directory="$sample_directory/tmp/cov_${cov_mode}_max_${max_seq_length}"

                printf "%d,%d,%d,%d,%s,%s,%s,%s\n" \
                       "$current_chunk" "$array_task_id_counter" "$cov_mode" "$max_seq_length" \
                       "$query_database_prefix" "$reference_database_prefix" "$results_directory" "$temporary_directory" \
                       >> "$manifest_file"
                
                ((array_task_id_counter++))
                if [[ $array_task_id_counter -ge $max_array_size ]]; then
                    ((current_chunk++))
                    array_task_id_counter=0
                fi
            done
        done
    done

    if [[ -f "$manifest_file" ]]; then
        write_log "Manifest file written successfully: $manifest_file" "INFO" 
    else
        write_log "Failed to write manifest file: $manifest_file" "ERROR"
    fi
}

start_slurm_runners() {
    local manifest_file="$1"

    local max_simultaneous_jobs="$2"
    local slurm_cpus_per_job="$3"
    local slurm_memory_per_job="$4"
    local slurm_partition="$5"

    local slurm_worker_script="$6"

    write_log "Starting SLURM runners with manifest file: $manifest_file" "INFO"

    numJobs=$(wc -l < "$manifest_file")

    #used Edwards implementation - Jon Slotved
    #find slurm array size based on the number of jobs and maximum simultaneous jobs
    if (( $numJobs % $max_simultaneous_jobs == 0 )); then
        Slurm_chunks=$(($numJobs / $max_simultaneous_jobs))
    else
        Slurm_chunks=$(($numJobs / $max_simultaneous_jobs + 1)) # This is a ceiling int calculation
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

    for ((i=1; i<=Slurm_chunks; i++))
    do
        current_chunk=$i

        array_start=$(cat "$manifest_file" | grep "^$current_chunk" | head -n 1 | awk -F',' '{print $2}')
        array_end=$(cat "$manifest_file" | grep "^$current_chunk" | tail -n 1 | awk -F',' '{print $2}')
        echo "Array start for chunk $current_chunk: $array_start"
        echo "Array end for chunk $current_chunk: $array_end"

        sbatch --array="$array_start-$array_end%$Slurm_CalcRunParallel" \
               --cpus-per-task="$slurm_cpus_per_job" \
               --mem="$slurm_memory_per_job" \
               --partition="$slurm_partition" \
               "$slurm_worker_script" -p "$pipeline_dir" -e "$conda_env_prefix" -m "$manifest_file" -c "$current_chunk"
               
    done
}

