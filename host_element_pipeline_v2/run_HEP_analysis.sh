#!/usr/bin/env bash
#SBATCH --job-name=HEP_analysis
#SBATCH --output=HEP_analysis.out
#SBATCH --error=HEP_analysis.err
#SBATCH --partition=project
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

#notes for slurm:
# use 16GB mem per cpu. too little will block the db being loaded properly

#author: Jon Slotved (JOSS@dksund.dk)

# The flow is sequential
# helper functions:
#   write_log is heavily used throughout the script
# main functions:
#   runs sequencially
#   run_mmseqs_search_and_convert and parallel_mmseqs_searches are 1 cohesive unit

#######################################
############# helper functions ########
#######################################

help() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -c <config_file> [-h]"
    echo "  -i <input_dir>       Input directory containing FASTA files"
    echo "  -o <output_dir>      Output directory for compiled results and logs"
    echo "  -c <config_file>     Configuration file specifying source directory and conda environment prefix"
    echo "  -h                   Display this help message"
}

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

write_version_info() {
    local conda_env_prefix="$1"
    local log_file="$2"
    
    write_log "mmseqs version: $(conda run -p "$conda_env_prefix" mmseqs version)" "INFO" "$log_file"
}

#######################################
############# main functions ##########
#######################################

validate_input() {
    local input_dir="$1"
    local output_dir="$2"
    local config_file="$3"
    if [[ ! -d "$input_dir" ]]; then
        echo
        write_log "No input directory found" "ERROR"
        echo
        help
        exit 1
    fi
    if [[ -z "$output_dir" ]]; then
        echo
        write_log "No output directory specified" "ERROR"
        echo
        exit 1
    fi
    if [[ ! -r "$config_file" ]]; then
        echo
        write_log "Configuration file not readable or does not exist" "ERROR"
        echo
        exit 1
    fi
    write_log "input validation passed" "INFO"
}

create_output_structure() {
    local output_dir="$1"
    local log_file="$2"

    mkdir -p "$output_dir"
    mkdir -p "$output_dir/processing_files"
    mkdir -p "$output_dir/500_bpTrimmed_fastas"
    mkdir -p "$output_dir/compiled_files"
    mkdir -p "$output_dir/logs"
    mkdir -p "$output_dir/tmp"

    write_log "wrote: $output_dir" "INFO" "$log_file"
    write_log "wrote: $output_dir/processing_files" "INFO" "$log_file"
    write_log "wrote: $output_dir/500_bpTrimmed_fastas" "INFO" "$log_file"
    write_log "wrote: $output_dir/compiled_files" "INFO" "$log_file"
    write_log "wrote: $output_dir/logs" "INFO" "$log_file"
    write_log "wrote: $output_dir/tmp" "INFO" "$log_file"
}

#loads into variables from the configuration file
load_config() {
    local config_file="$1"
    local log_file="$2"
    
    pipeline_dir="$(grep '^source_directory=' "$config_file" | awk -F'=' '{print $2}')"
    conda_env_prefix="$(grep '^conda_env_prefix=' "$config_file" | awk -F'=' '{print $2}')"
    max_sequence_lengths="$(grep '^max_seq_lengths_array=' "$config_file" | awk -F'=' '{print $2}')"
    coverage_modes="$(grep '^cov_modes_array=' "$config_file" | awk -F'=' '{print $2}')"
    execution_mode="$(grep '^mode=' "$config_file" | awk -F'=' '{print $2}')"

    #if slurm mode, load slurm-specific settings
    if [[ $execution_mode == "slurm" ]]; then
        slurm_cpus_per_job="$(grep '^slurm_cpus_per_job=' "$config_file" | awk -F'=' '{print $2}')"
        slurm_memory_per_job="$(grep '^slurm_memory_per_job=' "$config_file" | awk -F'=' '{print $2}')"
        slurm_partition="$(grep '^slurm_partition=' "$config_file" | awk -F'=' '{print $2}')"
        max_jobs_per_array="$(grep '^max_jobs_per_array=' "$config_file" | awk -F'=' '{print $2}')"
        write_log "slurm_cpus_per_job=$slurm_cpus_per_job" "INFO" "$log_file"
        write_log "slurm_memory_per_job=$slurm_memory_per_job" "INFO" "$log_file"
        write_log "slurm_partition=$slurm_partition" "INFO" "$log_file"
        write_log "max_jobs_per_array=$max_jobs_per_array" "INFO" "$log_file"
    fi

    #defining non config variables
    reference_fasta_file="$pipeline_dir/database/elementgeneList.fasta"

    write_log "pipeline_dir=$pipeline_dir" "INFO" "$log_file"
    write_log "conda_env_prefix=$conda_env_prefix" "INFO" "$log_file"
    write_log "max_sequence_lengths=$max_sequence_lengths" "INFO" "$log_file"
    write_log "coverage_modes=$coverage_modes" "INFO" "$log_file"
    write_log "execution_mode=$execution_mode" "INFO" "$log_file"
    write_log "reference_fasta_file=$reference_fasta_file" "INFO" "$log_file"
}

#sample ID "xxxx.fasta" per line
write_sample_id_list() {
    local input_dir="$1"
    local sample_id_list_dir="$2"
    local log_file="$3"
    local sample_id_list_file="$sample_id_list_dir/sample_ID_list.txt"
    local fasta_pattern="*.f*"

    find -L "$input_dir" -maxdepth 1 -name "$fasta_pattern" -exec basename {} ';' > "$sample_id_list_file"
    local command_exit_status=$?
    if [[ $command_exit_status -eq 0 ]]; then
        write_log "wrote sample ID list to $sample_id_list_file" "INFO" "$log_file"
        write_log "wrote $(wc -l < "$sample_id_list_file") sample IDs to $sample_id_list_file" "INFO" "$log_file"
    else
        write_log "Failed to write sample ID list to $sample_id_list_file" "ERROR" "$log_file"
    fi

}

#less than 500 bp sequences removal
remove_smalls() {
    local remove_smalls_script_file="$1"
    local input_fasta_dir="$2"
    local sample_id_list_file="$3"
    local trimmed_fasta_dir="$4"
    local log_file="$5"

    while read -r sample_filename; do
        perl "$remove_smalls_script_file" 500 "$input_fasta_dir/$sample_filename" > "$trimmed_fasta_dir/$sample_filename"
    done < "$sample_id_list_file"

    write_log "remove_smalls script processed $(ls "$trimmed_fasta_dir" | wc -l) files" "INFO" "$log_file"
}

parallel_mmseqs_searches() {
    # summary:
    #   runs [list_of_max_seq_lengths * list_of_cov_modes] mmseqs searches for each isolate.
    #
    # Arguments:
    #   1: conda_env_prefix        <path>   : path to the conda environment prefix
    #   2: processing_dir       <path>   : path to the processing files directory (should contain isolate directories with "query_nucl_db" subdirectories)
    #   3: reference_db_prefix  <path>   : path to the nucleotide reference database prefix
    #   4: max_sequence_lengths <STRING> : space-separated maximum sequence lengths
    #   5: coverage_modes       <STRING> : space-separated coverage modes
    #   6: log_file                <path>   : (optional)
    # Notes:
    #   Should decouple mmseqs and conversion steps

    local conda_env_prefix="$1"
    local processing_dir="$2"
    local reference_db_prefix="$3"
    local max_sequence_lengths="$4"
    local coverage_modes="$5"
    local log_file="$6"

    local -a max_sequence_lengths_array=()
    local -a coverage_modes_array=()
    read -r -a max_sequence_lengths_array <<< "$max_sequence_lengths"
    read -r -a coverage_modes_array <<< "$coverage_modes"

    write_log "starting isolate-level parallel mmseqs with conda bin $conda_env_prefix/bin/parallel" "INFO" "$log_file"

    "$conda_env_prefix/bin/parallel" --jobs "${SLURM_CPUS_PER_TASK:-1}" \
                                    --joblog parallel.joblog \
                                    run_mmseqs_search_and_convert "$conda_env_prefix" \
                                                                "{1}/query_nucl_db/{1/}_nucl_db_type_2" \
                                                                "{2}" \
                                                                "{1}/results_db" \
                                                                "{1}/tmp" \
                                                                "{3}" \
                                                                "{4}" \
                                                                "$log_file" \
                                                                ::: "$processing_dir"/* \
                                                                ::: "$reference_db_prefix" \
                                                                ::: "${coverage_modes_array[@]}" \
                                                                ::: "${max_sequence_lengths_array[@]}"
}

#combine mmseqs search results per isolate
combine_mmseqs_results() {
    # combines n=6 mmseqs convertalis tsv files into a single file
    # python script removes duplicate entries from the combined results
    #args:
    local conda_env_prefix="$1"
    local processing_dir="$2"
    local reference_fasta_file="$3"
    local results_combiner_script_file="$4"

    local reference_sequence_count=
    local sample_id=
    local results_dir=

    reference_sequence_count=$(grep -c "^>" "$reference_fasta_file")

    for sample_dir in "$processing_dir"/*; do
        sample_id=$(basename "$sample_dir")
        results_dir="$sample_dir/results_db"

        conda run -p "$conda_env_prefix" python3 "$results_combiner_script_file" \
                                                "$results_dir" \
                                                "$reference_fasta_file" \
                                                "$reference_sequence_count" \
                                                "$sample_dir/$sample_id"
    done
}

write_slurm_array_file() {
    local sample_id_list_file="$1"

}

#######################################
############# run script ##############
#######################################

#global variables, argparsing and input validation
input_dir=""
output_dir=""
config_file=""
while getopts ":i:o:c:h" opt; do
    case "$opt" in
        i) input_dir="$OPTARG";;
        o) output_dir="$OPTARG";;
        c) config_file="$OPTARG";;
        h) help; exit 0;;
        *) help; exit 1;;
    esac
done

validate_input "$input_dir" "$output_dir" "$config_file"

#setup
create_output_structure "$output_dir" "$output_dir/logs/run.log"
load_config "$config_file" "$output_dir/logs/run.log"
write_version_info "$conda_env_prefix" "$output_dir/logs/run.log"
write_sample_id_list "$input_dir" "$output_dir" "$output_dir/logs/run.log"

#source mmseqs functionality
source "$pipeline_dir/subscripts/mmseqs_functionality.sh"

#removing sequences shorter than 500bp
remove_smalls "$pipeline_dir/subscripts/removesmalls.pl" \
              "$input_dir" \
              "$output_dir/sample_ID_list.txt" \
              "$output_dir/500_bpTrimmed_fastas" \
              "$output_dir/logs/run.log"

#creating mmseqs databases for nucleotide sequences: 1 shared ref and 1 query for each isolate
write_nucl_reference_db "$conda_env_prefix" \
                        "$reference_fasta_file" \
                        "$output_dir/tmp/reference_db" \
                        "$output_dir/logs/run.log"

#run analysis (local or slurm)
case "$execution_mode" in
    local)
    write_log "Starting $execution_mode mmseqsmode" "INFO" "$output_dir/logs/run.log"
    #initiate search of (max_seq_lengths_array * cov_modes_array) combinations
    for trimmed_fasta in "$output_dir/500_bpTrimmed_fastas"/*; do

        sample_filename=$(basename "$trimmed_fasta")
        sample_id="${sample_filename%.*}"
        query_database_prefix="$output_dir/processing_files/$sample_id/query_nucl_db/${sample_id}_nucl_db_type_2"

        write_nucl_query_db "$conda_env_prefix" \
                            "$trimmed_fasta" \
                            "$query_database_prefix" \
                            "$output_dir/logs/run.log"
    done

    # exporting functions for GNU parallel
    export -f run_mmseqs_search_and_convert write_log
    parallel_mmseqs_searches "$conda_env_prefix" \
                            "$output_dir/processing_files" \
                            "$reference_db_prefix" \
                            "$max_sequence_lengths" \
                            "$coverage_modes" \
                            "$output_dir/logs/run.log"

    #combine the mmseqs search results per isolate
    combine_mmseqs_results "$conda_env_prefix" \
                        "$output_dir/processing_files" \
                        "$reference_fasta_file" \
                        "$pipeline_dir/subscripts/mmseq2_results_replicate_combine.py"
    ;;
    slurm)
    write_log "Starting $execution_mode mmseqs mode" "INFO" "$output_dir/logs/run.log"
    source "$pipeline_dir/subscripts/slurm_functionality.sh"

    write_manifest_file "$output_dir/processing_files" \
                        "$output_dir/500_bpTrimmed_fastas" \
                        "$output_dir/sample_ID_list.txt" \
                        "$reference_db_prefix" \
                        "$coverage_modes" \
                        "$max_sequence_lengths" \
                        "$output_dir/manifest.csv"

    #initate slurm runners
    start_slurm_runners "$output_dir/manifest.csv" \
                        "$max_jobs_per_array" \
                        "$slurm_cpus_per_job" \
                        "$slurm_memory_per_job" \
                        "$slurm_partition" \
                        "$pipeline_dir/subscripts/slurm_runner_worker.sh" 
    ;;
    *) write_log "Invalid mode: $execution_mode" "ERROR" "$output_dir/logs/run.log"; exit 1 ;;
esac
