#!/usr/bin/env bash
#SBATCH --job-name=HEP_analysis
#SBATCH --output=HEP_analysis.out
#SBATCH --error=HEP_analysis.err
#SBATCH --partition=standard
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
    echo "Optional:"
    echo "  -f <host_file>       Host file, tsv seperated, containing Genome_Ref and Host"
}

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

write_version_info() {
    local conda_env_prefix="$1"
    local log_file="${2:-}"
    
    write_log "mmseqs version: $(conda run -p "$conda_env_prefix" mmseqs version)" "INFO" "$log_file"
}

#######################################
############# main functions ##########
#######################################

validate_input() {
    local input_dir="$1"
    local output_dir="$2"
    local config_file="$3"
    local host_file="$4"
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
        help
        exit 1
    fi
    if [[ ! -r "$config_file" ]]; then
        echo
        write_log "Configuration file not readable or does not exist" "ERROR"
        echo
        help
        exit 1
    fi
    if [[ ! -r "$host_file" ]]; then
        echo
        write_log "Host .tsv file not provided. Assigning host from config" "WARNING"
        echo
        sleep 3
    fi
    if [[ -f "$host_file" ]]; then
        local header_line=$'Genome_Ref\tHost'
        if [[ $(grep -F "$header_line" "$host_file") ]]; then
            write_log "Host file contains the required header: $header_line" "INFO"
        else
            write_log "Host file missing the required header: $header_line" "ERROR"
            exit 1
        fi
    fi
    

}

create_output_structure() {
    local output_dir="$1"
    local log_file="${2:-}"

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
    local log_file="${2:-}"
    
    pipeline_dir="$(grep '^source_directory=' "$config_file" | awk -F'=' '{print $2}')"
    conda_env_prefix="$(grep '^conda_env_prefix=' "$config_file" | awk -F'=' '{print $2}')"
    reference_fasta_file="$(grep '^reference_fasta_file=' "$config_file" | awk -F'=' '{print $2}')"
    max_sequence_lengths="$(grep '^max_seq_lengths_array=' "$config_file" | awk -F'=' '{print $2}')"
    coverage_modes="$(grep '^cov_modes_array=' "$config_file" | awk -F'=' '{print $2}')"
    execution_mode="$(grep '^execution_mode=' "$config_file" | awk -F'=' '{print $2}')"
    base_host="$(grep '^base_host=' "$config_file" | awk -F'=' '{print $2}')"
    fasta_pattern="$(grep '^fasta_pattern=' "$config_file" | awk -F'=' '{print $2}')"
    mmseqs_min_seq_id="$(grep '^mmseqs_min_seq_id=' "$config_file" | awk -F'=' '{print $2}')"
    mmseqs_coverage="$(grep '^mmseqs_coverage=' "$config_file" | awk -F'=' '{print $2}')"
    mmseqs_sensitivity="$(grep '^mmseqs_sensitivity=' "$config_file" | awk -F'=' '{print $2}')"
    mmseqs_max_seqs="$(grep '^mmseqs_max_seqs=' "$config_file" | awk -F'=' '{print $2}')"

    #if slurm mode, load slurm-specific settings
    if [[ $execution_mode == "slurm" ]]; then
        slurm_cpus_per_job="$(grep '^slurm_cpus_per_job=' "$config_file" | awk -F'=' '{print $2}')"
        slurm_memory_per_job="$(grep '^slurm_memory_per_job=' "$config_file" | awk -F'=' '{print $2}')"
        slurm_partition="$(grep '^slurm_partition=' "$config_file" | awk -F'=' '{print $2}')"
        max_jobs_per_array="$(grep '^max_jobs_per_array=' "$config_file" | awk -F'=' '{print $2}')"
        max_parallel_jobs_per_array="$(grep '^max_parallel_jobs_per_array=' "$config_file" | awk -F'=' '{print $2}')"
        write_log "slurm_cpus_per_job=$slurm_cpus_per_job" "INFO" "$log_file"
        write_log "slurm_memory_per_job=$slurm_memory_per_job" "INFO" "$log_file"
        write_log "slurm_partition=$slurm_partition" "INFO" "$log_file"
        write_log "max_jobs_per_array=$max_jobs_per_array" "INFO" "$log_file"
        write_log "max_parallel_jobs_per_array=$max_parallel_jobs_per_array" "INFO" "$log_file"
    fi

    #defining non config variables
    if [[ -z "$reference_fasta_file" ]]; then
        reference_fasta_file="$pipeline_dir/database/elementgeneList.fasta"
    fi

    write_log "_______________configurations_______________" "INFO" "$log_file"
    write_log "pipeline_dir=$pipeline_dir" "INFO" "$log_file"
    write_log "conda_env_prefix=$conda_env_prefix" "INFO" "$log_file"
    write_log "max_sequence_lengths=$max_sequence_lengths" "INFO" "$log_file"
    write_log "coverage_modes=$coverage_modes" "INFO" "$log_file"
    write_log "execution_mode=$execution_mode" "INFO" "$log_file"
    write_log "base_host=$base_host" "INFO" "$log_file"
    write_log "fasta_pattern=$fasta_pattern" "INFO" "$log_file"
    write_log "mmseqs_min_seq_id=$mmseqs_min_seq_id" "INFO" "$log_file"
    write_log "mmseqs_coverage=$mmseqs_coverage" "INFO" "$log_file"
    write_log "mmseqs_sensitivity=$mmseqs_sensitivity" "INFO" "$log_file"
    write_log "mmseqs_max_seqs=$mmseqs_max_seqs" "INFO" "$log_file"
    write_log "reference_fasta_file=$reference_fasta_file" "INFO" "$log_file"
    write_log "_______________configurations_______________" "INFO" "$log_file"
}

validate_config() {
    local log_file="${1:-}"
    local config_values=()

    config_values=(
        pipeline_dir
        conda_env_prefix
        reference_fasta_file
        execution_mode
        fasta_pattern
        base_host
        mmseqs_min_seq_id
        mmseqs_coverage
        mmseqs_sensitivity
        mmseqs_max_seqs
        coverage_modes
        max_sequence_lengths
    )
    if [[ -z $execution_mode ]]; then
        write_log "Configuration value execution_mode is not set" "ERROR" "$log_file"
        exit 1
    fi
    if [[ $execution_mode == "slurm" ]]; then
        config_values+=(
            slurm_cpus_per_job
            slurm_memory_per_job
            slurm_partition
            max_jobs_per_array
            max_parallel_jobs_per_array
        )
    fi
    for config_value in "${config_values[@]}"; do
        if [[ -z "${!config_value}" ]]; then
            write_log "Configuration value $config_value is not set" "ERROR" "$log_file"
            exit 1
        fi
    done
}

#sample ID "xxxx.fasta" per line
write_sample_id_list() {
    local input_dir="$1"
    local sample_id_list_dir="$2"
    local fasta_pattern="$3"
    local log_file="${4:-}"
    local sample_id_list_file="$sample_id_list_dir/sample_ID_list.txt"

    find -L "$input_dir" -maxdepth 1 -name "$fasta_pattern" -exec basename {} ';' > "$sample_id_list_file"
    local command_exit_status=$?
    if [[ $command_exit_status -eq 0 ]]; then
        write_log "wrote sample ID list to $sample_id_list_file" "INFO" "$log_file"
        write_log "wrote $(wc -l < "$sample_id_list_file") sample IDs to $sample_id_list_file" "INFO" "$log_file"
    else
        write_log "Failed to write sample ID list to $sample_id_list_file" "ERROR" "$log_file"
    fi

}

write_host_file() {
    local base_host="$1"
    local sample_list="$2"
    local host_file_name="$3"
    local log_file="${4:-}"
    
    : > "$host_file_name"
    #header
    echo -e "Genome_Ref\tHost" > "$host_file_name"
    #loop over sample_list
    local sample_name
    while read -r sample; do
        sample_name=${sample%.*}
        echo -e "${sample_name}\t${base_host}" >> "$host_file_name"
    done < "$sample_list"

    write_log "Host file written to $host_file_name" "INFO" "$log_file"
}

#less than 500 bp sequences removal
remove_smalls() {
    local remove_smalls_script_file="$1"
    local input_fasta_dir="$2"
    local sample_id_list_file="$3"
    local trimmed_fasta_dir="$4"
    local log_file="${5:-}"

    while read -r sample_filename; do
        perl "$remove_smalls_script_file" 500 "$input_fasta_dir/$sample_filename" > "$trimmed_fasta_dir/$sample_filename"
    done < "$sample_id_list_file"

    write_log "remove_smalls script processed $(ls "$trimmed_fasta_dir" | wc -l) files" "INFO" "$log_file"
}



#######################################
############# run script ##############
#######################################
main() {
#global variables, argparsing and input validation
input_dir=""
output_dir=""
config_file=""
host_file=""
while getopts ":i:o:c:f:h" opt; do
    case "$opt" in
        i) input_dir="$OPTARG";;
        o) output_dir="$OPTARG";;
        c) config_file="$OPTARG";;
        f) host_file="$OPTARG";;
        h) help; exit 0;;
        *) help; exit 1;;
    esac
done

validate_input "$input_dir" "$output_dir" "$config_file" "$host_file"

#setup
create_output_structure "$output_dir" "$output_dir/logs/run.log"
load_config "$config_file" "$output_dir/logs/run.log"
validate_config "$output_dir/logs/run.log"
write_version_info "$conda_env_prefix" "$output_dir/logs/run.log"
write_sample_id_list "$input_dir" "$output_dir" "$fasta_pattern" "$output_dir/logs/run.log"
#create default host file if not provided
if [[ -z "$host_file" ]]; then
    write_log "Host file is not specified" "WARNING" "$output_dir/logs/run.log"
    sleep 2
    write_log "writing default host file with base_host=$base_host" "INFO" "$output_dir/logs/run.log"
    
    write_host_file "$base_host" \
                    "$output_dir/sample_ID_list.txt" \
                    "$output_dir/host_file.tsv" \
                    "$output_dir/logs/run.log"
    host_file="$output_dir/host_file.tsv"
    write_log "wrote host file: $(head -n 5 "$host_file")" "INFO" "$output_dir/logs/run.log"
fi


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
    write_log "Starting $execution_mode mode" "INFO" "$output_dir/logs/run.log"
    
    #TODO: encapsulate in a funtion later maybe?
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
                            "$mmseqs_min_seq_id" \
                            "$mmseqs_coverage" \
                            "$mmseqs_sensitivity" \
                            "$mmseqs_max_seqs" \
                            "$output_dir/logs/run.log"

    #combine the mmseqs search results per isolate
    for sample_dir in "$output_dir/processing_files"/*; do
        combine_mmseqs_results_per_isolate "$conda_env_prefix" \
                              "$sample_dir/results_db" \
                              "$reference_fasta_file" \
                              "$pipeline_dir/subscripts/mmseq2_results_replicate_combine.py" \
                              "$output_dir/logs/run.log"
    done

    # compile all per-isolate mmseqs results into 1 file
    compile_mmseqs_results \
        "$output_dir/processing_files" \
        "$output_dir/compiled_files" \
        "$output_dir/logs/run.log"

    run_host_element_screen_processor "$conda_env_prefix" \
                                     "$pipeline_dir/subscripts/host_element_screen_processor.py" \
                                     "$output_dir/compiled_files/result_compiled" \
                                     "$host_file" \
                                     "$reference_fasta_file" \
                                     "$output_dir/compiled_files/result_compiled"
    
    write_log "Finished local pipeline" "INFO" "$output_dir/logs/run.log"
    write_log " $(wc -l < "$output_dir/compiled_files/mmseq2_result_presence_absence.tsv") isolates compiled" "INFO" "$output_dir/logs/run.log"
    
    
    ;;
    slurm)


    write_log "Starting $execution_mode mode" "INFO" "$output_dir/logs/run.log"
    source "$pipeline_dir/subscripts/slurm_functionality.sh"

    #writes a metafile, which contains information about the samples and the parameters for the SLURM jobs
    write_slurm_meta_info_file "$output_dir/processing_files" \
                        "$output_dir/500_bpTrimmed_fastas" \
                        "$output_dir/sample_ID_list.txt" \
                        "$reference_db_prefix" \
                        "$coverage_modes" \
                        "$max_sequence_lengths" \
                        "$mmseqs_min_seq_id" \
                        "$mmseqs_coverage" \
                        "$mmseqs_sensitivity" \
                        "$mmseqs_max_seqs" \
                        "$max_jobs_per_array" \
                        "$output_dir/slurm_meta_info.csv" \
                        "$output_dir/logs/run.log"

    #initate slurm runners
    start_slurm_runners "$conda_env_prefix" \
                        "$pipeline_dir" \
                        "$output_dir/slurm_meta_info.csv" \
                        "$max_jobs_per_array" \
                        "$max_parallel_jobs_per_array" \
                        "$slurm_cpus_per_job" \
                        "$slurm_memory_per_job" \
                        "$slurm_partition" \
                        "$pipeline_dir/subscripts/slurm_runner_worker.sh" \
                        "$reference_fasta_file" \
                        "$output_dir/logs/run.log"

    start_slurm_compiler "$conda_env_prefix" \
                         "$pipeline_dir" \
                         "$output_dir" \
                         "$host_file" \
                         "$reference_fasta_file" \
                         "$slurm_memory_per_job" \
                         "$slurm_partition" \
                         "$pipeline_dir/subscripts/slurm_compiler_worker.sh" \
                         "$output_dir/logs/run.log"
    ;;
    *) write_log "Invalid mode: $execution_mode" "ERROR" "$output_dir/logs/run.log"; exit 1 ;;
esac
}
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi