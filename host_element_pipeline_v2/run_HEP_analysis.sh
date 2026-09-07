#!/bin/bash

#script author Jon Slotved (JOSS@dksund.dk)

#######################################
############# helper functions ########
#######################################

help() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -c <config_file> [-h]"
    echo "  -i <input_dir>       Input directory containing FASTA files"
    echo "  -o <output_dir>      Output directory for results and logs"
    echo "  -c <config_file>     Configuration file specifying source directory and conda environment prefix"
    echo "  -h                   Display this help message"
}

write_log() {
    local log_message="${1:-No log message provided}"
    local log_type="${2:-INFO}"
    local log_file="$3"
    local time=$(date +"%Y-%m-%d %H:%M:%S")
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
    mkdir -p "$output_dir/results"
    mkdir -p "$output_dir/logs"
    mkdir -p "$output_dir/tmp"

    write_log "wrote: $output_dir" "INFO" "$log_file"
    write_log "wrote: $output_dir/processing_files" "INFO" "$log_file"
    write_log "wrote: $output_dir/500_bpTrimmed_fastas" "INFO" "$log_file"
    write_log "wrote: $output_dir/results" "INFO" "$log_file"
    write_log "wrote: $output_dir/logs" "INFO" "$log_file"
    write_log "wrote: $output_dir/tmp" "INFO" "$log_file"
}

load_config() {
    local config_file="$1"
    local log_file="$2"
    source_directory=""
    conda_env_prefix=""
    reference_gene_list_name=""
    
    source_directory="$(grep '^source_directory=' "$config_file" | awk -F'=' '{print $2}')"
    conda_env_prefix="$(grep '^conda_env_prefix=' "$config_file" | awk -F'=' '{print $2}')"
    reference_gene_list_name="$(grep '^reference_gene_list_name=' "$config_file" | awk -F'=' '{print $2}')"

    write_log "source_directory=$source_directory" "INFO" "$log_file"
    write_log "conda_env_prefix=$conda_env_prefix" "INFO" "$log_file"
    write_log "reference_gene_list_name=$reference_gene_list_name" "INFO" "$log_file"
}

write_sample_ID_list() {
    local input_dir="$1"
    local sample_list_output_dir="$2"
    local log_file="$3"
    local id_list_file="$sample_list_output_dir/sample_ID_list.txt"
    local pattern="*.f*"

    find -L "$input_dir" -maxdepth 1 -name "$pattern" -exec basename {} ';' > "$id_list_file"

    write_log "Sample ID list 'based on $pattern pattern' written to $id_list_file" "INFO" "$log_file"
    write_log "wrote $(wc -l < "$id_list_file") sample IDs to $id_list_file" "INFO" "$log_file"
}

remove_smalls() {
    local remove_smalls_script="$1"
    local input_dir="$2"
    local sample_list="$3"
    local smalls_output_dir="$4"
    local log_file="$5"

    while read -r line; do
        perl "$remove_smalls_script" 500 "$input_dir/$line" > "$smalls_output_dir/$line"
    done < "$sample_list"

    write_log "remove_smalls script processed $(ls "$smalls_output_dir" | wc -l) files" "INFO" "$log_file"
}

write_reference_db() {
    local conda_env_prefix="$1"
    local reference_gene_list="$2"
    local db_output_dest="$3"
    local db_type="$4"
    local log_file="$5"

    mkdir -p "$db_output_dest"
    
    # db_type = 0 for protein, 1 for nucleotide
    conda run -p "$conda_env_prefix" mmseqs createdb "$reference_gene_list" \
                                                     "$db_output_dest/reference_dbtype_${db_type}" \
                                                     --dbtype "$db_type"
    
    if [[ -d "$db_output_dest" ]]; then
        write_log "write ref db {dbtype=$db_type}: $db_output_dest" "INFO" "$log_file"
    else
        write_log "failed to write ref db {dbtype=$db_type}: $db_output_dest" "ERROR" "$log_file"
    fi  
}

write_query_dbs() {
    local conda_env_prefix="$1"
    local input_dir="$2"
    local sample_list="$3"
    local db_output_dest="$4"
    local log_file="$5"
    local sample_ID=""

    while read -r line; do
        #remove extension and write dir
        sample_ID="${line%.*}"
        mkdir -p "$db_output_dest/${sample_ID}"

        conda run -p "$conda_env_prefix" mmseqs createdb "$input_dir/$line" \
                                                     "$db_output_dest/${sample_ID}/${sample_ID}_db" > /dev/null 2>&1
        if [[ -d "$db_output_dest/$sample_ID" ]]; then
            write_log "wrote query db: $sample_ID" "INFO"
        else
            write_log "Failed to write db: $sample_ID" "ERROR"
        fi

    done < "$sample_list"

    write_log "$(ls "$db_output_dest" | wc -l) directories generated in $db_output_dest" "INFO" "$log_file"
}




#######################################
############# run script ##############
#######################################

#default values, argparsing and input validation
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
write_sample_ID_list "$input_dir" "$output_dir" "$output_dir/logs/run.log"

#removing sequences shorter than 500bp
remove_smalls "$source_directory/subscripts/removesmalls.pl" "$input_dir" "$output_dir/sample_ID_list.txt" "$output_dir/500_bpTrimmed_fastas" "$output_dir/logs/run.log"

#creating mmseqs databases 1 ref n query
write_reference_db "$conda_env_prefix" "$source_directory/database/$reference_gene_list_name" "$output_dir/tmp/reference_db" 1 "$output_dir/logs/run.log"
write_query_dbs "$conda_env_prefix" "$output_dir/500_bpTrimmed_fastas" "$output_dir/sample_ID_list.txt" "$output_dir/processing_files" "$output_dir/logs/run.log"

#initiate analysis
