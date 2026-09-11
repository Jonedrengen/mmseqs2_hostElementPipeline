#!/usr/bin/env bash
#SBATCH --job-name=HEP_analysis
#SBATCH --output=HEP_analysis.log
#SBATCH --error=HEP_analysis.err
#SBATCH --partition=project
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

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
    source_directory=""
    conda_env_prefix=""
    reference_gene_list_name=""
    max_seq_lengths_array=""
    cov_modes_array=""
    
    source_directory="$(grep '^source_directory=' "$config_file" | awk -F'=' '{print $2}')"
    conda_env_prefix="$(grep '^conda_env_prefix=' "$config_file" | awk -F'=' '{print $2}')"
    reference_gene_list_name="$(grep '^reference_gene_list_name=' "$config_file" | awk -F'=' '{print $2}')"
    max_seq_lengths_array="$(grep '^max_seq_lengths_array=' "$config_file" | awk -F'=' '{print $2}')" 
    cov_modes_array="$(grep '^cov_modes_array=' "$config_file" | awk -F'=' '{print $2}')"

    write_log "source_directory=$source_directory" "INFO" "$log_file"
    write_log "conda_env_prefix=$conda_env_prefix" "INFO" "$log_file"
    write_log "reference_gene_list_name=$reference_gene_list_name" "INFO" "$log_file"
    write_log "max_seq_lengths_array=$max_seq_lengths_array" "INFO" "$log_file"
    write_log "cov_modes_array=$cov_modes_array" "INFO" "$log_file"
}

#sample ID "xxxx.fasta" per line
write_sample_ID_list() {
    local input_dir="$1"
    local sample_list_output_dir="$2"
    local log_file="$3"
    local id_list_file="$sample_list_output_dir/sample_ID_list.txt"
    local pattern="*.f*"

    find -L "$input_dir" -maxdepth 1 -name "$pattern" -exec basename {} ';' > "$id_list_file"
    local exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "wrote sample ID list to $id_list_file" "INFO" "$log_file"
        write_log "wrote $(wc -l < "$id_list_file") sample IDs to $id_list_file" "INFO" "$log_file"
    else
        write_log "Failed to write sample ID list to $id_list_file" "ERROR" "$log_file"
    fi

}

#less than 500 bp sequences removal
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

write_nucl_reference_db() {
    local conda_env_prefix="$1"
    local reference_gene_list="$2"
    local db_output_dest="$3"
    local log_file="$4"
    local db_type=2
    # ref db needed globally
    reference_db=""


    mkdir -p "$db_output_dest"
    reference_db="$db_output_dest/reference_nucl_db_type_${db_type}"
    write_log "reference_db=$reference_db" "INFO" "$log_file"

    # db_type = 0 for protein, 1 for nucleotide
    conda run -p "$conda_env_prefix" mmseqs createdb "$reference_gene_list" \
                                                     "$reference_db" \
                                                     --dbtype $db_type > /dev/null
    
    local exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "wrote ref db {dbtype=nucl}: $db_output_dest" "INFO" "$log_file"
    else
        write_log "failed to write ref db {dbtype=nucl}: $db_output_dest" "ERROR" "$log_file"
    fi
}

write_nucl_query_dbs() {
    local conda_env_prefix="$1"
    local input_dir="$2"
    local sample_list="$3"
    local db_output_dir="$4"
    local log_file="$5"
    local sample_ID=""
    local db_type=2
    
    # iterates over each sample in the sample list and create a query database for it in processing_files
    while read -r line; do
        #remove extension and write dir
        sample_ID="${line%.*}"
        mkdir -p "$db_output_dir/${sample_ID}/query_nucl_db"

        #make db (nucl = --dbtype 2)
        conda run -p "$conda_env_prefix" mmseqs createdb "$input_dir/$line" \
                                 "$db_output_dir/${sample_ID}/query_nucl_db/${sample_ID}_nucl_db_type_${db_type}" \
                                 --dbtype $db_type > /dev/null
        
        #log
        local exit_status=$?
        if [[ $exit_status -eq 0 ]]; then
            write_log "wrote query db nucl type_${db_type}: ${sample_ID} " "INFO"
        else
            write_log "Failed to write db nucl type_${db_type}: $sample_ID" "ERROR" "$log_file"
        fi

    done < "$sample_list"

    write_log "$(ls "$db_output_dir" | wc -l) directories generated in $db_output_dir" "INFO" "$log_file"
}

# controlled by parallel_mmseqs_searches function
run_mmseqs_search_and_convert() {
    # summary:
    #   Runs an mmseqs search and converts the results to a human-readable format.
    #   This function is not inteded to be called directly by users.
    #   This function is used internally to orchestrate mmseqs searches and conversions.
    # author: 
    #   Jon Slotved 09-09-2026
    # Arguments:
    #   1: conda_env_prefix <path>
    #   2: query_db         <path>
    #   3: reference_db     <path>
    #   4: results_db       <path>
    #   5: tmp_db           <path>
    #   6: search_type      <INT>   : made for and tested on 3 (nucleotide search)
    #   7: cov_mode         <INT>   : 1 (target/ref cov), 2 (query cov).
    #   8: max_seq_length   <INT>   : maximum sequence length for the search.
    #   9: log_file         <path>  : (optional)

    local conda_env_prefix="$1"
    
    #databases
    local query_db="$2"
    local reference_db="$3"
    local results_db="$4"
    local tmp_db="$5"

    # mmseqs search parameters
    local search_type="$6"
    local cov_mode="$7"
    local max_seq_length="$8"

    #log
    local log_file="$9"

    #get the sample name
    local sample_name=""
    sample_name=$(basename "$query_db")

    #search and convert outputs
    local mmseqs_search_output="$results_db/${sample_name}_db_cov_${cov_mode}_max_len_${max_seq_length}"
    local mmseqs_convert_output="$results_db/${sample_name}_db_cov_${cov_mode}_max_len_${max_seq_length}.tsv"

    # because mmseqs requires the results directory to exist beforehand
    mkdir -p "$results_db"
    mkdir -p "$tmp_db"

    #search
    conda run -p "$conda_env_prefix" mmseqs search "$query_db" \
                                                   "$reference_db" \
                                                   "$mmseqs_search_output" \
                                                   "$tmp_db" \
                                                   --search-type "$search_type" \
                                                   --cov-mode "$cov_mode" \
                                                   --max-seq-len "$max_seq_length" \
                                                   --min-seq-id 0.8 \
                                                   -c 0.8 \
                                                   -s 7.5 \
                                                   --max-seqs 1000 \
                                                   --threads 1 \
                                                   -a \
                                                   --mask 0 \
                                                   --comp-bias-corr 0 \
                                                   --strand 2 > /dev/null

    #logs
    local exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "mmseqs search completed: $sample_name against $(basename "$reference_db")" "INFO"
    else
        write_log "mmseqs search failed" "ERROR" "$log_file"
        exit $exit_status
    fi

    #convert to human readable
    conda run -p "$conda_env_prefix" mmseqs convertalis "$query_db" \
                                                       "$reference_db" \
                                                       "$mmseqs_search_output" \
                                                       "$mmseqs_convert_output" \
                                                       --format-output query,target,pident,qcov,tcov,alnlen,mismatch,gapopen,qlen,qstart,qend,tlen,tstart,tend,evalue,bits,cigar \
                                                       --search-type "$search_type" \
                                                       --threads 1 > /dev/null

    exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "mmseqs convertalis completed: $sample_name/$(basename "$mmseqs_convert_output")" "INFO"
    else
        write_log "mmseqs convertalis failed" "ERROR" "$log_file"
        exit $exit_status
    fi
    
    #add header. This is for compatibility with Edwards format (should change later -Jon) 
    #quick fix... Should change later
    cp "$mmseqs_convert_output" "$mmseqs_convert_output.tmp"
    printf "Query_Seq-id\tSubject_Seq-id\tPercent_Identity\tQuery_Coverage\tSubject_Coverage\tAlignment_Length\tMismatches\tGapOpenings\tQuery_Length\tQuery_Start\tQuery_End\tSubject_Length\tSubject_Start\tSubject_End\tE-Value\tBitscore\tCigar\n" > "$mmseqs_convert_output"
    cat "$mmseqs_convert_output.tmp" >> "$mmseqs_convert_output"
    rm -f "$mmseqs_convert_output.tmp"
}

#runs n=6 mmseqs searches per isolate found in processing_files
parallel_mmseqs_searches() {
    # summary:
    #   runs [list_of_max_seq_lengths * list_of_cov_modes] mmseqs searches for each isolate.
    #
    # Arguments:
    #   1: conda_env_prefix        <path>   : path to the conda environment prefix
    #   2: processing_files_dir    <path>   : path to the processing files directory (should contain isolate directories with "query_nucl_db" subdirectories)
    #   3: reference_db            <path>   : path to the reference database
    #   4: list_of_max_seq_lengths <STRING> : space-separated maximum sequence lengths
    #   5: list_of_cov_modes       <STRING> : space-separated coverage modes
    #   6: log_file                <path>   : (optional)

    local conda_env_prefix="$1"
    local processing_files_dir="$2"
    local reference_db="$3"
    local list_of_max_seq_lengths="$4"
    local list_of_cov_modes="$5"
    local log_file="$6"
    
    local max_seq_array=()
    read -r -a max_seq_array <<< "$list_of_max_seq_lengths"
    local cov_modes_array=()
    read -r -a cov_modes_array <<< "$list_of_cov_modes"

    write_log "defined arrays: max_seq_array=(${max_seq_array[*]}), cov_modes_array=(${cov_modes_array[*]})" "INFO" "$log_file"

    parallel --jobs "${SLURM_CPUS_PER_TASK:-1}" \
             --joblog parallel.joblog \
             run_mmseqs_search_and_convert "$conda_env_prefix" \
                                           "{1}/query_nucl_db/{1/}_nucl_db_type_2" \
                                           "{2}" \
                                           "{1}/results_db" \
                                           "{1}/tmp/cov_{3}_max_{4}" \
                                           3 \
                                           "{3}" \
                                           "{4}" \
                                           ::: "$processing_files_dir"/* \
                                           ::: "$reference_db" \
                                           ::: "${cov_modes_array[@]}" \
                                           ::: "${max_seq_array[@]}" 
}

#combine mmseqs search results per isolate
combine_mmseqs_results() {
    # combines n=6 mmseqs convertalis tsv files into a single file
    # python script removes duplicate entries from the combined results
    #args:
    local conda_env_prefix="$1"
    local processing_files_dir="$2"
    local reference_db="$3"
    local mmseqs_combiner_py_script="$4"

    local reference_db_ncontigs=
    local sample_name=
    local mmseqs_result_dir=
    
    reference_db_ncontigs=$(grep ">" "$reference_db" | wc -l)

    for dir in "$processing_files_dir"/*; do
        sample_name=$(basename "$dir")
        mmseqs_result_dir="$dir/results_db"

        conda run -p "$conda_env_prefix" python3 "$mmseqs_combiner_py_script" \
                                                "$mmseqs_result_dir" \
                                                "$reference_db" \
                                                "$reference_db_ncontigs" \
                                                "$dir/$sample_name"
    done
}

compile_fun() {
    echo
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
write_sample_ID_list "$input_dir" "$output_dir" "$output_dir/logs/run.log"

#removing sequences shorter than 500bp
remove_smalls "$source_directory/subscripts/removesmalls.pl" \
              "$input_dir" \
              "$output_dir/sample_ID_list.txt" \
              "$output_dir/500_bpTrimmed_fastas" \
              "$output_dir/logs/run.log"

#creating mmseqs databases for nucleotide sequences: 1 shared ref and 1 query for each isolate
write_nucl_reference_db "$conda_env_prefix" \
                        "$source_directory/database/$reference_gene_list_name" \
                        "$output_dir/tmp/reference_db" \
                        "$output_dir/logs/run.log" 
write_nucl_query_dbs "$conda_env_prefix" \
                     "$output_dir/500_bpTrimmed_fastas" \
                     "$output_dir/sample_ID_list.txt" \
                     "$output_dir/processing_files" \
                     "$output_dir/logs/run.log"

#initiate search of (max_seq_lengths_array * cov_modes_array) combinations
# exporting functions for GNU parallel
export -f run_mmseqs_search_and_convert write_log
parallel_mmseqs_searches "$conda_env_prefix" \
                        "$output_dir/processing_files" \
                        "$reference_db" \
                        "${max_seq_lengths_array[*]}" \
                        "${cov_modes_array[*]}" \
                        "$output_dir/logs/run.log"

#combine the mmseqs search results per isolate
combine_mmseqs_results "$conda_env_prefix" \
                        "$output_dir/processing_files" \
                        "$reference_db" \
                        "$source_directory/subscripts/mmseq2_results_replicate_combine.py"

