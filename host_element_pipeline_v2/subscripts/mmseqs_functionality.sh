#!/usr/bin/env bash

# Author: Jon Slotved
# Description: Script containing MMseqs2 functionality db creation, search and conversion for the host element pipeline
# why? - because of slurm parallelization optimization (and simplicity)

write_nucl_reference_db() {
    local conda_env_prefix="$1"
    local reference_fasta_file="$2"
    local reference_db_dir="$3"
    local log_file="$4"
    local database_type=2
    # ref db needed globally
    reference_db_prefix=""


    mkdir -p "$reference_db_dir"
    reference_db_prefix="$reference_db_dir/reference_nucl_db_type_${database_type}"
    write_log "reference_db_prefix=$reference_db_prefix" "INFO" "$log_file"

    # db_type = 0 for auto, 2 for nucleotide
    conda run -p "$conda_env_prefix" mmseqs createdb "$reference_fasta_file" \
                                                     "$reference_db_prefix" \
                                                     --dbtype $database_type > /dev/null
    
    local exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "wrote ref db {dbtype=nucl}: $reference_db_dir" "INFO" "$log_file"
    else
        write_log "failed to write ref db {dbtype=nucl}: $reference_db_dir" "ERROR" "$log_file"
    fi
}

write_nucl_query_db() {
    local conda_env_prefix="$1"
    local trimmed_fasta="$2"
    local query_database_prefix="$3"
    local log_file="${4:-}"
    local database_type=2

    mkdir -p "$(dirname "$query_database_prefix")"

    conda run -p "$conda_env_prefix" mmseqs createdb "$trimmed_fasta" \
                                                     "$query_database_prefix" \
                                                     --dbtype "$database_type" > /dev/null

    #log the exit status of the database creation
    local exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "wrote query db nucl type_${database_type}: $(basename "$query_database_prefix")" "INFO" "$log_file"
    else
        write_log "Failed to write query db nucl type_${database_type}: $(basename "$query_database_prefix")" "ERROR" "$log_file"
        return 1
    fi
}



# Runs all coverage-mode and maximum-length combinations for one isolate.
run_mmseqs_search_and_convert() {
    # Arguments:
    #   1: conda_env_prefix        <path>
    #   2: query_database_prefix   <path>
    #   3: reference_database_prefix <path>
    #   4: results_directory       <path>
    #   5: temporary_directory     <path>
    #   6: coverage_modes          <space-separated INTs>
    #   7: max_sequence_lengths    <space-separated INTs>
    #   8: log_file                <path> (optional)

    local conda_env_prefix="$1"
    local query_database_prefix="$2"
    local reference_database_prefix="$3"
    local results_directory="$4"
    local temporary_directory="$5"
    local coverage_modes="$6"
    local max_sequence_lengths="$7"
    local log_file="${8:-}"

    local query_database_name
    query_database_name=$(basename "$query_database_prefix")

    read -r -a coverage_modes_array <<< "$coverage_modes"
    read -r -a max_sequence_lengths_array <<< "$max_sequence_lengths"
    for coverage_mode in "${coverage_modes_array[@]}"; do
        for max_sequence_length in "${max_sequence_lengths_array[@]}"; do
            local search_database_prefix="$results_directory/${query_database_name}_db_cov_${coverage_mode}_max_len_${max_sequence_length}"
            local converted_results_file="$search_database_prefix.tsv"
            local search_temporary_directory="$temporary_directory/cov_${coverage_mode}_max_${max_sequence_length}"

            mkdir -p "$results_directory" "$search_temporary_directory"

            conda run -p "$conda_env_prefix" mmseqs search "$query_database_prefix" \
                                                           "$reference_database_prefix" \
                                                           "$search_database_prefix" \
                                                           "$search_temporary_directory" \
                                                           --search-type 3 \
                                                           --cov-mode "$coverage_mode" \
                                                           --max-seq-len "$max_sequence_length" \
                                                           --min-seq-id 0.8 \
                                                           -c 0.8 \
                                                           -s 7.5 \
                                                           --max-seqs 1000 \
                                                           --threads 1 \
                                                           -a \
                                                           --mask 0 \
                                                           --comp-bias-corr 0 \
                                                           --strand 2 > /dev/null

            local exit_status=$?
            if [[ $exit_status -eq 0 ]]; then
                write_log "mmseqs search completed: $query_database_name against $(basename "$reference_database_prefix")" "INFO" "$log_file"
            else
                write_log "mmseqs search failed" "ERROR" "$log_file"
                return "$exit_status"
            fi

            conda run -p "$conda_env_prefix" mmseqs convertalis "$query_database_prefix" \
                                                               "$reference_database_prefix" \
                                                               "$search_database_prefix" \
                                                               "$converted_results_file" \
                                                               --format-output query,target,pident,qcov,tcov,alnlen,mismatch,gapopen,qlen,qstart,qend,tlen,tstart,tend,evalue,bits,cigar \
                                                               --search-type 3 \
                                                               --threads 1 > /dev/null

            exit_status=$?
            if [[ $exit_status -eq 0 ]]; then
                write_log "mmseqs convertalis completed: $query_database_name/$(basename "$converted_results_file")" "INFO" "$log_file"
            else
                write_log "mmseqs convertalis failed" "ERROR" "$log_file"
                return "$exit_status"
            fi

            cp "$converted_results_file" "$converted_results_file.tmp"
            printf "Query_Seq-id\tSubject_Seq-id\tPercent_Identity\tQuery_Coverage\tSubject_Coverage\tAlignment_Length\tMismatches\tGapOpenings\tQuery_Length\tQuery_Start\tQuery_End\tSubject_Length\tSubject_Start\tSubject_End\tE-Value\tBitscore\tCigar\n" > "$converted_results_file"
            cat "$converted_results_file.tmp" >> "$converted_results_file"
            rm -f "$converted_results_file.tmp"        
        done
    done
}