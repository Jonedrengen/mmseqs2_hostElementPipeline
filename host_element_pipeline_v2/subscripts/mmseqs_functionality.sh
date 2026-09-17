#!/usr/bin/env bash

# Author: Jon Slotved
# Description: Script containing MMseqs2 functionality db creation, search and conversion for the host element pipeline
# why? - because of slurm parallelization optimization (and simplicity)

# writes a nucleotide reference database from the given reference FASTA file
write_nucl_reference_db() {
    local conda_env_prefix="$1"
    local reference_fasta_file="$2"
    local reference_db_dir="$3"
    local log_file="${4:-}"
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
#writes 1 nucleotide query database for a given isolate
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
        write_log "wrote query db nucl type_${database_type}: $(basename "$query_database_prefix")" "INFO"
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

    write_log "Starting mmseqs search and conversion for $query_database_name" "INFO"

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
                write_log "mmseqs search completed: $query_database_name against $(basename "$reference_database_prefix")" "INFO"
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
                write_log "mmseqs convertalis completed: $query_database_name/$(basename "$converted_results_file")" "INFO"
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

#for local mode. Runs GNU parallel to execute mmseqs searches and conversions in parallel
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
    local log_file="${6:-}"

    # strings to arrays, because GNU parallel need arrays for ::: expansion
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

#combine mmseqs search results per isolate (and cleans up temporary files)
combine_mmseqs_results_per_isolate() {
    # Combines one isolate's mmseqs convertalis files into a single file.
    # Arguments:
    #   1: conda_env_prefix
    #   2: isoalte_results_dir
    #   3: reference_fasta_file
    #   4: results_combiner_script_file
    #   5: log_file (optional)
    local conda_env_prefix="$1"
    local isoalte_results_dir="$2"
    local reference_fasta_file="$3"
    local results_combiner_script_file="$4"
    local log_file="${5:-}"


    local reference_sequence_count
    local sample_dir
    local sample_id

    reference_sequence_count=$(grep -c "^>" "$reference_fasta_file")
    sample_dir=$(dirname "$isoalte_results_dir")
    sample_id=$(basename "$sample_dir")

    conda run -p "$conda_env_prefix" python3 "$results_combiner_script_file" \
                                            "$isoalte_results_dir" \
                                            "$reference_fasta_file" \
                                            "$reference_sequence_count" \
                                            "$sample_dir/$sample_id"
    local exit_status=$?
    if [[ $exit_status -ne 0 ]]; then
        write_log "Failed to combine mmseqs results for $sample_id" "ERROR" "$log_file"
    else
        write_log "Successfully combined mmseqs results for $sample_id" "INFO"
    fi
}

# combine all per-isolate result files into one - Edwards compiler method
compile_mmseqs_results() {
    # Arguments:
    #   1: processing_files_dir
    #   2: compiled_files_dir
    #   3: log_file (optional)
    local processing_dir="$1"
    local compiled_dir="$2"
    local log_file="${3:-}"

    local result_compiled_dir="$compiled_dir/result_compiled"
    local result_presence_dir="$compiled_dir/result_presence_absence"
    local combined_results="$compiled_dir/mmseq2_result_compiled.tsv"
    local combined_presence="$compiled_dir/mmseq2_result_presence_absence.tsv"

    mkdir -p "$result_compiled_dir" "$result_presence_dir"
    : > "$combined_results"
    : > "$combined_presence"

    local sample_dir
    local sample_id
    local first_sample=true
    for sample_dir in "$processing_dir"/*; do
        sample_id=$(basename "$sample_dir")

        #add header (wish I thought of this stategy, but AI recommended it)
        if $first_sample; then
            head -n 1 "$sample_dir/${sample_id}_mmseq2_result_compiled.tsv" > "$combined_results"
            head -n 1 "$sample_dir/${sample_id}_mmseq2_result_presence_absence.tsv" > "$combined_presence"
            first_sample=false
        fi

        tail -n +2 "$sample_dir/${sample_id}_mmseq2_result_compiled.tsv" >> "$combined_results"
        tail -n +2 "$sample_dir/${sample_id}_mmseq2_result_presence_absence.tsv" >> "$combined_presence"
        cp "$sample_dir/${sample_id}_mmseq2_result_compiled.tsv" "$result_compiled_dir/"
        cp "$sample_dir/${sample_id}_mmseq2_result_presence_absence.tsv" "$result_presence_dir/"
    done
    write_log "$(ls "$result_compiled_dir" | wc -l) mmseqs result files compiled successfully" "INFO" "$log_file"
    write_log "$(ls "$result_presence_dir" | wc -l) mmseqs presence/absence result files compiled successfully" "INFO" "$log_file"
}

run_host_element_screen_processor() {
    local conda_env_prefix="$1"
    local host_element_screen_processor_script="$2"
    local compiled_dir="$3"
    local host_file="$4"
    local fasta_gene_file="$5"
    local output_dir="$6"

    "$conda_env_prefix/bin/python" "$host_element_screen_processor_script" "$compiled_dir" \
                                                                         "$host_file" \
                                                                         "$fasta_gene_file" \
                                                                         "$output_dir"
    
    local exit_status=$?
    if [[ $exit_status -ne 0 ]]; then
        write_log "Failed to run host element screen processor" "ERROR" "$log_file"
    else
        write_log "Successfully ran host element screen processor" "INFO" "$log_file"
    fi
}