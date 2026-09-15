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

    # db_type = 0 for protein, 1 for nucleotide
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

write_nucl_query_dbs() {
    local conda_env_prefix="$1"
    local trimmed_fasta_dir="$2"
    local sample_id_list_file="$3"
    local processing_dir="$4"
    local log_file="$5"
    local sample_id=""
    local database_type=2
    
    # iterates over each sample in the sample list and create a query database for it in processing_files
    while read -r sample_filename; do
        #remove extension and write dir
        sample_id="${sample_filename%.*}"
        mkdir -p "$processing_dir/${sample_id}/query_nucl_db"

        #make db (nucl = --dbtype 2)
        conda run -p "$conda_env_prefix" mmseqs createdb "$trimmed_fasta_dir/$sample_filename" \
                                 "$processing_dir/${sample_id}/query_nucl_db/${sample_id}_nucl_db_type_${database_type}" \
                                 --dbtype $database_type > /dev/null
        
        #log
        local exit_status=$?
        if [[ $exit_status -eq 0 ]]; then
            write_log "wrote query db nucl type_${database_type}: ${sample_id} " "INFO"
        else
            write_log "Failed to write db nucl type_${database_type}: $sample_id" "ERROR" "$log_file"
        fi

    done < "$sample_id_list_file"

    write_log "$(ls "$processing_dir" | wc -l) directories generated in $processing_dir" "INFO" "$log_file"
}

# controlled by parallel_mmseqs_searches function
run_mmseqs_search_and_convert() {
    # summary:
    #   Runs an mmseqs search and converts the results to a human-readable format.
    #   This function is not inteded to be called directly by users.
    #   This function is used internally to orchestrate mmseqs searches and conversions.
    #   The search type is fixed to 3 (nucleotide search) for this pipeline.
    # author: 
    #   Jon Slotved 09-09-2026
    # Arguments:
    #   1: conda_env_prefix <path>
    #   2: query_db         <path>
    #   3: reference_database_prefix <path>
    #   4: results_directory        <path>
    #   5: temporary_directory      <path>
    #   6: coverage_mode            <INT> : 1 (target/ref cov), 2 (query cov).
    #   7: max_sequence_length      <INT> : maximum sequence length for the search.
    #   8: log_file                 <path> : (optional)

    local conda_env_prefix="$1"
    
    #databases
    local query_database_prefix="$2"
    local reference_database_prefix="$3"
    local results_directory="$4"
    local temporary_directory="$5"

    local coverage_mode="$6"
    local max_sequence_length="$7"

    #log
    local log_file="$8"

    #get the sample name
    local query_database_name=""
    query_database_name=$(basename "$query_database_prefix")

    #search and convert outputs
    local search_database_prefix="$results_directory/${query_database_name}_db_cov_${coverage_mode}_max_len_${max_sequence_length}"
    local converted_results_file="$search_database_prefix.tsv"

    # because mmseqs requires the results directory to exist beforehand
    mkdir -p "$results_directory"
    mkdir -p "$temporary_directory"

    #search
    conda run -p "$conda_env_prefix" mmseqs search "$query_database_prefix" \
                                                   "$reference_database_prefix" \
                                                   "$search_database_prefix" \
                                                   "$temporary_directory" \
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

    #logs
    local exit_status=$?
    if [[ $exit_status -eq 0 ]]; then
        write_log "mmseqs search completed: $query_database_name against $(basename "$reference_database_prefix")" "INFO"
    else
        write_log "mmseqs search failed" "ERROR" "$log_file"
        exit $exit_status
    fi

    #convert to human readable
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
        exit $exit_status
    fi
    
    #add header. This is for compatibility with Edwards format (should change later -Jon) 
    #quick fix... Should change later
    cp "$converted_results_file" "$converted_results_file.tmp"
    printf "Query_Seq-id\tSubject_Seq-id\tPercent_Identity\tQuery_Coverage\tSubject_Coverage\tAlignment_Length\tMismatches\tGapOpenings\tQuery_Length\tQuery_Start\tQuery_End\tSubject_Length\tSubject_Start\tSubject_End\tE-Value\tBitscore\tCigar\n" > "$converted_results_file"
    cat "$converted_results_file.tmp" >> "$converted_results_file"
    rm -f "$converted_results_file.tmp"
}

