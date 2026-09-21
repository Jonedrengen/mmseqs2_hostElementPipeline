#!/usr/bin/env bash
#SBATCH --job-name=mmseqs_compiler
#SBATCH --output=mmseqs_compiler_%j.out
#SBATCH --error=mmseqs_compiler_%j.err

# Author: Jon Slotved
# Description: Worker script for running SLURM jobs in the host element pipeline

write_log() {
	local log_message="${1:-No log message provided}"
	local log_type="${2:-INFO}"
	local log_file="${3:-}"
	local time=""
	time=$(date +"%Y-%m-%d %H:%M:%S")
	echo "[$log_type] [$time] $log_message"
	if [[ -n "$log_file" ]]; then
		echo "[$log_type] [$time] $log_message" >> "$log_file"
	fi
}

while getopts "p:e:o:h:r:m:l:" opt; do
	case $opt in
		p) pipeline_dir="$OPTARG" ;;
		e) conda_env_prefix="$OPTARG" ;;
		o) output_dir="$OPTARG" ;;
		h) host_file="$OPTARG" ;;
		r) reference_fasta_file="$OPTARG" ;;
		m) module_to_load="$OPTARG" ;;
		l) log_file="$OPTARG" ;;
		*) echo "you should not be passing anything here" ;;
	esac
done

source "$pipeline_dir/subscripts/mmseqs_functionality.sh"
source "$pipeline_dir/subscripts/slurm_functionality.sh"
load_module "$module_to_load" "$log_file"

compile_mmseqs_results \
	"$output_dir/processing_files" \
	"$output_dir/compiled_files" \
	"$log_file"

run_host_element_screen_processor \
	"$conda_env_prefix" \
	"$pipeline_dir/subscripts/host_element_screen_processor.py" \
	"$output_dir/compiled_files/result_compiled" \
	"$host_file" \
	"$reference_fasta_file" \
	"$output_dir/compiled_files/result_compiled"

#moving slurm stuff
mv "$SLURM_SUBMIT_DIR/mmseqs_compiler_${SLURM_JOB_ID}.out" "$output_dir/compiled_files"
mv "$SLURM_SUBMIT_DIR/mmseqs_compiler_${SLURM_JOB_ID}.err" "$output_dir/compiled_files"

write_log "Finished SLURM compiler" "INFO" "$log_file"


