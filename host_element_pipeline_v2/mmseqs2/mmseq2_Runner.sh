#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p project
#SBATCH -e mmseq2_Runner_%A_%a.err
#SBATCH -o mmseq2_Runner_%A_%a.out
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2


# Created by Edward Sung (edward.sung@gwu.edu) on 2/18/2025
#modified high level by Jon Slotved (JOSS@dksund.dk)

STARTTIMER="$(date +%s)"

# User Input
Data_Folder_input=$1
Data_Folder_Samplelist_SLURM_ARRAY_READY_input=$2
Data_Folder_Samplelist_SLURM_ARRAY_READY_index_set_input=$3
Reference_File_input=$4
SearchType_input=$5
Percent_Identity_input=$6
Percent_Coverage_input=$7
main_output_folder_input=$8
#config (defaults to local path)
config_file_local="/dpssi/data/Projects/mtg_host_elements_files_and_output/proj/general_JonThesis/Food-epidemiology/host_element_V3/config/config.env"
config_file=${9:-${config_file_local}}


# Conda Enviroment - Please change to load your conda environment

conda_source=$(grep '^GLOBAL__CONDA_SH__=' "$config_file" | awk -F'__=' '{print $2}' | xargs)
conda_env=$(grep '^HEP__MMSEQ2_ENV__=' "$config_file" | awk -F'__=' '{print $2}' | xargs)
. "$conda_source"
conda activate "$conda_env"

# Script and Tool Locations - Include any additional script path as needed
project_root=$(grep '^GLOBAL__PROJECT_ROOT__=' "$config_file" | awk -F'__=' '{print $2}' | xargs)
mmseq2_python_scripts="$project_root/pipeline_modules/host_element_pipeline/mmseq2/scripts/python_scripts"

threads=20

# Set Search Type
if [[ "$SearchType_input" == "prot" ]]; then
   dbtype=1
   searchType=1

elif [[ "$SearchType_input" == "nucl" ]]; then
   dbtype=2
   searchType=3
else
    echo "Invalid input. Please use 'prot' or 'nucl' for SearchType."
    exit 1
fi


# Temp Flags for Short sequences
#tmp_flags="--mask 0 --prefilter-mode 2 --comp-bias-corr 0"
tmp_flags=""


# Obtains the filename indicated by slurm_array_id
fileInput="$(cat $Data_Folder_Samplelist_SLURM_ARRAY_READY_input | grep "^${Data_Folder_Samplelist_SLURM_ARRAY_READY_index_set_input}__@__${SLURM_ARRAY_TASK_ID}__@__" | awk -F "__@__" '{print $3}')"

reference_file_basename=$(basename "$Reference_File_input")

# Remove extensions to just get filename
filename=${fileInput%.*}
reference_filename=${reference_file_basename%.*}

# Create sample folder for outputs
mkdir $main_output_folder_input/processing_files/$filename
mkdir $main_output_folder_input/processing_files/$filename/slurm_outputs
mkdir $main_output_folder_input/processing_files/$filename/query_db
mkdir $main_output_folder_input/processing_files/$filename/reference_db
mkdir $main_output_folder_input/processing_files/$filename/search_results
mkdir $main_output_folder_input/processing_files/$filename/tmp_space
mkdir $main_output_folder_input/processing_files/$filename/results_replicates


# mmseqs createdb databases of both the query and reference
mmseqs createdb $Data_Folder_input/$fileInput $main_output_folder_input/processing_files/$filename/query_db/${filename}_DB --dbtype $dbtype
mmseqs createdb $Reference_File_input $main_output_folder_input/processing_files/$filename/reference_db/${reference_filename}_DB --dbtype $dbtype




run_mmseqs_replicate() {
   cov_mode_input=$1

   # Loop through the variable max-seq-len
   for val in 10000 35283 65535; do
      mmseqs search \
      $main_output_folder_input/processing_files/$filename/query_db/${filename}_DB \
      $main_output_folder_input/processing_files/$filename/reference_db/${reference_filename}_DB \
      $main_output_folder_input/processing_files/$filename/search_results/${filename}_result_DB_cov${cov_mode_input}_len${val} \
      $main_output_folder_input/processing_files/$filename/tmp_space \
      --search-type $searchType \
      -s 7.5 \
      --min-seq-id $Percent_Identity_input \
      --cov-mode $cov_mode_input \
      -c $Percent_Coverage_input \
      --max-seq-len $val \
      --threads $threads \
      -a \
      --mask 0 \
      --comp-bias-corr 0 \
      --strand 2 \
      --max-seqs 1000
      
 #     --search-type $searchType \               # Search type 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment
 #     -s 7.5 \                                  # Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive
 #     --min-seq-id $Percent_Identity_input \    # List matches above this sequence identity (for clustering) (range 0.0-1.0)
 #     --cov-mode $cov_mode_input \              # 1: coverage of target ; 2: coverage of query
 #     -c $Percent_Coverage_input \              # List matches above this fraction of aligned (covered) residues (see --cov-mode)
 #     --max-seq-len $val \                      # Maximum sequence length
 #     --threads $threads \
 #     -a \                                      # Add backtrace string (convert to alignments with mmseqs convertalis module)
 #     --mask 0 \                                # Mask sequences in prefilter stage with tantan: 0: w/o low complexity masking, 1: with low complexity masking
 #     --comp-bias-corr 0 \                      # Correct for locally biased amino acid composition (range 0-1)
 #     --strand 2                                # Strand selection only works for DNA/DNA search 0: reverse, 1: forward, 2: both
 #     --max-seqs INT		    		 # Maximum results per query sequence allowed to pass the prefilter (affects sensitivity) [300]

      # Convert to readable text
      mmseqs convertalis \
      $main_output_folder_input/processing_files/$filename/query_db/${filename}_DB \
      $main_output_folder_input/processing_files/$filename/reference_db/${reference_filename}_DB \
      $main_output_folder_input/processing_files/$filename/search_results/${filename}_result_DB_cov${cov_mode_input}_len${val} \
      $main_output_folder_input/processing_files/$filename/results_replicates/${filename}_mmseq2_result_cov${cov_mode_input}_len${val}.tsv \
      --format-output query,target,pident,qcov,tcov,alnlen,mismatch,gapopen,qlen,qstart,qend,tlen,tstart,tend,evalue,bits,cigar \
      --search-type $searchType


      # Add header
      sed -i "1i\Query_Seq-id\tSubject_Seq-id\tPercent_Identity\tQuery_Coverage\tSubject_Coverage\tAlignment_Length\tMismatches\tGapOpenings\tQuery_Length\tQuery_Start\tQuery_End\tSubject_Length\tSubject_Start\tSubject_End\tE-Value\tBitscore\tCigar" $main_output_folder_input/processing_files/$filename/results_replicates/${filename}_mmseq2_result_cov${cov_mode_input}_len${val}.tsv

   done
}



# PASS ONE -- cov-mode 1 -- coverage of target, where it screen query against target at coverage -c %
# Covers if query is long and target is short in comparison

run_mmseqs_replicate 1

# PASS TWO -- cov-mode 2 -- coverage of query, where it screen target against query at coverage -c %
# Covers if query is short and target is long in comparison

run_mmseqs_replicate 2



# Python helper function to combine the result replicates and remove any duplicates
ref_file_numsize=$(grep -c "^>" $Reference_File_input)
python $mmseq2_python_scripts/mmseq2_results_replicate_combine.py $main_output_folder_input/processing_files/$filename/results_replicates $Reference_File_input $ref_file_numsize $filename

mv ${filename}_mmseq2_result_compiled.tsv $main_output_folder_input/processing_files/$filename

# This will only generate the presence absence matrix generates, which is currently based on a max size of less than a 1000 seqs.
if [[ "$ref_file_numsize" -lt 5000 ]]; then
   mv ${filename}_mmseq2_result_presence_absence.tsv $main_output_folder_input/processing_files/$filename
fi


# Clean-up File System
rm -r $main_output_folder_input/processing_files/$filename/tmp_space
rm -r $main_output_folder_input/processing_files/$filename/query_db
rm -r $main_output_folder_input/processing_files/$filename/reference_db
rm -r $main_output_folder_input/processing_files/$filename/search_results
mv mmseq2_Runner_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err $main_output_folder_input/processing_files/$filename/slurm_outputs
mv mmseq2_Runner_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $main_output_folder_input/processing_files/$filename/slurm_outputs


# Script Timer
ENDTIMER="$(date +%s)"
DURATION=$[${ENDTIMER} - ${STARTTIMER}]
HOURS=$((${DURATION} / 3600))
MINUTES=$(((${DURATION} % 3600)/ 60))
SECONDS=$(((${DURATION} % 3600) % 60))

echo "RUNTIMER: $HOURS:$MINUTES:$SECONDS (hh:mm:ss)"
