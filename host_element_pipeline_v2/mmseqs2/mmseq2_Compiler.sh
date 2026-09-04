#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p project
#SBATCH -o mmseq2_Compiler_%j.out
#SBATCH -e mmseq2_Compiler_%j.err

# -p tiny,small-gpu,debug

#modified high level by Jon Slotved (JOSS@dksund.dk)


STARTTIMER="$(date +%s)"

# Inputs
main_output_folder_input=$1
reference_input=$2
jobname_input=$3
#config (defaults to local path)
config_file_local="/dpssi/data/Projects/mtg_host_elements_files_and_output/proj/general_JonThesis/Food-epidemiology/host_element_V3/config/config.env"
config_file=${4:-${config_file_local}}


# Conda Enviroment - Please change to load your conda environment

conda_source=$(grep '^GLOBAL__CONDA_SH__=' "$config_file" | awk -F'__=' '{print $2}' | xargs)
conda_env=$(grep '^HEP__CONDA_ENV__=' "$config_file" | awk -F'__=' '{print $2}' | xargs)
. "$conda_source"
conda activate "$conda_env"




# Create file system
mkdir $main_output_folder_input/compiled_files
mkdir $main_output_folder_input/compiled_files/result_compiled
mkdir $main_output_folder_input/compiled_files/result_presence_absence
mkdir $main_output_folder_input/compiled_files/slurmFiles


# Create files to store youre results
echo -e "Query_Seq-id\tSubject_Seq-id\tPercent_Identity\tQuery_Coverage\tSubject_Coverage\tAlignment_Length\tMismatches\tGapOpenings\tQuery_Length\tQuery_Start\tQuery_End\tSubject_Length\tSubject_Start\tSubject_End\tE-Value\tBitscore\tCigar\tGenomeName" > $main_output_folder_input/compiled_files/${jobname_input}_mmseq2_result_compiled.tsv

echo -e "GenomeName\t$(grep "^>" $reference_input | sed 's/>//' | paste -s -d '\t')" > $main_output_folder_input/compiled_files/${jobname_input}_mmseq2_result_presence_absence.tsv


# Create a list of your output folders
ls $main_output_folder_input/processing_files > tmplist_output_folders


# Loop through your output folders to compile your results
while read -r line
do
   echo "Extracting results from: $line"
   tail -n +2 $main_output_folder_input/processing_files/$line/${line}_mmseq2_result_compiled.tsv >> $main_output_folder_input/compiled_files/${jobname_input}_mmseq2_result_compiled.tsv

   tail -n +2 $main_output_folder_input/processing_files/$line/${line}_mmseq2_result_presence_absence.tsv >> $main_output_folder_input/compiled_files/${jobname_input}_mmseq2_result_presence_absence.tsv

   cp $main_output_folder_input/processing_files/$line/${line}_mmseq2_result_compiled.tsv $main_output_folder_input/compiled_files/result_compiled
   cp $main_output_folder_input/processing_files/$line/${line}_mmseq2_result_presence_absence.tsv $main_output_folder_input/compiled_files/result_presence_absence

done < tmplist_output_folders


# Clean-up file system
rm tmplist_output_folders
mv mmseq2_Compiler_${SLURM_JOB_ID}.out $main_output_folder_input/compiled_files/slurmFiles
mv mmseq2_Compiler_${SLURM_JOB_ID}.err $main_output_folder_input/compiled_files/slurmFiles


# Script Timer
ENDTIMER="$(date +%s)"
DURATION=$[${ENDTIMER} - ${STARTTIMER}]
HOURS=$((${DURATION} / 3600))
MINUTES=$(((${DURATION} % 3600)/ 60))
SECONDS=$(((${DURATION} % 3600) % 60))

echo "RUNTIMER: $HOURS:$MINUTES:$SECONDS (hh:mm:ss)"
