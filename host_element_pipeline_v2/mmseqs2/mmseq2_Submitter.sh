#!/bin/sh

# Created by Edward Sung (edward.sung@gwu.edu) on 02/18/2025
# Updated by Edward - 5/22/25
# Version v1.2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -p project

#modified high level by Jon Slotved (JOSS@dksund.dk)

# User Inputs
Data_Folder_input=$1
Data_Folder_Samplelist_input=$2
Reference_File_input=$3
SearchType_input=$4		# {'prot' or 'nucl'}
Percent_Identity_input=$5	# {0.8}
Percent_Coverage_input=$6	# {0.8}
Job_Name_input=$7
partition=${8:-project}
dependency=${9:-}
#config (defaults to local path)
config_file_local="/dpssi/data/Projects/mtg_host_elements_files_and_output/proj/general_JonThesis/Food-epidemiology/host_element_V3/config/config.env"
config_file=${10:-${config_file_local}}

if [ -z "$config_file" ] || [ ! -f "$config_file" ]; then
   echo "ERROR: Config file missing or not provided: $config_file"
   exit 1
fi

# Script Locations (Path to where all slurm-array scripts live, use `pwd` to find path.
#HEP_mmseq2_scripts="/scratch/liu_price_lab/ehsung/github/Development/ehsung/microbiome/mmseq2/scripts"
project_root=$(grep '^GLOBAL__PROJECT_ROOT__=' "$config_file" | awk -F'__=' '{print $2}' | xargs)
HEP_mmseq2_scripts="$project_root/pipeline_modules/host_element_pipeline/mmseq2/scripts"


dep_flag=""
[ -n "$dependency" ] && dep_flag="--dependency=afterok:${dependency}"

# Error for required number of inputs
if [ $# -lt 7 ]
then
   echo "Please give all 6 arugments: (Data_Folder_input) (Data_Folder_Samplelist_input) (Reference_File_input) (SearchType_input ['prot' or 'nucl']) (Percent_Identity_input [0 - 1.0]) (Percent_Coverage_input [0 - 1.0]) (Job_Name_input)"
   exit 1
fi


# Start-Up / Generate SLURM-ARRAY-READY samplelist
# Convert samplelist to SLURM-ARRAY-READY format: appends indexing and __@__ to each sample name in file
$HEP_mmseq2_scripts/Slurm_Array_SampleListReady.sh $Data_Folder_Samplelist_input
samplelist_filename=${Data_Folder_Samplelist_input%.*} # Strip extensions
echo 
echo


# Create File System
mkdir ${Job_Name_input}_output
mkdir ${Job_Name_input}_output/processing_files



# Version Tracker File
touch ${Job_Name_input}_output/version_tracker.txt
echo "mmseq2 : v17.b804f" >> ${Job_Name_input}_output/version_tracker.txt



# Slurm Array Settings
numFiles=$(cat ${samplelist_filename}_SLURM-ARRAY-READY.txt | wc -l) # Gets the total number of samples/files
Slurm_MaxArraySize=1000 # Maximum HPC MaxArraySize

# Calculates the SlurmJobArray parallel submissions
if (( $numFiles % $Slurm_MaxArraySize == 0 ))
then
   Slurm_chunks=`expr $numFiles / $Slurm_MaxArraySize`
else
   Slurm_chunks=`expr $numFiles / $Slurm_MaxArraySize + 1` # This is a ceiling int calculation
fi

# Most used nodes 12, minimum is 1, it will only use 1 if its submits more than 6 batches
# Please adjust these numbers accordingly to your specifications or HPC needs
# Example: if samplelist contains 1000 files, it will submit 1 SlurmArray job that will run use 12 compute nodes at a time.
# Example: if samplelist contains 10000 files, it will submit 10 SlurmArray jobs that will run only 1 compute node per SlurmArray job at a time.
if [ $Slurm_chunks == 1 ]
then
   Slurm_CalcRunParallel=12

elif [ $Slurm_chunks == 2 ]
then
   Slurm_CalcRunParallel=6

elif [ $Slurm_chunks == 3 ]
then
   Slurm_CalcRunParallel=4

elif [ $Slurm_chunks == 4 ]
then
   Slurm_CalcRunParallel=3

elif [ $Slurm_chunks == 5 ]
then
   Slurm_CalcRunParallel=2

elif [ $Slurm_chunks == 6 ]
then
   Slurm_CalcRunParallel=2

else
   Slurm_CalcRunParallel=1
fi

# Manually Set Slurm_CalcRunParallel for this special case (to speed up process)
# If you need to speed-up the process, you can manually select the number of nodes to be used per SlurmArray job.
# Uncomment (delete the #) and type in your int number.
#Slurm_CalcRunParallel=2

# Splits up the samplelist by index_set and runs the SlurmArray jobs based on the start and end of each set
for ((i=0; i<$Slurm_chunks; i++))
do
   # Tracks the index_set number for the samplelist (First number)
   index_set=$i

   # Identifes the start and end of the array to be submitted
   array_start="$(cat ${samplelist_filename}_SLURM-ARRAY-READY.txt | grep "^${i}__" | head -1 | awk -F "__@__" '{print $2}')"
   array_end="$(cat ${samplelist_filename}_SLURM-ARRAY-READY.txt | grep "^${i}__" | tail -1 | awk -F "__@__" '{print $2}')"
   
   # Submit the jobs to HPC
   echo "sbatch -p $partition --array=${array_start}-${array_end}%${Slurm_CalcRunParallel} -J $Job_Name_input $HEP_mmseq2_scripts/mmseq2_Runner.sh $Data_Folder_input ${samplelist_filename}_SLURM-ARRAY-READY.txt $index_set $Reference_File_input $SearchType_input $Percent_Identity_input $Percent_Coverage_input ${Job_Name_input}_output $config_file"
   sbatch -p $partition $dep_flag --array=$array_start-$array_end%$Slurm_CalcRunParallel -J $Job_Name_input $HEP_mmseq2_scripts/mmseq2_Runner.sh $Data_Folder_input ${samplelist_filename}_SLURM-ARRAY-READY.txt $index_set $Reference_File_input $SearchType_input $Percent_Identity_input $Percent_Coverage_input ${Job_Name_input}_output "$config_file"
done

# Compile the results data and Clean-up file system script
sbatch -p $partition --dependency=singleton -J $Job_Name_input $HEP_mmseq2_scripts/mmseq2_Compiler.sh ${Job_Name_input}_output $Reference_File_input $Job_Name_input "$config_file"

echo "---------- Your jobs have been submitted to HPC, thank you. ----------"
