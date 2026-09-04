#!/bin/bash

# This Slurm_Array_SampleListReady.sh is a helper script for Edward's (edward.sung@gwu.edu) slurm array scripts.
# It adds the indexing numbers and "__@__" separator/identifier in front of each sample name in the sample list.
# The slurm array scripts uses these indexing numbers in association with the slurm task ID call to parses out the filenames.

# Created by Edward Sung (edward.sung@gwu.edu) on 2/25/2024

# User Inputs
samplelist_input=$1

# Adds two indexes (first number tracks set per 1000, second number tracks the sample in the set (1-1000) to the start of each row for the sample names
samplelist_filename=${samplelist_input%.*}

# This Surm_MaxArraySize to be adjusted based on the HPC system's MaxArraySize allowed. For Pegasus, this is 1000
Slurm_MaxArraySize=1000

# i is the set index for groups of 1000s
i=0

# j is the sample index ranging from 1 up to 1000 per set index
j=1

# Initialize SLURM-ARRAY-READY.txt file
touch ${samplelist_filename}_SLURM-ARRAY-READY.txt

while read -r line
do
   # Appends my indexes to the filenames in the provided samplelist
   str="${i}__@__${j}__@__${line}"
   echo $str >> ${samplelist_filename}_SLURM-ARRAY-READY.txt
 
   if [ $j -lt $Slurm_MaxArraySize ]
   then
      ((j+=1)) # increments j by 1 until 1000
   else
      ((i+=1)) # increments set index by 1
      j=1 # resets the j indexing back to 1 after hitting 1000
   fi
done < $samplelist_input

echo "Indexing numbers and __@__ has been appended to each line in the sample list."
echo "${samplelist_filename}_SLURM-ARRAY-READY.txt is ready, thank you!"
