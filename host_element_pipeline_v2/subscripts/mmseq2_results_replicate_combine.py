# mmseq2_results_replicate_combine.py
# Version: 1.2
# Updated: Edward Sung (edward.sung@gwu.edu) on 05/22/25

# Created by: Edward Sung (edward.sung@gwu.edu) on 02/23/25

# updated by: Jon Slotved (JOSS@ssi.dk) on 10/09/2026

# Purpose is to combine the replicates of a mmseq2 (mmseq2_Submitter.sh) by removing duplicates and merging unique hits into one result file.

#%% Import Packages
# Packages
import pandas as pd
import numpy as np
import os
import argparse

'''
The ref_file_numsize variable is necessary for the presence absence matrix. If you are doing a very large reference search, 
ie has more than 1000 fasta sequences, then the presence absence matrix portion will not run or generate.
This is because it will not only take long, but it will generate a large matrix that is not feasible to navigate.
If you still want a presence absence matrix, please utilize the generated mmseq2 table to generate a presence absence matrix.
'''

#%% Argument Parsers
parser = argparse.ArgumentParser()
parser.add_argument("dir_path", help="User provided directory path to the mmseq replicates")
parser.add_argument("ref_file", help="User provided reference file used for mmseq")
parser.add_argument("ref_file_numsize", help="User provided reference file number of fasta in reference")
parser.add_argument("output_name", help="User provided output_file_name")
args = parser.parse_args()

dir_path = args.dir_path + "/"
ref_file = args.ref_file
ref_file_numsize = args.ref_file_numsize
output_name = args.output_name

#%% Testing Purposes
# os.chdir("/Users/edward.sung/Desktop/blast_update_workspace")
#
# dir_path = "/Users/edward.sung/Desktop/dropstyle/colv_mmseq2/testrun_output/processing_files/File1562_genes_genes/results_replicates/"
# ref_file = "/Users/edward.sung/Desktop/dropstyle/colv_mmseq2/colv_db.fasta"
# ref_file_numsize= "999"
# output_name = "ExampleTest"


#%% Import and generate compiled unique mmseq2 results
cols = ["Query_Seq-id",
        "Subject_Seq-id",
        "Percent_Identity",
        "Query_Coverage",
        "Subject_Coverage",
        "Alignment_Length",
        "Mismatches",
        "GapOpenings",
        "Query_Length",
        "Query_Start",
        "Query_End",
        "Subject_Length",
        "Subject_Start",
        "Subject_End",
        "E-Value",
        "Bitscore",
        "Cigar"]

# Initialize combined dataframe
main_data_df = pd.DataFrame(columns=cols)

# Add in each replicate data
for file in os.listdir(dir_path):
    file_path = os.path.join(dir_path, file)

    if file.endswith(".tsv") and os.path.getsize(file_path) > 0:
        rawData = pd.read_csv(file_path, sep="\t", usecols=cols) # Load the data table
        rawData["GenomeName"] = file.split("_")[0] # Attach Genome file name

        main_data_df = pd.concat((main_data_df, rawData), axis=0)

# Drop full row duplicates
main_data_df_unique = main_data_df.drop_duplicates()

# Reset index
main_data_df_unique.reset_index(drop=True, inplace=True)

# Export Table
main_data_df_unique.to_csv(output_name + "_mmseq2_result_compiled.tsv", sep='\t', index=False)
#%% Generate a presence absence matrix row based on what is in reference file.

# Refer to my note at the top of code, but this is to not generate a matrix if there are too many fasta sequences in ref
if int(ref_file_numsize) < 1000:
    reference_fasta_names = []

    # Obtain a list of all fasta name in reference file
    with open(ref_file, "r", encoding="utf-8", errors="replace") as ref:
        for line in ref:
            if line.startswith(">"):  # Identify FASTA headers
                reference_fasta_names.append(line.strip()[1:])  # Remove '>' and store only the name

    # Initialize presence absence matrix
    presence_absence_list = []

    base_dict = {
        "GenomeName": os.path.basename(output_name)
    }

    # Loop through and create presence absence matrix
    for fasta_name in reference_fasta_names:
        if fasta_name in main_data_df_unique["Subject_Seq-id"].values:
            subset_df = main_data_df_unique[main_data_df_unique["Subject_Seq-id"] == fasta_name]

            if subset_df.empty:
                base_dict[fasta_name] = "0"

            else:
                # Filter out for only 80 percent identity and 80 coverage of subject (Reference)
                subset_filter_df = subset_df[(subset_df['Percent_Identity'] >= 80) & (subset_df['Subject_Coverage'] >= 0.8)]

                if subset_filter_df.empty:
                    base_dict[fasta_name] = "0"

                else:
                    subset_filter_df = subset_filter_df.sort_values(by="Bitscore", ascending=False)  # Sort be decending bitscore
                    subset_filter_df.reset_index(drop=True, inplace=True)


                    best_percent_identity = -1
                    best_subject_coverage = -1

                    # But we will loop through to ensure we have the top hit
                    for index in range(len(subset_filter_df)):
                        current_percent_identity = subset_filter_df["Percent_Identity"][index]
                        current_subject_coverage = subset_filter_df["Subject_Coverage"][index]

                        if (current_percent_identity > best_percent_identity and current_subject_coverage > best_subject_coverage):
                            best_percent_identity = current_percent_identity
                            best_subject_coverage = current_subject_coverage

                        elif (current_percent_identity > best_percent_identity and current_subject_coverage == best_subject_coverage):
                            best_percent_identity = current_percent_identity
                            best_subject_coverage = current_subject_coverage

                        elif (current_percent_identity == best_percent_identity) and (current_subject_coverage > best_subject_coverage):
                            best_percent_identity = current_percent_identity
                            best_subject_coverage = current_subject_coverage

                    # Multipled by 100 because the default query coverage is [0-1], while percent_identity is [0-100]
                    best_subject_coverage = best_subject_coverage * 100

                    base_dict[fasta_name] = str(round(best_percent_identity, 2)) + "__" + str(round(best_subject_coverage, 2))
        else:
            base_dict[fasta_name] = "0"

    presence_absence_list.append(base_dict)
    presence_absence_matrix = pd.DataFrame(presence_absence_list)

    presence_absence_matrix.to_csv(output_name + "_mmseq2_result_presence_absence.tsv", sep='\t', index=False)
