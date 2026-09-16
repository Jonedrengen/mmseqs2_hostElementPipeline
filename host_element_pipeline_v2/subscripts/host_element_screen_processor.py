# host_element_screen_processor.py script works as part of host_element_pipeline.sh.
# Version: 4.0
# Updated: Edward Sung (edward.sung@gwu.edu) on 5/27/25

# Processes each of the compiled mmseq2 result_compiled files to generate:

# Summary Files/DataFrames - summary_df, threshold_count_df, presence_df, proportional_df, host_element_prevalence_df, element_presence_df(s).


# Current host_element_screen_processor.py
# Inputs: python host_element_screen_processor.py (path_to_mmseq2_results) (host_file.txt) (element_gene_names.txt) (output_file_names)

#%% Import Packages and Dependencies
# Packages
import pandas as pd
import numpy as np
import xlsxwriter
import os
import argparse


#%% Argument Parsers
parser = argparse.ArgumentParser()
parser.add_argument("dir_path", help="User provided directory path to the mmseq2 compiled results files files")
parser.add_argument("host_file", help="User provided host file for genome inputs")
parser.add_argument("element_gene_names_file", help="Automatically provided file from pipeline for list of element gene names used in mmseq2")
parser.add_argument("output_name", help="User provided output_file_name")
args = parser.parse_args()

dir_path = args.dir_path + "/"
host_file = args.host_file
element_gene_names_file = args.element_gene_names_file
output_name = args.output_name


#%% Testing Purposes
# os.chdir("/Users/edward.sung/Desktop/blast_update_workspace/host_element_processing/")
#
# dir_path = "result_compiled/"
# host_file = "BigFUTI_futimapALL_Updated_HostLabels.txt"
# element_gene_names_file = "20250527_elementgeneList.fasta"
# output_name = "testrun_hostelement"


#%% Locate mmseq2 result files
# List Data Files, searches for _mmseq2_result_compiled.tsv files
fileList = [file for file in os.listdir(dir_path) if file.endswith("_mmseq2_result_compiled.tsv")]


def base_filename_func(file):
    '''
    A function to obtain the basename for different types of files.
    Currently uses if else manually to identify the suffix to remove to obtain basename.
    '''
    suffix_patterns = [
        "__scaffolds-500bpTrimmed_mmseq2_result_compiled.tsv",
        "_scaffolds-500bpTrimmed_mmseq2_result_compiled.tsv",
        "-500bpTrimmed_mmseq2_result_compiled.tsv",
        "_mmseq2_result_compiled.tsv"
    ]
    # If file has a suffix in suffix_patterns, remove it and return (hopefully) basename
    for suffix in suffix_patterns:
        if file.endswith(suffix):
            return file.replace(suffix, "")

    # If not in pattern, just return the filename
    return file


#%% host_info_df DataFrame (For creating host label input file, column names are not important, but column position is important. First column are the Genome_Ref and the second column is Host ID.)
print("Importing Host Labels...")

host_info_df = pd.read_csv(host_file, sep="\t", dtype={'Genome_Ref':str})
host_info_df.columns = ["Genome_Ref", "Host"]

# Remove any suffix from file in fileList to allow for matching between the lists
# Subset down to what is only in fileList
host_info_df = host_info_df[host_info_df['Genome_Ref'].isin([base_filename_func(file) for file in fileList])]

if len(host_info_df) != len(fileList):
    print("Sample(s) are missing a host label.")

HostList = list(host_info_df.Host.unique())

# Adds Poultry as a host (PENDING Expansion to include multiple host/lists)
# Poultry - combination of Chicken and Turkey
poultryAdded = False

# Normalize for case-insensitive comparison
hostlist_lower = [host.lower() for host in HostList]

if "poultry" not in hostlist_lower: # Checks is poultry exists in the HostList
    if all(host in hostlist_lower for host in ["chicken", "turkey"]): # Checks if chicken AND turkey exists in HostList. If do, create Poultry host, which will be used as a new group that is the sum of Chicken and Turkey data
        HostList.append("Poultry")
        poultryAdded = True

print("Host Labels imported.")


#%% # Import Target Element Gene List used for gene screening
print("Importing Target Element Gene List...")

ElementGeneList_cleaned = []

with open(element_gene_names_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            ElementGeneList_cleaned.append(line.strip()[1:])  # remove ">" and any trailing newline

print("Target Element Gene List imported.")


#%% summary_df, threshold_count_df and presence_df
print("Creating summary_df, threshold_count_df and presence_df...")

# summary_df dataframe setup
noData_summary = [] # Tracks 0 hit genomes
summary_dict = {}

# threshold_count_df dataframe setup
threshold_count_dict = {}

# presence_df dataframe setup
presence_dict = {}


# Loop through all the genome files
for file in fileList:
    print("Processing: ", file)
    filepath = dir_path + file

    # Compiles the mmseq2 data into the dict
    Genome_Ref = base_filename_func(file) # Removes suffix to obtain genome reference name

    # Identifies the Host
    Genome_Host = host_info_df.loc[host_info_df["Genome_Ref"] == Genome_Ref, "Host"].values[0]

    # 0 hit mmseq2 result files will have a file size of 0 (Also means all element_genes are not present)
    if os.path.getsize(filepath) == 0:

        # summary_dict data
        summary_dict[Genome_Ref] = [0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    Genome_Host,
                                    "No_Hits"]

        # threshold_count_dict
        gene_hit_count = 0
        called_difference = len(ElementGeneList_cleaned) - gene_hit_count

        threshold_count_dict[Genome_Ref] = [gene_hit_count, called_difference, len(ElementGeneList_cleaned)]

        # presence_dict
        presence_dict[Genome_Ref] = [0 for ElementGene in ElementGeneList_cleaned]

        # Collect all the genomes that are 0 hits mmseq2 files (To be calculated from summary_dict in the future Update)
        noData_summary.append(Genome_Ref)
    else:
        # Load in mmseq2 results data
        rawData_df = pd.read_csv(filepath, sep="\t", header=0, usecols=[0, 1, 2, 4, 14, 15])

        rawData_df["Genome_Ref"] = Genome_Ref
        rawData_df["Host"] = Genome_Host

        # Filter results to only 80% identity and 80% coverage
        filteredData_df = rawData_df[(rawData_df["Percent_Identity"] >= 80) &
                                     (rawData_df["Subject_Coverage"] >= 0.8)]   # Note, mmseq2 returns coverage values at a decimal

        filteredData_gene_hit_list = list(filteredData_df["Subject_Seq-id"].unique())

        # summary_dict data
        summary_dict[Genome_Ref] = [filteredData_df["Percent_Identity"].min(),
                                    filteredData_df["Percent_Identity"].max(),
                                    filteredData_df["Percent_Identity"].mean(),
                                    filteredData_df["Subject_Coverage"].min(),
                                    filteredData_df["Subject_Coverage"].max(),
                                    filteredData_df["Subject_Coverage"].mean(),
                                    filteredData_df["E-Value"].min(),
                                    filteredData_df["E-Value"].max(),
                                    Genome_Host,
                                    "Hits"]

        # threshold_count_dict
        gene_hit_count = len(filteredData_gene_hit_list) # Counts all the element_genes called by mmseq2
        called_difference = len(ElementGeneList_cleaned) - gene_hit_count
        called_proportional = round(gene_hit_count / len(ElementGeneList_cleaned), 4)

        threshold_count_dict[Genome_Ref] = [gene_hit_count, called_difference, called_proportional, len(ElementGeneList_cleaned)]

        # presence_dict
        presence_dict[Genome_Ref] = [1 if ElementGene in filteredData_gene_hit_list else 0 for ElementGene in ElementGeneList_cleaned]



# dict to dataframe function
def dict_to_df_func(dict_input, colnames_input): # colnames_input is a list of column names
    # Function converts dict to dataframe with appropriate column names
    df = pd.DataFrame.from_dict(dict_input, orient = "index")
    df.reset_index(inplace=True)
    df.columns = colnames_input
    return df



# summary_df
summary_df_colnames = ["Genome_Ref", "Percent_Identity_Min", "Percent_Identity_Max", "Percent_Identity_Mean",
                       "Percent_Coverage_Min", "Percent_Coverage_Max", "Percent_Coverage_Mean",
                       "E-Value_Min", "E-Value_Max", "Host", "Query_Hits"]
summary_df = dict_to_df_func(summary_dict, summary_df_colnames)

# summary_df Overview
samples_summary_total = len(fileList)
elementGenes_summary_total = len(ElementGeneList_cleaned)
noData_summary_total = len(noData_summary)
samples_analyzed_summary_total = samples_summary_total - noData_summary_total

summary_overview_df = pd.DataFrame({"Stats": ["Total Number of Genomes",
                                              "Total Number of Genomes with Hits",
                                              "Total Number of Genomes with No Hits",
                                              "Total Number of Element Genes Screened"],
                                    "Count": [samples_summary_total,
                                              samples_analyzed_summary_total,
                                              noData_summary_total,
                                              elementGenes_summary_total]})
print("summary_df created.")



# threshold_count_df
threshold_count_df_colnames = ["Genome_Ref", "Gene_Count_Hits", "Gene_Count_Miss", "Gene_Count_Hits_Prop", "Total_Element_Genes"]
threshold_count_df = dict_to_df_func(threshold_count_dict, threshold_count_df_colnames)
print("threshold_count_df created.")



# presence_df
presence_df_colnames = ElementGeneList_cleaned.copy() # Copy ElementGeneList to insert Element_Gene name into the list for column names
presence_df_colnames.insert(0, "Genome_Ref")
presence_df = dict_to_df_func(presence_dict, presence_df_colnames)

print("presence_df created.")


#%% proportional_df and host_element_prevalence_df based on Host Groups (Derived from presence_df)
print("Creating proportional_df and host_element_prevalence_df...")

# Initalize proportional_df
proportional_df = pd.DataFrame({"Element_Gene" : ElementGeneList_cleaned})

# Initalize host_element_prevalence_df
host_element_prevalence_df = pd.DataFrame(pd.Series(ElementGeneList_cleaned).str.split('_').str[0].value_counts()).reset_index()
host_element_prevalence_df.columns = ["Element", "Total_Num_Genes_in_Element"]
# Sort by count (descending) then by element name (ascending)
host_element_prevalence_df.sort_values(by=["Total_Num_Genes_in_Element", "Element"],
                                       ascending=[False, True],
                                       inplace=True)

for host in HostList:
    if host == "Poultry" and poultryAdded: # This is true if we add the Poultry host, else if Poultry already existed, will be treated as a stand-alone host like the other hosts
        # Combines all the Chicken and Turkey genomes as one list
        chicken_ref_genomes = list(host_info_df[host_info_df.Host == "Chicken"].Genome_Ref)
        turkey_ref_genomes = list(host_info_df[host_info_df.Host == "Turkey"].Genome_Ref)
        host_ref_genomes = chicken_ref_genomes + turkey_ref_genomes
    else:
        host_ref_genomes = list(host_info_df[host_info_df.Host == host].Genome_Ref)

    # Totals up the number of genomes (pior to adding Element_Gene column)
    total_host_ref_genomes = len(host_ref_genomes)

    # Subset the presence_df for only host specific genomes
    host_presence_df = presence_df[presence_df["Genome_Ref"].isin(host_ref_genomes)]

    host_prevalence_df = host_presence_df.copy().T
    host_prevalence_df.columns = host_prevalence_df.iloc[0]
    host_prevalence_df = host_prevalence_df[1:].reset_index()
    host_prevalence_df.rename(columns={"index": "Element_Gene"}, inplace=True)

    # proportional_df calculations
    # Column Names
    host_Total_Num_Genome_Refs = host + "__Total_Num_Genomes"
    host_Num_Genome_Refs_Called = host + "__Num_Genomes_Called"
    host_Proportional_Called = host + "__Proportional_Called"

    # Generate host specific prop_df (this is to be merged with the main proportional_df)
    host_prop_df = pd.DataFrame(host_presence_df.sum())
    host_prop_df = host_prop_df.iloc[1:].reset_index() # Drop the first row, which is the subset genome names

    # Rename columns
    host_prop_df.columns = ["Element_Gene", host_Num_Genome_Refs_Called]

    # Insert totals genomes column
    host_prop_df[host_Total_Num_Genome_Refs] = total_host_ref_genomes

    # Calculate the proportional abundance per element gene
    host_prop_df[host_Proportional_Called] = host_prop_df[host_Num_Genome_Refs_Called] / total_host_ref_genomes
    host_prop_df[host_Proportional_Called] = host_prop_df[host_Proportional_Called].astype(float)
    host_prop_df[host_Proportional_Called] = host_prop_df[host_Proportional_Called].round(4)

    # Column re-organize
    host_prop_df = host_prop_df[["Element_Gene", host_Total_Num_Genome_Refs, host_Num_Genome_Refs_Called, host_Proportional_Called]]

    # Merge with the main proportional_df
    proportional_df = pd.merge(proportional_df, host_prop_df, on="Element_Gene", how="left")



    # element_prevalence_df calculations
    host_prevalence_df["Element"] = host_prevalence_df["Element_Gene"].str.split("_").str[0] # Extract Element
    host_prevalence_df.drop("Element_Gene", axis=1, inplace=True)
    host_prevalence_df_sum = host_prevalence_df.groupby("Element").sum() # Sum up all genes per element

    # Convert number to be 1 or 0 for called Element, as long as there is a called gene within
    host_prevalence_df_called = host_prevalence_df_sum.where(~(host_prevalence_df_sum > 0), 1)

    host_prevalence_df_called[host_Num_Genome_Refs_Called] = host_prevalence_df_called.sum(axis=1) # Sum across the row (Element)
    host_prevalence_df_called.reset_index(inplace=True)

    host_prevalence_df_called[host_Total_Num_Genome_Refs] = total_host_ref_genomes

    host_prevalence_df_called[host_Proportional_Called] = host_prevalence_df_called[host_Num_Genome_Refs_Called] / total_host_ref_genomes
    host_prevalence_df_called[host_Proportional_Called] = host_prevalence_df_called[host_Proportional_Called].astype(float)
    host_prevalence_df_called[host_Proportional_Called] = host_prevalence_df_called[host_Proportional_Called].round(4)

    host_prevalence_df_called_ToMerge = host_prevalence_df_called[["Element", host_Total_Num_Genome_Refs, host_Num_Genome_Refs_Called, host_Proportional_Called]]

    host_element_prevalence_df = pd.merge(host_element_prevalence_df, host_prevalence_df_called_ToMerge, on="Element", how="left")


print("proportional_df created.")
print("host_element_prevalence_df created.")


#%% element_presence_df DataFrame (Derived from presence_df)
print("Creating element_presence_df...")

# Initialize
element_presence_dict = {"Genome_Ref" : presence_df["Genome_Ref"].tolist()}
element_gene_dict = {}

# Exception Rule to Element Calling; Requires at least 2 genes called to call element present
exception_element_calling = ["EL18", "EL35", "EL37", "EL38", "EL40"]

# Fill in the dict for each element, all the associated genes
for gene in ElementGeneList_cleaned:
    ele = gene.split('_')[0] # Extract the element

    if ele in element_gene_dict:
        element_gene_dict[ele].append(gene)
    else:
        element_gene_dict[ele] = [gene]

# Subset presence_df and populate element_presence_dict
for ele, gene_list in element_gene_dict.items():
    subset_list = ["Genome_Ref"] + gene_list

    # Subset presence_df to only selected genes for an element
    subset_presence_df = presence_df[subset_list].copy()

    # sum (row) all called genes for each genome
    subset_presence_df["Genes_Called"] = subset_presence_df.iloc[:, 1:].sum(axis=1)

    # Make element calls for each genome
    if ele in exception_element_calling:
        element_presence_dict[ele] = [1 if val >= 2 else 0 for val in subset_presence_df["Genes_Called"]]
    else:
        element_presence_dict[ele] = [1 if val >= 1 else 0 for val in subset_presence_df["Genes_Called"]]

# Convert dict to dataframe
element_presence_df = pd.DataFrame(element_presence_dict)

print("element_presence_df created.")


#%% Exporting DataFrames
print("Exporting data...")

# Main_Data Excel
file_output_main = output_name + "_Main_Data.xlsx"
wb = pd.ExcelWriter(file_output_main, engine='xlsxwriter')

summary_df.to_excel(wb, sheet_name="Summary_Data")
summary_overview_df.to_excel(wb, sheet_name="Summary_Data", index=False, startcol=13)
threshold_count_df.to_excel(wb, sheet_name='Threshold_Count_Data')
presence_df.to_excel(wb, sheet_name='Presence_Data')
proportional_df.to_excel(wb, sheet_name='Proportional_Data')
host_element_prevalence_df.to_excel(wb, sheet_name='Host_Element_Prevalence')

#JON_IMPLEMENTED_______________________________
#wb.save()
wb.close() # save is deprecated in newer versions
#JON_IMPLEMENTED_______________________________

print("Main_Data exported.")

# Element_Presence DataFrames
element_presence_output_name = output_name + "_element_presence.tsv"
element_presence_df.to_csv(element_presence_output_name, sep='\t', index=False)
print(element_presence_output_name + " exported.")

print("All Data exported.")
print("host_element_screen_processor.py script complete.")
