an extension of the mmseqs2 tool: https://github.com/soedinglab/mmseqs2 


## Contacts: 
Jon Sztuk Slotved (JOSS@ssi.dk)
Maliha Aziz (xxx@gwu.edu)



## Conda environment

Create the environment from the included file:

```bash
conda env create -f mmseq2_env.yml
```

## Usage

make a new config, based on template configuration file. Just copy it:

```bash
cp host_element_pipeline_v2/template_config.env host_element_pipeline_v2/config.env
```

set the paths and execution mode in your config.env. Then run the pipeline:

- `-f /path/to/host_file.tsv` is optional. **NOTE:** if used, see below for file struture
- For SLURM, set `execution_mode=slurm` and submit the same script:

```bash
bash host_element_pipeline_v2/run_HEP_analysis.sh \
    -i /path/to/input_fastas \
    -o /path/to/output \
    -c host_element_pipeline_v2/config.env \
    -f /path/to/host_file.tsv
```



## warnings
- usually requires minimum 16GB of system mem.
- for the referece/target fasta database. Avoid more than 1000 genes, as it could overload the system. -Edward Sung


## Pipeline output tree

The pipeline writes this structure under the output directory:
`result_compiled_Main_Data.xlsx` contains summary stats over gene presence
`result_compiled_element_presence.tsv` is the MGE element presence/absence matrix

```text
output/
├── sample_ID_list.txt
├── host_file.tsv
├── 500_bpTrimmed_fastas/
│   └── <sample>.fasta
├── tmp/
│   └── reference_db/
│       └── reference_nucl_db_type_2.*
├── processing_files/
│   └── <sample>/
│       ├── <sample>_mmseq2_result_compiled.tsv
│       └── <sample>_mmseq2_result_presence_absence.tsv
├── compiled_files/
│   ├── mmseq2_result_compiled.tsv
│   ├── mmseq2_result_presence_absence.tsv
│   ├── result_compiled_element_presence.tsv
│   ├── result_compiled_Main_Data.xlsx
│   ├── result_compiled/
│   │   └── <sample>_mmseq2_result_compiled.tsv
│   └── result_presence_absence/
│       └── <sample>_mmseq2_result_presence_absence.tsv
└── logs/
    └── run.log
```

## Some of the output files

### `sample_ID_list.txt`

One filename per line:

```text
Animal5.fasta
Human2.fasta
```

### `host_file.tsv`

The host file is tab-separated. If no host file is supplied, the pipeline creates one using `base_host` from the configuration (Important that header is exactly as below):

```text
Genome_Ref	Host
Animal5	Animal
Human2	Human
```

### `result_compiled_element_presence.tsv`

This is the main result file. It contains element presence or absence calls for each sample:

```text
Genome_Ref	EL01	EL35	EL40 ...
Animal5	1	0	1 ...
Human2	1	1	1 ...
```