# Project Scaffolding:
#### This project consists of directories pertaining to the steps performed for RNA sequencing data processing differential expression analysis and gene module enrichment
1. [snakefile](snakefile):
    - file that designs and executes the pipeline for the establishing raw expression counts from a ras fastQ files
    - consists of rules that are supported by either bash command bioinformatics tools or custom scripts to improve the  the results of standard rna sequencing pipelines
1. "scripts" directory:
    - contains scripts that were designed to support jobs in the snakemake file
    - more info [here](scripts/README_scriptsFolder.md)
1. "coding_notebooks" directory:
    - folder that is reserved for coding notebooks where more narration and rationalization can articulated for decisions, results and justifications for next steps in the file
    - more info [here](coding_notebooks/ReadMe_notebooks_.md)

#### Additional Directories not tracked by this repository
1. bams/
    - contains the files produced from alignment that are utilized to generate counts per gene for each experimental sample
1. hisat_index/
    - contains index files that are used for alignment
    - the memory demands for generating the indices exceeded available resources so the index files were downloaded ans stored here
    - these indices are splice aware allowing RNA fragments that span introns to be aligned effectively to the genes they are expressions of based on neighboring exons for that gene
    - information for getting these files con be found [here](project notes/README_pipelineResources.md)
1. raw_fastqc_reports
1. raw_fastqs
1. trimmed_fastp_reports
1. trimmed_fastqc_reports
1. trimmed_fastqs

#### Top level files:
- all support the snakemake pipeline
1. .gitignore
1. config.yml
1. counts.txt
1. counts.tst.summary
1. custom_adapters.fa
1. overrepresented_sequences.csv
1. snakfile
1. strand_consensus.txt
1. uodated_samplesNames_map.csv

This readme file is in the process of being written but the computational aspects of the pipeline are finished.