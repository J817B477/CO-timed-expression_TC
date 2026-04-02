# Project Scaffolding:
#### This project consists of directories pertaining to the steps performed for RNA sequencing data processing differential expression analysis and gene module enrichment
1. [snakefile](snakefile):
    - while not a directory, it is responsible for most of the pipeline execution
    - file that designs and executes the pipeline for the establishing raw expression counts from a ras fastQ files
    - consists of rules that are supported by either bash command bioinformatics tools or custom scripts to improve the  the results of standard rna sequencing pipelines
1. "scripts" directory:
    - contains scripts that were designed to support jobs in the snakemake file
    - more info [here](scripts/README_scriptsFolder.md)
1. "coding_notebooks" directory:
    - folder that is reserved for coding notebooks where more narration and rationalization can articulated for decisions, results and justifications for next steps in the file
    - more info [here](coding_notebooks/ReadMe_notebooks.md)

#### Additional Directories not tracked by this repository
1. bams/
    - contains the files produced from alignment that are utilized to generate counts per gene for each experimental sample
1. hisat_index/
    - contains index files that are used for alignment
    - the memory demands for generating the indices exceeded available resources so the index files were downloaded ans stored here
    - these indices are splice aware allowing RNA fragments that span introns to be aligned effectively to the genes they are expressions of based on neighboring exons for that gene
    - information for getting these files con be found [here](project notes/README_pipelineResources.md)
1. raw_fastqc_reports
    - contains standard folders outputted by fastQC
    - *_data.txt files used to establish overrepresented sequences with "No Hit" values for possible source that are then confirmed to not be technical artifacts using BLAST
1. raw_fastqs
    - primary source data from experiment
    - includes reads of RNA libraries
    - sequencing performed by Plasmidsaurus
    - UMI based deduplication was performed on reads prior to receiving them from Plasimidsaurus
1. trimmed_fastp_reports
    - include json and html formats that share quantity and quality states before and after trimming using fstp
1. trimmed_fastqc_reports
    - same format as raw fastQC reports
    - provide before and after perspective when compared with raw fastQC reports for the same sample
1. trimmed_fastqs
    - updated versions of the original fastq after performing trimming
    - they are the fastq files aligned with the GTF (Gene Transfer Format) annotation file that maps features to reads used for generating expression counts

#### Top level files:
- all support the snakemake pipeline (except gitignore)
1. [.gitignore](.gitignore)
    - manages what project files/folders get tracked and synced to the remote version of the repository on github
1. [config.yml](config.yml)
    - provides parameters to snakefile for source data
    - used to establish 
1. [counts.txt](counts.txt)
    - output of feature counts
    - provides matrix of counts data with samples defining columns and genes defining 
    - main resource that downstream analysis stems from
1. [counts.txt.summary](counts.txt.summary)
    - includes metadata from featureCounts indicating the status of different reads
    - multimapped and assigned reads collectively account for sufficient number of reads
    - fractional accounting for reads mapped to multiple genes was not performed to avoid ambiguity in reporting
1. [custom_adapters.fa](custom_adapters.fa)
    - output of rule that produces fasta file with over-repressented sequences from raw fastq files that are indicated to be non-biological (non-human in this particular experiment)  by BLAST
    - in this particular context none of the sequences checked were established as non-human by BLAST resulting in an empty fasta file
1. [overrepresented_sequences.csv](overrepresented_sequences.csv)
    - compilation of overrepresented sequences from the raw fastqc text reports
1. [snakefile](snakefile)
    - the main file for designing and executing the workflow from raw reads processing through DEGs detected across experimental group contrasts in differential analysis
1. [strand_consensus.txt](strand_consensus.txt)
    - file used to provide the strand direction argument used in the featureCounts call for generating the expression counts matrix (counts.txt)
    - its content is generated programmatically, establishing the overall consensus of read direction across RNA seq libraries that were used to generate the sequencing reads as established by the infer_experiment.py script from RSeQC
1. [updated_samplesNames_map.csv](updated_samplesNames_map.csv)
    - the fastq file's naming protocol was inconsistent and not parsable programmatically for generating the "coldata" metadata utilized in DESeq2's differential analysis
    - this is an artifact of the name update that maps the new name to the old name for each sample to confirm agreement (ie, that each new name updated the correct sample)

This readme file is in the process of being written but the computational aspects of the pipeline are finished. Some of the content is overstated with the purpose of demonstrating understanding and does not necessarily reflect the level of communication that would be exercised when collaborating professionals in the field. I will be adding logging to the clarify the standard output that was informative for analysis of efficacy of pipelines procedures (such as percentage of successful alignment). 