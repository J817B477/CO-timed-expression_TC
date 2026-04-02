snakefile directives resources
- provides information for all scripts and programs used to support rules/jobs of [snakemake](../snakefile) pipeline

1. raw reads quality check: FastQC/MultiQC
	- supports `raw_qc`, `raw_multiqc`, `trimmed_qc`, `trimmed_multiqc` rules 
1. [scripts/generate_overrep_fa.py](../scripts/generate_overrep_fa.py)
	- supports 'get_overrep_fa' Rule
	- output custom fasta file to use in trimming with fastp
	- is empty if all overrepresented sequences are biological
1. trimming: fastp
	- uses custom_adapters.fa as --adapter_fa argument
1. alignment: hisat2
	- supports align_reads rule 
	- relies on pre-produced splice aware indexes downloaded as stated here
	```bash
    wget --content-disposition ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz
	```
1. strand direction inference: infer_experiment.py (RSeQC)
	- supports infer_strand rule
	- significantly improves FeatureCounts results
	- outputs inference results to files of strand\ folder
1. [scripts/infer_strand.py](../scripts/infer_strand.py)
	- supports summarize_strand rule
	- generates strand direction consensus
	- output file contents is argument value for eventual featureCounts call's "-s" flag
	- potential values: "1" means forward, "2" means reverse, "0" means unstranded
1. counts generation: featureCounts
	- supports get_counts rule 
	- uses strand_consensus.txt for clarifying strandedness
1. [scripts/clean_counts.R](../scripts/clean_counts.R)
	- supports clean_counts rule
	- updates columns names to parsable experimental design encoded strings
	- exports column renaming map for proper renaming validation
	- appends transcripts per million standardized counts per sample for latter reporting purposes
1. [scripts/DEanalysis_expandedContrasts.R](../scripts/DEanalysis_expandedContrasts.R)
	- supports 
	- generates DEG datasets for secondary analysis on jupyter notebook
		- includes adjusted p-values, $log_{2}(fc)$ values, and perturbation direction in independent and comprehensive output csvs 
	- generates volcano plots as pdfs and svgs for reporting

### More info about supporting scripts [Here](../scripts/README_scriptsFolder.md)