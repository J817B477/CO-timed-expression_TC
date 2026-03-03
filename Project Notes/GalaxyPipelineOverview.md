### 02-07-2026
All raw fastq files were reprocessed based on single read protocols:
- the reads were organized into a single collection so that the analyses and transformations were outputted in their own collection
1. **fastQ file quality check**
	- analyses of read quality were conducted for the sake of defining appropriate parameters for read trimming
	1. `FastQC` was run in Galaxy, generating a html report with statistics for each fastq file and providing adapter sequences and over represented sequences detected in individual files
		- they consistently indicated that each fastQ file had inconsistent representation for the first 15 positions from the 5' end of the reads
	2. `MultiQC `was used to create an aggregated report for all of the fastQ files clarifying the adapters and over-represented sequences detected collectively
		- several sequences were indicated as over represented across the full set of fastQ files but there were no adapters found in the files
	3. `NCBI BLAST` was run on the over-represented sequences to confirm that they were native to the human genome ensuring that they represented true biological signal and were not artifacts of the procedure
		- all sequences were native to the human Genome
2. **Trimming**
	- Trimming was conducted using Trimmomatic with parameters chosen for trimming procedures base on the information 
	-   the following parameters were assessed for deployment and the decision to use or to not use were rationalized as follows:
		1. `ILLUMINACLIP` **<span style="color:rgb(163, 25, 25)">was not performed</span>** due to the lack of adapter presence
		2. `SLIDINGWINDOW` **<span style="color:rgb(163, 25, 25)">was not performed</span>** due to quality issues all being positional and falling within the first 15 positions of the 5' end of the reads
		3. `LEADING / TRAILING` trimming **<span style="color:rgb(163, 25, 25)">were also not performed</span>** for the same reason as `SLIDINGWINDOW`
		4. `CROP` **<span style="color:rgb(163, 25, 25)">was not performed</span>** as there was no indication of 3' end quality issues
		5. `HEADCROP` **<span style="color:rgb(0, 176, 80)">was performed</span>** to remove the first 15 positions of each read from their 5' ends to address quality issues perceived in the 
		6. `MINLEN` <span style="color:rgb(0, 176, 80)">was performed</span> subsequently to remove short reads (< 30 bp) after positional trimming to retain only the most reliable signal
		> **generally:** Trimming algorithms based on quality were avoided due to the fact that the reports did not indicate significant issues with quality that could not be handled by targeting the early positions of each read. While it may seem harmless, quality based trimming would increase length variability that is not related to the artifacts introduced in the laboratory procedure. The retention of reads of lengths greater than or equal 30 bp is to optimize the data for processing by `BWA MEM`
3. **Alignments**
	- were initially intended to be run using BWA MEM which is most effective on linger reads
	- unlike the original trimmed fastQ files that assumed paired reads, the new unpaired read fastQ files were significantly larger resulting in Out of Memory errors when alignment was run on them using BWA MEM in Galaxy - possibly due to memory restrictions for the web-based implementation
	- instead BWA was used for alignment which ran effectively but is less ideal for longer reads
	- to address this, the trimmed fastQ files will be downloaded and aligned locally with BWA MEM locally with more memory allocation