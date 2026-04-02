computation resources for pipeline as configured in snakemake:
  - 32 core i9 intel processor
  - 32 gb of ram 
  - 16 cores allocated to snakemake call
  - max of 8 thread allocated to supportive jobs
  - min of 1 thread for alignment using star 

source files : 
- gtf: ~/genomeRefs/hg38/gencode.v47.annotation.gtf
- bed: ~/genomeRefs/hg38/gencode.v47.annotation.bed
  - generated from hg38 gencode version 47 gtf using:
    ```bash
      gtfToGenePred Homo_sapiens.GRCh38.113.gtf temp.genePred
      genePredToBed temp.genePred Homo_sapiens.GRCh38.113.bed
    ```
- hisat splice-aware indices: 
  ```bash
    wget --content-disposition ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz
  ```