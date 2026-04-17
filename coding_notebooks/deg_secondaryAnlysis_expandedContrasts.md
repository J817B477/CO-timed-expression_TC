<script>  
window.MathJax = {  
tex: {  
inlineMath: [['$', '$'], ['\\(', '\\)']],  
displayMath: [['$$', '$$'], ['\\[', '\\]']]  
}  
};  
</script>  
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

The purpose of this notebook is to generate a pipeline to providing evidence for which gene groups and functional pathways played a role in the gene expression associated with the RNA-sequencing samples generated from treatments of hepatocyte cells with an undisclosed phytochemical treatment and with acetaminophen. The experimental treatment groups from which RNA libraries were produced consisted of hepatocyte cell line cultures with varying combinations of exposure to each treatment. This combination followed a preconditioning protocol where several of the treatment groups were primed with acetaminophen before exposure to the undisclosed phytochemical with a washout step in between. Time exposures were of 0, 4, 8, or 24 hour intervals. One control group consisted of replicates exposed to no treatments. The second consisted of replicates with 24 hours of exposure to only the culture medium only. 

Clarifying interactivity and independence between the two treatments was a central question in this experiment. Because of restrictions of the available data, the standard practice of defining an interactive testing design in DESeq2 was limited to experimental groups with a maximum of 4 hours of exposure to either treatment. This was because the inclusion of treatment groups with exposures beyond exposures of 4 hours resulted in full rank deficiency in the experimental design matrix upon which statistical testing relies in Differential Expression Analysis. Additionally, there was no 4 hour control for the medium available for contrasts. The results for hypotheses testing of interaction models derived from 0 to 4 hour exposure groups resulted in little indication of significance due to interaction. However, the analysis performed here indicates that the number of DEGs varied directly with the length of time the cultures were exposed to treatment prior to the generation of libraries indicating that a maximum 4 hour exposure is likely insufficient for assessing treatment interactivity. As a result of this, secondary analysis was performed in this notebook on data produced from differential expression testing results for contrasts generated between all of the individual experimental groups. 

#### Analytical Procedure:

The data used in this notebook's analysis consists of records for genes that were differentially expressed (referred to as DEGs) in one of the 28 experimental group contrasts generated in analysis of sample gene expression counts using R programming language's DESeq2 library (part of the Bioconductor suite of libraries for bioinformatics analysis). The hypothesis test design performed with DESeq2 used the distinct experimental groups as predictors in order to produce the measures for each gene record's attribute utilized in this analysis. 

The counts data used for differential analysis testing was initially produced from processing of RNA-sequencing fastq files that contain sequence reads along with metadata for each read. The fastq files were assessed for quality using FastQC and trimmed using fastp. The trimmed reads were aligned to the HISAT2 splice-aware index files for the GRCH38 human reference transcriptome (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz). Gene expression counts and relevant metadata were generated using FeatureCounts, and additional counts in the reportable transcripts per million format was generated using a custom script in R. **Only raw counts were used for differential expression analysis**. The RNA sequencing processing and differential analysis pipeline was navigated using Snakemake with supporting custom R and Python scripts and Conda coding environment for high level reproducibility. The results of the Differential Expression Analysis across all contrasts was separated into multiple comma-separated variable (csv) files: one containing the Benjamini-Hochberg adjusted p-values for each contrast, one containing the $log_{2}(fold change)$ values for each contrast, and a data frame of derived expression direction data. (This data was generated for volcano plots during differential expression analysis and maintained for use in this analysis). Each of these data frames served as source data for analysis performed.

**The steps of the analysis:**

1. The data pertaining to differentially expressed genes are separated into three perspective levels, each of which informs on the generation factors attributable to each gene. The significance level captures whether or not each DEG is found to be significant in each contrast. Here contrast-specific significance for each gene is defined at the 95% confidence level (adjusted p-value of less than .05). The direction level captures whether or not each DEG is up-regulated (found to be significantly increased respective to the base experimental group in each contrast) or down-regulated (found to be significantly decreased respective to each contrast). Finally, the effect size level perspective on the data uses each contrasts $log_{2}(Fold Change)$ generated for each gene in DESeq2 to quantify the magnitude that each DEG is altered by. (Fold change is just the factor by which the gene's expression level changes relative to the base of the contrast). 

1.  Atomic factors are generated to associate the impacts of the differential expression to the different treatments in the experiment. Because of the nature of the sets of contrasts available, there is confounding in both treatment and time. To mitigate this, contrasts and there differential expression profiles were allocated to informing on one the atomic factors so that these factors serve as proxy's for the different treatments in the experiment. Contrasts where there is deviation in time exposure to only one of the independent treatments (relative to the base of the contrast) is allocated to the factor serving as a proxy for that treatment. Contrasts where cultures are exposed to both treatments prior to RNA library development are allocated to the factor for dual exposure. 

1. Measures are generated respective to significance and direction levels of data including significance/up-regulation/down-regulation consensus scores. Additionally, entropy and sign balance are calculated along with a sign switching indicator. Each gene's consensus scores are simple binomial proportions for the instances of the gene's significance/up-regulation/down-regulation across the set of contrasts that inform on each factor. Each gene's significance and direction entropy scores measure the heterogeneity of the gene's expression characteristics across each factor's contrasts. Sign balance captures the proportionality of up and down regulation for a gene within each factor's contrast set indicating an overall direction profile for the factor, and sign switching indicates whether each gene has an instance of directional change in each factor's contrast set.

1. These measures are used, along with z-score-scaled effect size data, to establish a feature space for clustering the set of genes that were found to have statistically significant expression for at least one contrast articulated in the hypothesis testing performed in DESeq2. The clustering is conducted over two proximity measures for comparison of performance relative to the clustering algorithm's objective measures.

1. The better of the proximity measures, based on the cohesion/separation of clusters produced by each and the membership size distribution for each, is then selected for gene cluster enrichment analysis using Gene Ontology annotations for biological processes as stated in the "Genome wide annotation for Human" library (org.hs.eg.db version 3.20.0) that was queried using the Cluster Profiler library. 

1. The results of enrichment analysis is then concatenated into a data frame which is then exported as a csv for future clarification of per-cluster annotations. Several static visual tools are generated to connect the gene clusters to experimental factors, contrasts, and most significant enrichment profiles for each cluster. 

1. Association testing analysis is finally performed to map functions clarified by enrichment analysis to  the independent and combination treatment factor proxies. The testing determines if a function's association with an experimental factor is statistically significant. Also, consensus of promotion or reduction are evaluated for each function found to be associated with each atomic factor.

#### Limitations of Analysis:
The primary limitation of this analysis is that it can only approximate latent structure in the available data. The goal then is to establish factor associated patterns within the gene groupings through clustering analysis and suggestions of what functional behavior might be occurring in response to the experimental variables explored. In spite of the use of statistical testing, this analysis is limited to establishing similar response patterns across overlapping experimental contexts and cannot definitively establish separation between the effect of time and the effect of treatment type to a specific statistical significance threshold.


#### Technical information for running this notebook:  

The notebook requires the aforementioned csv's stored in this repository. The production of the csv's are as described above and can be ascertained in greater detail by reading the rules of the snakemake file and exploring the supportive scripts called by the rules. It is not possible to reproduce the source files for this notebook without the original fastq files. Access to original fastq files can only be accomplished through permission given by the principle investigator who assigned this analysis project to me. Currently, any inquiry for this access must be facilitated by the author of this notebook and its containing repository. (Reach out to me and I will reach out on your behalf). 

The notebook requires that you have an environment established which can run R in jupyter notebook. Generating the Conda environment from the provided conda_envs/jupyter.yaml will provide the exact capacity, but an environment with both jupyter notebooks and irkernel installed should be able to run this notebook after all dependencies below are installed. (This notebook required a separate conda environment due to version conflicts between its dependencies and those of the conda environment used for the snakemake pipeline, so be sure to use the one specified). All of the dependencies were established using mamba to facilitate agreement within the environment, and it is recommended that this is used to generate an environment that supports all of the dependencies below. I also recommend an IDE with strong support options for both R and python. The full workflow supporting this notebook and the notebook itself were developed in VSCode. 

This notebook uses relative pathways to reference read and written and files. It may be necessary to adjust these pathways depending on how you decide to store the source data. If you clone the full repository, you should have the project directory structure as needed to use existing path references, but remember that snake make and supporting files should be treated as read only, since they likely will be missing source files. 

# Contents

1. [Generation of Data Matrices](#i-generation-of-data-matrices)
1. [Extraction of Atomic and Experimental Factors Contrasts](#ii-extraction-of-atomic-and-experimental-factors-contrasts)
1. [Generation of Measures for Feature Space Dimensions](#iii-generation-of-measures-for-feature-space-dimensions)
1. [Generation of Feature Matrix](#iv-generation-of-feature-matrix)
1. [Clustering Analysis of Feature Space and Gene Enrichment Analysis](#v-clustering-analysis-of-feature-space-and-gene-enrichment-analysis)
1. [Static Visualizations for Navigating Characteristics of Gene Clusters](#vi-static-visualizations-for-navigating-characteristics-of-gene-clusters)
1. [Factor-Function Association Hypothesis Testing ](#vii-factor-function-association-hypothesis-testing)

_____
# I. Generation of Data Matrices
_____

##### [contents](#contents)

### Significance Matrix

    Number of Genes prior to Subsetting to DEGs:  78932 
    Number of Genes after to Subsetting to DEGs:  3782 
    
    Peek at significance Matrix: 



<table class="dataframe">
<caption>A matrix: 5 × 28 of type lgl</caption>
<thead>
	<tr><th></th><th scope=col>CO_4hr_APAP_24hr_vs_Control</th><th scope=col>CO_24hr_vs_Control</th><th scope=col>Control_24hr_vs_Control</th><th scope=col>CO_4hr_APAP_4hr_vs_Control</th><th scope=col>CO_8hr_APAP_8hr_vs_Control</th><th scope=col>CO_4hr_vs_Control</th><th scope=col>APAP_4hr_vs_Control</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_Control_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr</th><th scope=col>⋯</th><th scope=col>CO_4hr_APAP_4hr_vs_Control_24hr</th><th scope=col>CO_8hr_APAP_8hr_vs_Control_24hr</th><th scope=col>CO_4hr_vs_Control_24hr</th><th scope=col>APAP_4hr_vs_Control_24hr</th><th scope=col>CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr</th><th scope=col>CO_4hr_APAP_4hr_vs_CO_4hr</th><th scope=col>CO_4hr_APAP_4hr_vs_APAP_4hr</th><th scope=col>CO_8hr_APAP_8hr_vs_CO_4hr</th><th scope=col>CO_8hr_APAP_8hr_vs_APAP_4hr</th><th scope=col>APAP_4hr_vs_CO_4hr</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003</th><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001167</th><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001460</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001461</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>⋯</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001497</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
</tbody>
</table>



#### The significance matrix consists of boolean values indicating where each gene in the DEG is differentially expressed relative to the base of each contrast. True values represent all instances where the $\text{adjusted p-value} \leq 0.05$. The data informs on measures created for downstream gene clustering analysis. The data was subset to include only those genes that were found to have significant differential expression in at least one contrast. This reduced the overall number of transcription from 78,932 to 3782. All of the subsequent analysis is based on this subset, referred to as the DEG set.

### Directional Matrix


<table class="dataframe">
<caption>A data.frame: 6 × 28</caption>
<thead>
	<tr><th></th><th scope=col>CO_4hr_APAP_24hr_vs_Control</th><th scope=col>CO_24hr_vs_Control</th><th scope=col>Control_24hr_vs_Control</th><th scope=col>CO_4hr_APAP_4hr_vs_Control</th><th scope=col>CO_8hr_APAP_8hr_vs_Control</th><th scope=col>CO_4hr_vs_Control</th><th scope=col>APAP_4hr_vs_Control</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_Control_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr</th><th scope=col>⋯</th><th scope=col>CO_4hr_APAP_4hr_vs_Control_24hr</th><th scope=col>CO_8hr_APAP_8hr_vs_Control_24hr</th><th scope=col>CO_4hr_vs_Control_24hr</th><th scope=col>APAP_4hr_vs_Control_24hr</th><th scope=col>CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr</th><th scope=col>CO_4hr_APAP_4hr_vs_CO_4hr</th><th scope=col>CO_4hr_APAP_4hr_vs_APAP_4hr</th><th scope=col>CO_8hr_APAP_8hr_vs_CO_4hr</th><th scope=col>CO_8hr_APAP_8hr_vs_APAP_4hr</th><th scope=col>APAP_4hr_vs_CO_4hr</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003</th><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0 </td><td>0 </td><td>0 </td><td>⋯</td><td>0 </td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000001167</th><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0 </td><td>0 </td><td>-1</td><td>⋯</td><td>0 </td><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000001460</th><td>1 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1 </td><td>1 </td><td>1 </td><td>⋯</td><td>0 </td><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000001461</th><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>-1</td><td>-1</td><td>0 </td><td>⋯</td><td>-1</td><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000001497</th><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0 </td><td>0 </td><td>0 </td><td>⋯</td><td>0 </td><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000001617</th><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0 </td><td>0 </td><td>0 </td><td>⋯</td><td>0 </td><td>0 </td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




3782


#### The directional matrix provides information about the direction of change for genes shown to be significant in the significance matrix. Its values consist of $1$ indicating that expression was increased relative to each contrast base, $-1$  indicating that expression was decreased relative to each base, and $0$ indicating that the gene didn't have significant expression change relative to each base.  

### Effect Size Matrix


<table class="dataframe">
<caption>A data.frame: 6 × 28</caption>
<thead>
	<tr><th></th><th scope=col>CO_4hr_APAP_24hr_vs_Control</th><th scope=col>CO_24hr_vs_Control</th><th scope=col>Control_24hr_vs_Control</th><th scope=col>CO_4hr_APAP_4hr_vs_Control</th><th scope=col>CO_8hr_APAP_8hr_vs_Control</th><th scope=col>CO_4hr_vs_Control</th><th scope=col>APAP_4hr_vs_Control</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_Control_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr</th><th scope=col>⋯</th><th scope=col>CO_4hr_APAP_4hr_vs_Control_24hr</th><th scope=col>CO_8hr_APAP_8hr_vs_Control_24hr</th><th scope=col>CO_4hr_vs_Control_24hr</th><th scope=col>APAP_4hr_vs_Control_24hr</th><th scope=col>CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr</th><th scope=col>CO_4hr_APAP_4hr_vs_CO_4hr</th><th scope=col>CO_4hr_APAP_4hr_vs_APAP_4hr</th><th scope=col>CO_8hr_APAP_8hr_vs_CO_4hr</th><th scope=col>CO_8hr_APAP_8hr_vs_APAP_4hr</th><th scope=col>APAP_4hr_vs_CO_4hr</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003</th><td> 0.9718274</td><td> 0.7770204</td><td> 1.4152657</td><td> 0.70568538</td><td>-0.009971626</td><td> 0.2596122</td><td> 0.285636505</td><td> 0.194807</td><td>-0.4434383</td><td> 0.2661420</td><td>⋯</td><td>-0.70958028</td><td>-1.42523728</td><td>-1.15565344</td><td>-1.1296291</td><td>-0.71565700</td><td> 0.44607316</td><td> 0.42004887</td><td>-0.26958384</td><td>-0.2956081</td><td> 0.02602429</td></tr>
	<tr><th scope=row>ENSG00000001167</th><td>-2.6757986</td><td> 1.0579827</td><td> 0.1807754</td><td> 1.91842513</td><td> 0.796376095</td><td> 1.1181610</td><td> 0.581552462</td><td>-3.733781</td><td>-2.8565740</td><td>-4.5942237</td><td>⋯</td><td> 1.73764971</td><td> 0.61560068</td><td> 0.93738554</td><td> 0.4007770</td><td>-1.12204903</td><td> 0.80026417</td><td> 1.33687267</td><td>-0.32178486</td><td> 0.2148236</td><td>-0.53660849</td></tr>
	<tr><th scope=row>ENSG00000001460</th><td> 3.3602216</td><td> 0.5098319</td><td>-1.1367617</td><td> 0.00447294</td><td> 0.700851942</td><td> 1.2548651</td><td> 0.055667483</td><td> 2.850390</td><td> 4.4969833</td><td> 3.3557487</td><td>⋯</td><td> 1.14123465</td><td> 1.83761365</td><td> 2.39162677</td><td> 1.1924292</td><td> 0.69637900</td><td>-1.25039212</td><td>-0.05119454</td><td>-0.55401312</td><td> 0.6451845</td><td>-1.19919758</td></tr>
	<tr><th scope=row>ENSG00000001461</th><td>-4.2693866</td><td> 0.1707683</td><td> 0.8599573</td><td>-1.87896157</td><td>-0.699149394</td><td>-0.6235310</td><td>-0.494020428</td><td>-4.440155</td><td>-5.1293440</td><td>-2.3904251</td><td>⋯</td><td>-2.73891891</td><td>-1.55910674</td><td>-1.48348831</td><td>-1.3539778</td><td> 1.17981217</td><td>-1.25543060</td><td>-1.38494114</td><td>-0.07561843</td><td>-0.2051290</td><td> 0.12951054</td></tr>
	<tr><th scope=row>ENSG00000001497</th><td>-4.3205230</td><td>-0.5098104</td><td>-0.6231656</td><td>-0.64657533</td><td>-0.586729007</td><td>-0.7131603</td><td>-0.007529912</td><td>-3.810713</td><td>-3.6973575</td><td>-3.6739477</td><td>⋯</td><td>-0.02340977</td><td> 0.03643655</td><td>-0.08999471</td><td> 0.6156356</td><td> 0.05984633</td><td> 0.06658493</td><td>-0.63904542</td><td> 0.12643126</td><td>-0.5791991</td><td> 0.70563035</td></tr>
	<tr><th scope=row>ENSG00000001617</th><td> 1.0305880</td><td>-0.3980008</td><td>-0.5706939</td><td> 0.51066571</td><td>-0.853496482</td><td>-0.0308911</td><td> 0.323684410</td><td> 1.428589</td><td> 1.6012819</td><td> 0.5199223</td><td>⋯</td><td> 1.08135959</td><td>-0.28280260</td><td> 0.53980278</td><td> 0.8943783</td><td>-1.36416219</td><td> 0.54155681</td><td> 0.18698130</td><td>-0.82260538</td><td>-1.1771809</td><td> 0.35457551</td></tr>
</tbody>
</table>




3782


#### The effect size matrix includes each DEG set gene's $log_{2} \text{ fold change}$ measure for each contrast. This measure indicates the amount each gene changed espective to the base of each contrast. Fold change refers to the multiple that each expression count changed by and the sign of the fold change indicates the direction of the change (i.e. whether the gene was up or down regulated). Due to the range in magnitude of the fold changes they are reported as base 2 logarithms. This matrix served as the foundation of the feature space used for the clustering of DEG set genes.

## 1. Table and Bar Plot: frequency distribution of GEDs per contrast

    
    
    contrast                               count
    ------------------------------------  ------
    CO_4hr_APAP_24hr_vs_Control             1316
    CO_24hr_vs_Control                       410
    Control_24hr_vs_Control                  332
    CO_4hr_APAP_4hr_vs_Control               202
    CO_8hr_APAP_8hr_vs_Control               272
    CO_4hr_vs_Control                         31
    APAP_4hr_vs_Control                      112
    CO_4hr_APAP_24hr_vs_CO_24hr             1201
    CO_4hr_APAP_24hr_vs_Control_24hr        1095
    CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr     1308
    CO_4hr_APAP_24hr_vs_CO_8hr_APAP_8hr     1417
    CO_4hr_APAP_24hr_vs_CO_4hr              1391
    CO_4hr_APAP_24hr_vs_APAP_4hr            1084
    CO_24hr_vs_Control_24hr                  318
    CO_24hr_vs_CO_4hr_APAP_4hr               863
    CO_24hr_vs_CO_8hr_APAP_8hr               740
    CO_24hr_vs_CO_4hr                        549
    CO_24hr_vs_APAP_4hr                      668
    CO_4hr_APAP_4hr_vs_Control_24hr          697
    CO_8hr_APAP_8hr_vs_Control_24hr          738
    CO_4hr_vs_Control_24hr                   479
    APAP_4hr_vs_Control_24hr                 484
    CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr        97
    CO_4hr_APAP_4hr_vs_CO_4hr                 65
    CO_4hr_APAP_4hr_vs_APAP_4hr               30
    CO_8hr_APAP_8hr_vs_CO_4hr                 99
    CO_8hr_APAP_8hr_vs_APAP_4hr              141
    APAP_4hr_vs_CO_4hr                        68



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_16_1.png)
    


The above table and corresponding bar chart clarify how many of the DEG set genes were found to be significant in each of the differential expression contrasts. The disparity in the number of differentially expressed genes in contrasts with long exposures vs short exposures suggests that exposure time influences the gene expression profiles of each contrasts.

## 2. Heat maps of DEG candidates per contrast (significance, direction, effect size).

    `use_raster` is automatically set to TRUE for a matrix with more than
    2000 columns You can control `use_raster` argument by explicitly
    setting TRUE/FALSE to it.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    
    'magick' package is suggested to install to give better rasterization.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_19_1.png)
    


Figure 1. Is a hierarchically clustered heatmap of the full set DEGs that were found in at least one of the contrasts established through differential expression analysis hypothesis testing. The clustering algorithm groups genes by similarity in the binary distribution of significant expression across all of the samples. Because the data is binary, the default euclidean distance proximity measure was not used. Cosine similarity was chosen over Jaccard and Hamming as a proximity measure for genes to be clustered on (dendrogram not shown). It ultimately defines the  difference between genes as the cosine of the angle between each gene when they are represented in the space formed by axes defined by the contrasts. It is directly proportional to the angle itself but can be derived without finding the angle using: 

$$cos(\theta) = \frac{\vec{x}\cdot\vec{y}}{\|\vec{x}\|\|\vec{y}\|}$$

where the vector is effectively the direction values across the contrasts for the rows representing each gene.

Cosine ensures that genes are compared in their overall pattern of direction of expression while minimizing the influence of the number of times that they are perturbed and in doing so makes the direction perturbation more impactful. Grouping the genes allows for trends in DEGs instances of significance to be be observed generally across all contrasts. This in turn illustrates some trends in expression by genes that are grouped by similarity in expression pattern. This provides perspective for the experimental groups in the contrasts. 

Jaccard similarity is used for the clustering of contrasts. This ensures contrasts are grouped by agreement on what genes are DEGs rather than what genes are not so that lack of significant expression does not drive similarity assessed between contrasts. 

The information from the heatmap derived from data at this level of specificity is restricted to a broad perspective on significant expression trends across contrasts. While this might be anticipated logically, the heatmap clarifies that exposure time is a significant determining factor in the number of differentially expressed genes. The contrasts containing experimental groups with longer exposures have more significant differential expression than the contrasts containing only experimental groups with shorter exposures. 

    `use_raster` is automatically set to TRUE for a matrix with more than
    2000 columns You can control `use_raster` argument by explicitly
    setting TRUE/FALSE to it.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    
    'magick' package is suggested to install to give better rasterization.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_21_1.png)
    


Figure 2. highlights the direction of the change in a genes expression captured by the contrasts. This level of specificity provides indications of trends for the manner in which the expression of a gene is significantly differentiated for a given contrast. For instance, long exposures to experimental treatments tend to have an overall suppressive effect in the contrasts in which they are a factor relative to shorter exposures. However, shorter exposures tend to up-regulate more relative to the 24 hour control. This suggests that the medium itself might have an overall suppressive effect on gene expression. In fact when looking at the contrasts with long treatment exposures in contrast to the 24 hour control, the significance in the contrast seems to be less overall than when contrasts with experimental groups with shorter overall exposures. The contrasts of experimental groups with lower exposure times show a higher level of up-regulation than down regulation within the set of genes that tend to be differentially expressed at these intervals of exposure. 

Regions within the set of DEGs also offer insight into overall directional expression differentiation. There is a region of genes (at the end of the heatmap) that has a fundamentally different significance profile in contrasts with long exposure to acetaminophen than is seen in contrasts containing long exposure to the phytochemical. Specifically, when there is high instances of significance seen in this region for long exposure to one of the treatments, there is low instances of significance seen with long exposure to the other treatment. 

    `use_raster` is automatically set to TRUE for a matrix with more than
    2000 columns You can control `use_raster` argument by explicitly
    setting TRUE/FALSE to it.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    
    'magick' package is suggested to install to give better rasterization.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_23_1.png)
    


Figure 3. provides another level of specificity in the information provided through the heat mapping of the fold change of each gene within each contrast. This representation obscures the statistical significance of gene expression that was explicit Figures 1 and 2, as seen in the light coloration in regions that are grey in the prior heatmaps. However, insight into the magnitude of change introduced by of visual encoding of the fold change provides some additional context. For instance, even though the differential expression captured in the contrast between the 4 hour bi-treatment exposure and the 4 hour exposure to only acetaminophen produced few genes that were significantly differentially expressed, the heatmap shows that there was a consistent trend of down-regulation for genes in this contrast. 

_____
# II. Extraction of Atomic and Experimental Factors Contrasts
_____

##### [contents](#contents)

## 3. Metadata tables for atomic factorization of contrast groups

| CO | APAP | CONDITIONED | CO_4hr_APAP_24hr | CO_24hr | Control_24hr | CO_4hr_APAP_4hr | CO_8hr_APAP_8hr | CO_4hr | APAP_4hr | Control |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |

<!-- | CO | APAP || CO_4hr | APAP_4hr | CO_4hr_APAP_4hr | CO_8hr_APAP_8hr | CO_4hr_APAP_24hr | CO_24hr | medium_24hr |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- |

| Control 
| -------- | -------- | -------- |  -->


#### This table shows the results of systematic allocation of contrasts to each factor. Contrasts where the exposure between the base and compared experimental group only had differences in exposure time to acetaminophen was allocated to the apap factor and contrasts with the same characteristic for the phytochemical treatment is allocated to the co factor. The remaining contrasts were allocated to the conditioned factor to inform on the effect of the dual treatment exposure of some the experimental groups. 

#### In order to ensure that there was no replicate signatures for the derived factors across the contrasts, experimental groups were also included along with the set of derived atomic factors. The atomic factors represent the contrasts most likely to reveal information about gene expression related to one of the specific treatments or their combination where acetaminophen primes cell line cultures prior to the phytochemical treatment. This factorization forms the foundation of the analysis. Any inclusion of the exerimental group factors are for potential comparative reference.


<table class="dataframe">
<caption>A data.frame: 28 × 11</caption>
<thead>
	<tr><th></th><th scope=col>CO</th><th scope=col>APAP</th><th scope=col>CONDITIONED</th><th scope=col>CO_4hr_APAP_24hr</th><th scope=col>CO_24hr</th><th scope=col>Control_24hr</th><th scope=col>CO_4hr_APAP_4hr</th><th scope=col>CO_8hr_APAP_8hr</th><th scope=col>CO_4hr</th><th scope=col>APAP_4hr</th><th scope=col>Control</th></tr>
	<tr><th></th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th><th scope=col>&lt;list&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_Control</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>CO_24hr_vs_Control</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>Control_24hr_vs_Control</th><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_4hr_vs_Control</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>CO_8hr_APAP_8hr_vs_Control</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>CO_4hr_vs_Control</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>APAP_4hr_vs_Control</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_CO_24hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_Control_24hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_CO_8hr_APAP_8hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_CO_4hr</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr_vs_APAP_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_24hr_vs_Control_24hr</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_24hr_vs_CO_4hr_APAP_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_24hr_vs_CO_8hr_APAP_8hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_24hr_vs_CO_4hr</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_24hr_vs_APAP_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_4hr_vs_Control_24hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_8hr_APAP_8hr_vs_Control_24hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_vs_Control_24hr</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>APAP_4hr_vs_Control_24hr</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_4hr_vs_CO_4hr</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_4hr_APAP_4hr_vs_APAP_4hr</th><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_8hr_APAP_8hr_vs_CO_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CO_8hr_APAP_8hr_vs_APAP_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>APAP_4hr_vs_CO_4hr</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
</tbody>
</table>



### This table provides readable collections of contrasts associated with each factor as established by the true/false sequences in the rows of the above table.


<table class="dataframe">
<caption>A data.frame: 11 × 3</caption>
<thead>
	<tr><th scope=col>factor</th><th scope=col>contrasts</th><th scope=col>length</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>CO              </td><td>CO_24hr_vs_Control, CO_4hr_vs_Control, CO_24hr_vs_Control_24hr, CO_24hr_vs_CO_4hr, CO_4hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_APAP_4hr                                                                                                                                                                                                                                                                                                                                              </td><td> 6</td></tr>
	<tr><td>APAP            </td><td>APAP_4hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr, CO_4hr_APAP_24hr_vs_CO_4hr, APAP_4hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_CO_4hr                                                                                                                                                                                                                                                                                                                                           </td><td> 5</td></tr>
	<tr><td>CONDITIONED     </td><td>CO_4hr_APAP_24hr_vs_Control, CO_4hr_APAP_4hr_vs_Control, CO_8hr_APAP_8hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_24hr, CO_4hr_APAP_24hr_vs_Control_24hr, CO_4hr_APAP_24hr_vs_CO_8hr_APAP_8hr, CO_4hr_APAP_24hr_vs_APAP_4hr, CO_24hr_vs_CO_4hr_APAP_4hr, CO_24hr_vs_CO_8hr_APAP_8hr, CO_24hr_vs_APAP_4hr, CO_4hr_APAP_4hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr, CO_8hr_APAP_8hr_vs_CO_4hr, CO_8hr_APAP_8hr_vs_APAP_4hr, APAP_4hr_vs_CO_4hr</td><td>16</td></tr>
	<tr><td>CO_4hr_APAP_24hr</td><td>CO_4hr_APAP_24hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_24hr, CO_4hr_APAP_24hr_vs_Control_24hr, CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr, CO_4hr_APAP_24hr_vs_CO_8hr_APAP_8hr, CO_4hr_APAP_24hr_vs_CO_4hr, CO_4hr_APAP_24hr_vs_APAP_4hr                                                                                                                                                                                                                                                      </td><td> 7</td></tr>
	<tr><td>CO_24hr         </td><td>CO_24hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_24hr, CO_24hr_vs_Control_24hr, CO_24hr_vs_CO_4hr_APAP_4hr, CO_24hr_vs_CO_8hr_APAP_8hr, CO_24hr_vs_CO_4hr, CO_24hr_vs_APAP_4hr                                                                                                                                                                                                                                                                                                            </td><td> 7</td></tr>
	<tr><td>Control_24hr    </td><td>Control_24hr_vs_Control, CO_4hr_APAP_24hr_vs_Control_24hr, CO_24hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_Control_24hr, CO_4hr_vs_Control_24hr, APAP_4hr_vs_Control_24hr                                                                                                                                                                                                                                                                              </td><td> 7</td></tr>
	<tr><td>CO_4hr_APAP_4hr </td><td>CO_4hr_APAP_4hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr, CO_24hr_vs_CO_4hr_APAP_4hr, CO_4hr_APAP_4hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr, CO_4hr_APAP_4hr_vs_CO_4hr, CO_4hr_APAP_4hr_vs_APAP_4hr                                                                                                                                                                                                                                                            </td><td> 7</td></tr>
	<tr><td>CO_8hr_APAP_8hr </td><td>CO_8hr_APAP_8hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_8hr_APAP_8hr, CO_24hr_vs_CO_8hr_APAP_8hr, CO_8hr_APAP_8hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr, CO_8hr_APAP_8hr_vs_CO_4hr, CO_8hr_APAP_8hr_vs_APAP_4hr                                                                                                                                                                                                                                                            </td><td> 7</td></tr>
	<tr><td>CO_4hr          </td><td>CO_4hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_4hr, CO_24hr_vs_CO_4hr, CO_4hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_CO_4hr, CO_8hr_APAP_8hr_vs_CO_4hr, APAP_4hr_vs_CO_4hr                                                                                                                                                                                                                                                                                                                  </td><td> 7</td></tr>
	<tr><td>APAP_4hr        </td><td>APAP_4hr_vs_Control, CO_4hr_APAP_24hr_vs_APAP_4hr, CO_24hr_vs_APAP_4hr, APAP_4hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_APAP_4hr, CO_8hr_APAP_8hr_vs_APAP_4hr, APAP_4hr_vs_CO_4hr                                                                                                                                                                                                                                                                                                      </td><td> 7</td></tr>
	<tr><td>Control         </td><td>CO_4hr_APAP_24hr_vs_Control, CO_24hr_vs_Control, Control_24hr_vs_Control, CO_4hr_APAP_4hr_vs_Control, CO_8hr_APAP_8hr_vs_Control, CO_4hr_vs_Control, APAP_4hr_vs_Control                                                                                                                                                                                                                                                                                                            </td><td> 7</td></tr>
</tbody>
</table>



_____
# III. Generation of Measures for Feature Space Dimensions
_____

##### [contents](#contents)

#### The measures generated here are designed to reveal the nature of each DEG set gene's expression, respective to each factor. These measures for the atomic factors are designed to help biological function enrichment in downstream analysis. The extension of these measures to all factors is to provide additional dimensions for performing clustering analysis.

## 4. Measure Definitions
1. **Significance Consensus (Factor Level)**
$$
  f(g) = \frac{\Sigma_{c=1}^{N}I(g\in S_{c})}{N}
$$

  - $I$ is a logical function that indicates if the condition is met
  - $\Sigma$ effectively counts the number of successes
  - $c$ refers to one of the contrasts in the set that defined a specific factor
  - $g$ is a gene
  - $N$ is $N(C_{f})$ where $C_{f}$ is the set of contrasts for a given factor 

2. **Entropy (Gene Level: not Factor Level)**
  $$
    H(g) = -\Sigma_{f=1}^{N}p_{g,f}log(p_{g,f})
  $$
- $N$ is $N(f\in F)$ (size of the full set of factors)
- $p_{g,f}$ is $\frac{f_{i} * n_{i}}{\Sigma f_{i} * n_{i}}$ (proportion weighted by number of factor-associated contrasts)

3. **Sign Balance (Factor level: Not Gene Level)**
  $$
    B(g) = \frac{N(g = 1) - N(g = -1)}{N(g = 1) + N(g = -1)}
  $$

4. **Sign Switch (Factor Level: not Gene Level)**
  $$
  I(g_{c}) = \left\{ 
  \begin{array}{ll} 
  1 & \text{if } \{-1, 1\} \subseteq c \\ 
  0 & \text{otherwise} 
  \end{array} 
  \right.
  $$


## 5. Generation of Consensus Values (significance and direction)

    [1] "Five Point Summary Distribution of Significance Consensus per Factor"



           CO               APAP         CONDITIONED     CO_4hr_APAP_24hr
     Min.   :0.00000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
     1st Qu.:0.00000   1st Qu.:0.0000   1st Qu.:0.0625   1st Qu.:0.0000  
     Median :0.00000   Median :0.2000   Median :0.1250   Median :0.1429  
     Mean   :0.08007   Mean   :0.1777   Mean   :0.1768   Mean   :0.3329  
     3rd Qu.:0.16667   3rd Qu.:0.4000   3rd Qu.:0.2500   3rd Qu.:0.5714  
     Max.   :1.00000   Max.   :1.0000   Max.   :0.8750   Max.   :1.0000  
        CO_24hr        Control_24hr    CO_4hr_APAP_4hr  CO_8hr_APAP_8hr 
     Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
     1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
     Median :0.1429   Median :0.1429   Median :0.1429   Median :0.1429  
     Mean   :0.1794   Mean   :0.1565   Mean   :0.1232   Mean   :0.1324  
     3rd Qu.:0.2857   3rd Qu.:0.1429   3rd Qu.:0.1429   3rd Qu.:0.1429  
     Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
         CO_4hr          APAP_4hr          Control      
     Min.   :0.0000   Min.   :0.00000   Min.   :0.0000  
     1st Qu.:0.0000   1st Qu.:0.00000   1st Qu.:0.0000  
     Median :0.1429   Median :0.00000   Median :0.1429  
     Mean   :0.1013   Mean   :0.09772   Mean   :0.1010  
     3rd Qu.:0.1429   3rd Qu.:0.14286   3rd Qu.:0.1429  
     Max.   :1.0000   Max.   :0.85714   Max.   :1.0000  



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_40_2.png)
    


The gene consensus scores distribution per contrasts factor indicate the extent to which each gene is significantly expressed in the of contrasts that share a factor. High scores for that factor indicate a high level of significant DEGs associated with that factor across contrasts while low scores indicate few of the set of the total set of differentially expressed genes are significantly expressed in contrasts associated with that factor. The vertical boundaries of the boxes indicate the inter-quartile range for each factor's distribution with the lower bound indicating the 25th percentile value (25% of all values in the distribution are found to be below the value at the bottom of the box) percentile and the higher bound indicating the 75th percentile value. The line across the middle indicates the boundary of 50% of the data. If the line is high it means that half of the data is concentrated in the high values found in the inner quartile range while lower concentrates more data in the lower values of the IQR. Longer boxes indicate a wider IQR, meaning that the distribution is more heterogeneous while shorter boxes are more homogeneous. The over all height of the boxes clarify the amount of overall significance associated with that contrast. Because data reflects frequency of significance in contrasts associated with these factors, there is not a significant amount of insight to be gained regarding biological signal but it provides some perspective on the factors' ability to define a feature space that supports successful clustering of genes as measured by distinctive functional associations- higher consensus score indicate higher association between a gene's expression and the factor associated with it. It still remains to be seen if the specific genes with higher scores will reveal that treatment parameters are associated with specific functional pathways. 

    [1] "Five Point Summary Distribution of Down-regulation Consensus per Factor"



           CO               APAP          CONDITIONED     CO_4hr_APAP_24hr
     Min.   :0.00000   Min.   :0.00000   Min.   :0.0000   Min.   :0.0000  
     1st Qu.:0.00000   1st Qu.:0.00000   1st Qu.:0.0000   1st Qu.:0.0000  
     Median :0.00000   Median :0.00000   Median :0.0625   Median :0.0000  
     Mean   :0.03962   Mean   :0.09122   Mean   :0.0940   Mean   :0.1922  
     3rd Qu.:0.00000   3rd Qu.:0.20000   3rd Qu.:0.1250   3rd Qu.:0.2857  
     Max.   :0.66667   Max.   :0.80000   Max.   :0.6250   Max.   :1.0000  
        CO_24hr        Control_24hr    CO_4hr_APAP_4hr   CO_8hr_APAP_8hr  
     Min.   :0.0000   Min.   :0.0000   Min.   :0.00000   Min.   :0.00000  
     1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.00000   1st Qu.:0.00000  
     Median :0.0000   Median :0.0000   Median :0.00000   Median :0.00000  
     Mean   :0.1069   Mean   :0.0672   Mean   :0.06346   Mean   :0.06463  
     3rd Qu.:0.1429   3rd Qu.:0.1429   3rd Qu.:0.14286   3rd Qu.:0.14286  
     Max.   :0.8571   Max.   :0.8571   Max.   :0.57143   Max.   :0.85714  
         CO_4hr           APAP_4hr          Control       
     Min.   :0.00000   Min.   :0.00000   Min.   :0.00000  
     1st Qu.:0.00000   1st Qu.:0.00000   1st Qu.:0.00000  
     Median :0.00000   Median :0.00000   Median :0.00000  
     Mean   :0.05092   Mean   :0.04941   Mean   :0.04971  
     3rd Qu.:0.14286   3rd Qu.:0.14286   3rd Qu.:0.14286  
     Max.   :0.42857   Max.   :0.42857   Max.   :0.57143  



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_42_2.png)
    


    [1] "Five Point Summary Distribution of Up-regulation Consensus per Factor"



       downreg_CO       downreg_APAP     downreg_CONDITIONED
     Min.   :0.00000   Min.   :0.00000   Min.   :0.0000     
     1st Qu.:0.00000   1st Qu.:0.00000   1st Qu.:0.0000     
     Median :0.00000   Median :0.00000   Median :0.0625     
     Mean   :0.03962   Mean   :0.09122   Mean   :0.0940     
     3rd Qu.:0.00000   3rd Qu.:0.20000   3rd Qu.:0.1250     
     Max.   :0.66667   Max.   :0.80000   Max.   :0.6250     
     downreg_CO_4hr_APAP_24hr downreg_CO_24hr  downreg_Control_24hr
     Min.   :0.0000           Min.   :0.0000   Min.   :0.0000      
     1st Qu.:0.0000           1st Qu.:0.0000   1st Qu.:0.0000      
     Median :0.0000           Median :0.0000   Median :0.0000      
     Mean   :0.1922           Mean   :0.1069   Mean   :0.0672      
     3rd Qu.:0.2857           3rd Qu.:0.1429   3rd Qu.:0.1429      
     Max.   :1.0000           Max.   :0.8571   Max.   :0.8571      
     downreg_CO_4hr_APAP_4hr downreg_CO_8hr_APAP_8hr downreg_CO_4hr   
     Min.   :0.00000         Min.   :0.00000         Min.   :0.00000  
     1st Qu.:0.00000         1st Qu.:0.00000         1st Qu.:0.00000  
     Median :0.00000         Median :0.00000         Median :0.00000  
     Mean   :0.06346         Mean   :0.06463         Mean   :0.05092  
     3rd Qu.:0.14286         3rd Qu.:0.14286         3rd Qu.:0.14286  
     Max.   :0.57143         Max.   :0.85714         Max.   :0.42857  
     downreg_APAP_4hr  downreg_Control  
     Min.   :0.00000   Min.   :0.00000  
     1st Qu.:0.00000   1st Qu.:0.00000  
     Median :0.00000   Median :0.00000  
     Mean   :0.04941   Mean   :0.04971  
     3rd Qu.:0.14286   3rd Qu.:0.14286  
     Max.   :0.42857   Max.   :0.57143  



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_43_2.png)
    


The distributions of both up and down regulation consensus measures are very similar with the most meaningful distinction between distributions found in the outliers and the 0-hour control. The common distribution of the atomic factors suggests a balanced distribution of how and to what extent genes are shown to be significantly expressed across factors.

## 6. Generation of Entropy Scores (Cross Factor Weighted Aggregated for each Gene)


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      1.000   1.585   2.419   2.340   2.808   3.373 



       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0.000   0.000   1.585   1.469   2.522   3.318 



       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0.000   0.000   1.585   1.372   2.419   3.322 



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_47_0.png)
    


Entropy captures the dispersion of the factor scores for each gene. If the area of the boxes is small, the gene's factor score entropy is itself homogeneous, meaning the range of scores across genes are closer in value score than for larger area boxes. High scores indicate broad distribution across factors meaning it is not highly associated with any one factor. In other words, it is difficult to pinpoint any factor or subset of factors that explain its significant expression. Low entropy on the other hand indicates a much higher level of association between its expression and specific factors. Entropy of factor scores for significance indicates uniformity of expression across test contrasts indicating that the effect of experimental factors is higher for direction of significant change than for significance alone. The median entropy is high within the inner quartile range of both induced and repressed DEGs with the implication being that 50% of genes have high entropy scores indicating a consistent directional pattern across contrasts. An even higher percentage are found to express consistently with respect to significance and direction, regardless of the experimental conditions. 

## 7. Generation of Sign Balance Score (for each factor of each gene)


    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_50_0.png)
    


The factor response bias plot shows consistent balance across factors with a slight skew towards down regulation for all factors except the phytochemical atomic factor. The potential range for sign balance is $1$, indicating that DEGs were up-regulated in all contrasts, and $-1$, indicating that DEGs were down-regulated in all cases. The maximum observed magnitude of $|-0.07|$ in the CO 24 hour and APAP 24 hour experimental group factor indicates that while there is a consistent imbalance toward repression of gene expression. There is seemingly contradictory information here but it may reflect how the factors are aggregated where the experimental group factors serve primarily as background information. Downstream analysis for biological function associations supports this overall sign balance for phytochemical influence. 

    `use_raster` is automatically set to TRUE for a matrix with more than
    2000 columns You can control `use_raster` argument by explicitly
    setting TRUE/FALSE to it.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    
    'magick' package is suggested to install to give better rasterization.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_52_1.png)
    


##### The heat map shows each DEG's sign balance for the contrasts belonging to each factor. Many genes are shown to have very strong bias for expression induction or repression clarifying that the net sign balance across factors is due in large part to the balance of frequency in which a gene is consistently up-regulated or down-regulated. There are many genes that vary in magnitude of sign bias across the contrasts associated with each factor and others that experience directional shifts in expression perturbation. These genes will likely provide the most insight into what influence the different factors have on the hepatocyte cells being studied. 

## 8. Generates Switching Indicators (for each Factor of each Gene)


<table class="dataframe">
<caption>A data.frame: 11 × 5</caption>
<thead>
	<tr><th scope=col>factor</th><th scope=col>proportion_switch</th><th scope=col>Number of Contrasts</th><th scope=col>expected_switch</th><th scope=col>bio_significance_ratio</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>switch_CO              </td><td>0.07588577</td><td> 6</td><td>0.056319408</td><td> 1.3474178</td></tr>
	<tr><td>switch_APAP            </td><td>0.24034902</td><td> 5</td><td>0.048387097</td><td> 4.9672131</td></tr>
	<tr><td>switch_CONDITIONED     </td><td>0.18085669</td><td>16</td><td>0.238762559</td><td> 0.7574751</td></tr>
	<tr><td>switch_CO_4hr_APAP_24hr</td><td>0.09941830</td><td> 7</td><td>0.004230566</td><td>23.5000000</td></tr>
	<tr><td>switch_CO_24hr         </td><td>0.09968271</td><td> 7</td><td>0.096774194</td><td> 1.0300546</td></tr>
	<tr><td>switch_Control_24hr    </td><td>0.08699101</td><td> 7</td><td>0.089106293</td><td> 0.9762611</td></tr>
	<tr><td>switch_CO_4hr_APAP_4hr </td><td>0.10206240</td><td> 7</td><td>0.130618720</td><td> 0.7813765</td></tr>
	<tr><td>switch_CO_8hr_APAP_8hr </td><td>0.09836066</td><td> 7</td><td>0.135642517</td><td> 0.7251462</td></tr>
	<tr><td>switch_CO_4hr          </td><td>0.08143839</td><td> 7</td><td>0.084875727</td><td> 0.9595016</td></tr>
	<tr><td>switch_APAP_4hr        </td><td>0.07905870</td><td> 7</td><td>0.092543628</td><td> 0.8542857</td></tr>
	<tr><td>switch_Control         </td><td>0.07403490</td><td> 7</td><td>0.015864622</td><td> 4.6666667</td></tr>
</tbody>
</table>



The table above indicates proportion of direction switching events associated with each factor across all genes. The proportion switch is what is what was measured and the expected switch is found as the means of 100 permutations of redistributing the direction of significant expression among the contrasts that inform on each factor for each gene. This reshuffling is limited to the factor contrasts where the gene already had significant expression, reducing the extent of randomization of the redistributions. The significance ratio indicates the fold change of the actual switching relative to the expected, indicating where if the observes change was biologically meaningful.

There are limitations associated with the information gained from direction switching due to the confounding. Of the atomic factors, the phytochemical factor had the closest number of switches to what was expected while the conditioned factor had less and the acetaminophen had significantly more than expected. While this is highly dependent on the systemic approach to defining what contrasts inform on what factors, the results are consistent with downstream analysis. 

_____
# IV. Generation of Feature Matrix
_____

##### [contents](#contents)

## 9. Standardization of Effect Size Matrix 
Standard z-score definition used to normalize the values of the effect size data. 

$$
  z(x_{i}) = \frac{x_{i} - \mu}{\sigma}
$$

- $\mu$ is the per gene mean
- $\sigma$ is the per gene standard deviation

#### Before and after of z-score transformation


<table class="dataframe">
<caption>A data.frame: 5 × 5</caption>
<thead>
	<tr><th></th><th scope=col>CO_4hr_APAP_24hr_vs_Control</th><th scope=col>CO_24hr_vs_Control</th><th scope=col>Control_24hr_vs_Control</th><th scope=col>CO_4hr_APAP_4hr_vs_Control</th><th scope=col>CO_8hr_APAP_8hr_vs_Control</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003</th><td> 0.9718274</td><td> 0.7770204</td><td> 1.4152657</td><td> 0.70568538</td><td>-0.009971626</td></tr>
	<tr><th scope=row>ENSG00000001167</th><td>-2.6757986</td><td> 1.0579827</td><td> 0.1807754</td><td> 1.91842513</td><td> 0.796376095</td></tr>
	<tr><th scope=row>ENSG00000001460</th><td> 3.3602216</td><td> 0.5098319</td><td>-1.1367617</td><td> 0.00447294</td><td> 0.700851942</td></tr>
	<tr><th scope=row>ENSG00000001461</th><td>-4.2693866</td><td> 0.1707683</td><td> 0.8599573</td><td>-1.87896157</td><td>-0.699149394</td></tr>
	<tr><th scope=row>ENSG00000001497</th><td>-4.3205230</td><td>-0.5098104</td><td>-0.6231656</td><td>-0.64657533</td><td>-0.586729007</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A matrix: 5 × 5 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>CO_4hr_APAP_24hr_vs_Control</th><th scope=col>CO_24hr_vs_Control</th><th scope=col>Control_24hr_vs_Control</th><th scope=col>CO_4hr_APAP_4hr_vs_Control</th><th scope=col>CO_8hr_APAP_8hr_vs_Control</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003</th><td> 1.197897</td><td> 0.9255121</td><td> 1.8179241</td><td> 0.8257695</td><td>-0.1748817</td></tr>
	<tr><th scope=row>ENSG00000001167</th><td>-1.143254</td><td> 0.8177606</td><td> 0.3570435</td><td> 1.2696726</td><td> 0.6803624</td></tr>
	<tr><th scope=row>ENSG00000001460</th><td> 1.479748</td><td>-0.3642708</td><td>-1.4295111</td><td>-0.6912056</td><td>-0.2406932</td></tr>
	<tr><th scope=row>ENSG00000001461</th><td>-1.590014</td><td> 0.7466161</td><td> 1.1093017</td><td>-0.3320538</td><td> 0.2888221</td></tr>
	<tr><th scope=row>ENSG00000001497</th><td>-1.915255</td><td> 0.3294956</td><td> 0.2627222</td><td> 0.2489324</td><td> 0.2841856</td></tr>
</tbody>
</table>



## 10. Feature Matrix with Effect Sizes and Atomic Factor Measures


#### Peek at feature matrix


<table class="dataframe">
<caption>A data.frame: 6 × 86</caption>
<thead>
	<tr><th></th><th scope=col>CO_4hr_APAP_24hr_vs_Control</th><th scope=col>CO_24hr_vs_Control</th><th scope=col>Control_24hr_vs_Control</th><th scope=col>CO_4hr_APAP_4hr_vs_Control</th><th scope=col>CO_8hr_APAP_8hr_vs_Control</th><th scope=col>CO_4hr_vs_Control</th><th scope=col>APAP_4hr_vs_Control</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_Control_24hr</th><th scope=col>CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr</th><th scope=col>⋯</th><th scope=col>switch_APAP</th><th scope=col>switch_CONDITIONED</th><th scope=col>switch_CO_4hr_APAP_24hr</th><th scope=col>switch_CO_24hr</th><th scope=col>switch_Control_24hr</th><th scope=col>switch_CO_4hr_APAP_4hr</th><th scope=col>switch_CO_8hr_APAP_8hr</th><th scope=col>switch_CO_4hr</th><th scope=col>switch_APAP_4hr</th><th scope=col>switch_Control</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003</th><td> 0.9718274</td><td> 0.7770204</td><td> 1.4152657</td><td> 0.70568538</td><td>-0.009971626</td><td> 0.2596122</td><td> 0.285636505</td><td> 0.194807</td><td>-0.4434383</td><td> 0.2661420</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001167</th><td>-2.6757986</td><td> 1.0579827</td><td> 0.1807754</td><td> 1.91842513</td><td> 0.796376095</td><td> 1.1181610</td><td> 0.581552462</td><td>-3.733781</td><td>-2.8565740</td><td>-4.5942237</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001460</th><td> 3.3602216</td><td> 0.5098319</td><td>-1.1367617</td><td> 0.00447294</td><td> 0.700851942</td><td> 1.2548651</td><td> 0.055667483</td><td> 2.850390</td><td> 4.4969833</td><td> 3.3557487</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001461</th><td>-4.2693866</td><td> 0.1707683</td><td> 0.8599573</td><td>-1.87896157</td><td>-0.699149394</td><td>-0.6235310</td><td>-0.494020428</td><td>-4.440155</td><td>-5.1293440</td><td>-2.3904251</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001497</th><td>-4.3205230</td><td>-0.5098104</td><td>-0.6231656</td><td>-0.64657533</td><td>-0.586729007</td><td>-0.7131603</td><td>-0.007529912</td><td>-3.810713</td><td>-3.6973575</td><td>-3.6739477</td><td>⋯</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ENSG00000001617</th><td> 1.0305880</td><td>-0.3980008</td><td>-0.5706939</td><td> 0.51066571</td><td>-0.853496482</td><td>-0.0308911</td><td> 0.323684410</td><td> 1.428589</td><td> 1.6012819</td><td> 0.5199223</td><td>⋯</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
</tbody>
</table>



_____
# V. Clustering Analysis of Feature Space and Gene Enrichment Analysis
_____

##### [contents](#contents)

## 11. Hierarchical Clustering Models

#### distance based clustering dendrogram


    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_74_0.png)
    


#### correlation based clustering dendrogram


    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_76_0.png)
    


##### Dendrograms show a balance of merger events in the clustering architecture. The even density in the subtrees of the correlation proximity clustering results dendrogram shows a more heterogeneous merger profile. This means that tree cutting to establish individual clusters will generally produce clusters with strongly related genes with respect to their response patterns while not producing clusters with disparate clustering sizes. The dendrogram for the distance proximity clustering results is slightly more heterogeneous with respect to the merger density in the subtrees. This is not inherently an indication of bad clustering but the tendency is for gene clusters to have a relatively stable range. This is particularly true when working with a smaller gene set like that produced from differential expression analysis. 

## 12. Hclust Modeling Evaluation

#### Bar Plots of cluster membership size

     ..cutHeight not given, setting it to 0.536  ===>  99% of the (truncated) height range in dendro.
     ..done.



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_80_1.png)
    


     ..cutHeight not given, setting it to 1.93  ===>  99% of the (truncated) height range in dendro.
     ..done.



    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_81_1.png)
    


##### Between the two different matrices that were clustered on, the distance matrix has a better face value distribution than the correlation matrix. There is a much more stable distribution without uniformity. Also the distribution for the correlation matrix has many more cluster sizes near the maximum size which transcends what is often a limitation of the results of hierarchical clustering.

## 13. Cohesion and Separation Analysis of Distance and Correlation Based Clustering


<table class="dataframe">
<caption>A data.frame: 6 × 9</caption>
<thead>
	<tr><th scope=col>method</th><th scope=col>min</th><th scope=col>q1</th><th scope=col>median</th><th scope=col>q3</th><th scope=col>max</th><th scope=col>mean</th><th scope=col>sd</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>corr_cluster_corr</td><td> 0.52176100</td><td> 0.71264907</td><td>0.788514942</td><td>0.87357782</td><td>0.9727855</td><td>0.787983455</td><td>0.1043378</td><td>  95</td></tr>
	<tr><td>dist_cluster_corr</td><td> 0.31037405</td><td> 0.58168673</td><td>0.694331547</td><td>0.79363508</td><td>0.9384319</td><td>0.682860970</td><td>0.1437380</td><td>  83</td></tr>
	<tr><td>corr_cluster_var </td><td> 0.04686435</td><td> 0.12149310</td><td>0.215255845</td><td>0.32582875</td><td>0.6790033</td><td>0.239550994</td><td>0.1494043</td><td>  95</td></tr>
	<tr><td>dist_cluster_var </td><td> 0.13851950</td><td> 0.23419399</td><td>0.351193920</td><td>0.44314095</td><td>0.8882836</td><td>0.369144847</td><td>0.1650366</td><td>  83</td></tr>
	<tr><td>corr_silhouette  </td><td>-0.46800920</td><td>-0.08700259</td><td>0.009555757</td><td>0.09843747</td><td>0.4376006</td><td>0.001404587</td><td>0.1353559</td><td>3782</td></tr>
	<tr><td>dist_silhouette  </td><td>-0.48962889</td><td>-0.05349108</td><td>0.169765288</td><td>0.35977017</td><td>0.6992970</td><td>0.147340063</td><td>0.2547224</td><td>3782</td></tr>
</tbody>
</table>



##### Based on standard evaluation metrics for assessing hierarchical clustering results the correlation provided more stable clusters. Clustering on the correlation matrix yielded better cluster cohesion and similar separation between clusters. The lower variance scores shows that the clustering algorithm is grouping genes so that there is lower variance between them. It also is producing more clusters that may lead to more fidelity in gene enrichment. 

### sample of data object returned from enrichment analysis


    #
    # over-representation test
    #
    #...@organism 	 Homo sapiens 
    #...@ontology 	 BP 
    #...@keytype 	 ENSEMBL 
    #...@gene 	 chr [1:34] "ENSG00000006327" "ENSG00000011143" "ENSG00000067057" ...
    #...pvalues adjusted by 'BH' with cutoff <0.05 
    #...22 enriched terms found
    'data.frame':	22 obs. of  9 variables:
     $ ID         : chr  "GO:0022618" "GO:0071826" "GO:0002181" "GO:0042254" ...
     $ Description: chr  "protein-RNA complex assembly" "protein-RNA complex organization" "cytoplasmic translation" "ribosome biogenesis" ...
     $ GeneRatio  : chr  "6/33" "6/33" "4/33" "5/33" ...
     $ BgRatio    : chr  "266/21261" "274/21261" "178/21261" "365/21261" ...
     $ pvalue     : num  3.02e-06 3.59e-06 1.61e-04 2.32e-04 2.42e-04 ...
     $ p.adjust   : num  0.00138 0.00138 0.02907 0.02907 0.02907 ...
     $ qvalue     : num  0.00112 0.00112 0.0236 0.0236 0.0236 ...
     $ geneID     : chr  "ENSG00000130741/ENSG00000136937/ENSG00000168028/ENSG00000168066/ENSG00000183207/ENSG00000242372" "ENSG00000130741/ENSG00000136937/ENSG00000168028/ENSG00000168066/ENSG00000183207/ENSG00000242372" "ENSG00000130741/ENSG00000136937/ENSG00000161016/ENSG00000168028" "ENSG00000104626/ENSG00000132341/ENSG00000145912/ENSG00000168028/ENSG00000242372" ...
     $ Count      : int  6 6 4 5 2 2 2 4 2 3 ...
    #...Citation
     T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
     clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
     The Innovation. 2021, 2(3):100141 



#### A data frame was generated containing structured results of Gene Enrichment Analysis testing. Tests were performed using Ensemble identifiers that mapped to Gene Ontology's Biological Process annotations.  In all, 43 clusters were enriched.Further analysis must be performed to determine if the annotations are cohesive within each cluster. 

_____
# VI. Static Visualizations for Navigating Characteristics of Gene Clusters
_____

##### [contents](#contents)

## 14. Aggregated Factor x Feature Matrix and Heatmap



<table class="dataframe">
<caption>A matrix: 11 × 5 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>significance_consensus</th><th scope=col>upreg_consensus</th><th scope=col>downreg_consensus </th><th scope=col>sign_balance</th><th scope=col>direction_switching</th></tr>
</thead>
<tbody>
	<tr><th scope=row>CO</th><td>0.4804336</td><td>0.2427287</td><td>0.2377049</td><td> 0.09439450</td><td>0.4553146</td></tr>
	<tr><th scope=row>APAP</th><td>0.8884188</td><td>0.4323109</td><td>0.4561079</td><td>-0.04926846</td><td>1.2017451</td></tr>
	<tr><th scope=row>CONDITIONED</th><td>2.8286621</td><td>1.3246959</td><td>1.5039662</td><td>-0.67691704</td><td>2.8937070</td></tr>
	<tr><th scope=row>CO_4hr_APAP_24hr</th><td>2.3299841</td><td>0.9849286</td><td>1.3450555</td><td>-0.49270227</td><td>0.6959281</td></tr>
	<tr><th scope=row>CO_24hr</th><td>1.2556848</td><td>0.5071391</td><td>0.7485457</td><td>-0.44509078</td><td>0.6977790</td></tr>
	<tr><th scope=row>Control_24hr</th><td>1.0954521</td><td>0.6250661</td><td>0.4703860</td><td>-0.12915565</td><td>0.6089371</td></tr>
	<tr><th scope=row>CO_4hr_APAP_4hr</th><td>0.8625066</td><td>0.4182972</td><td>0.4442094</td><td>-0.20358717</td><td>0.7144368</td></tr>
	<tr><th scope=row>CO_8hr_APAP_8hr</th><td>0.9264939</td><td>0.4740878</td><td>0.4524061</td><td>-0.11821787</td><td>0.6885246</td></tr>
	<tr><th scope=row>CO_4hr</th><td>0.7091486</td><td>0.3527234</td><td>0.3564252</td><td>-0.14347788</td><td>0.5700687</td></tr>
	<tr><th scope=row>APAP_4hr</th><td>0.6840296</td><td>0.3381809</td><td>0.3458488</td><td>-0.15195664</td><td>0.5534109</td></tr>
	<tr><th scope=row>Control</th><td>0.7072977</td><td>0.3593337</td><td>0.3479640</td><td>-0.29052530</td><td>0.5182443</td></tr>
</tbody>
</table>




    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_97_0.png)
    


## 15. Cluster centroid behavior in the experimental effect space as defined by contrasts. 

#### The Heatmap indicates what differential expression pattern defines this cluster across pairwise contrasts of composite experimental states. It encodes the direction and magnitude of the within-cluster effect size of expression (aggregated as the cluster mean) respective to the individual contrasts in the experiment.


    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_101_0.png)
    


## 16. Cluster centroid behavior in the experimental effect space as defined by atomic factors. 

#### Indicates the aggregated effect size of each cluster across the factors informed on by experimental contrasts. 


    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_104_0.png)
    


## 17. Discrete Bivariate Plot: Gene Cluster vs Gene Ontology Biological Functions They are Enriched For

#### Provides the the (maximum) top five Annotations each cluster was enriched for of the 41 of clusters that were found to be enriched. The the clustering analysis successfully grouped 43% of the DEG set into distinct functionally enriched clusters. Furthermore, the largely consistent logical continuity of the top annotations found for each cluster indicates that the clustering results are characterized by strong biological signal further validating the feature space construction utilized in this pipeline. The discovered annotations are utilized in the final stages of this pipeline where functional profiles of the atomic factors in the experiment are constructed.


    
![png](deg_secondaryAnlysis_expandedContrasts_files/deg_secondaryAnlysis_expandedContrasts_106_0.png)
    


_____
# VII. Factor-Function Association Hypothesis Testing  
_____

##### [contents](#contents)

#### Hypothesis testing for the association between biological functions identified via cluster-based GO enrichment and experimental factors is performed using the camera() function from **Bioconductor's limma library**. Gene-level statistics are constructed from the difference between up- and down-regulation consensus scores for each factor derived for the experimental treatments (including interaction-derived factors), yielding a signed summary statistic per gene.

#### camera() tests whether genes annotated to a given biological function exhibit a systematic shift in this gene-level statistic relative to all other genes in the analysis universe. This corresponds to a competitive gene set test comparing the mean statistic of genes inside the gene set against genes outside the set.

#### The method utilized by camera() accounts for inter-gene correlation within gene sets when estimating the variance of the test statistic, reducing inflation of significance due to co-expression structure.

## 18. Table of Biological Functions Associated with Acetaminophen

#### Includes the results of the hypothesis testing. The results shown are above the chosen $p = 0.01$ threshold. It is important not to over interpret the p-values, but the 99% significance threshold is meant to ensure a clean association table.  


<table class="dataframe">
<caption>A data.frame: 55 × 4</caption>
<thead>
	<tr><th></th><th scope=col>NGenes</th><th scope=col>Direction</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>cytoplasmic translation</th><td>61</td><td>Up  </td><td>1.348321e-13</td><td>2.143830e-11</td></tr>
	<tr><th scope=row>ribosome biogenesis</th><td>31</td><td>Up  </td><td>2.293456e-09</td><td>1.823297e-07</td></tr>
	<tr><th scope=row>ribosomal small subunit biogenesis</th><td>15</td><td>Up  </td><td>1.091673e-08</td><td>5.785868e-07</td></tr>
	<tr><th scope=row>rRNA processing</th><td>16</td><td>Up  </td><td>1.080292e-07</td><td>4.294160e-06</td></tr>
	<tr><th scope=row>establishment of organelle localization</th><td>11</td><td>Down</td><td>2.685820e-07</td><td>8.540908e-06</td></tr>
	<tr><th scope=row>embryonic organ development</th><td> 6</td><td>Up  </td><td>1.195478e-06</td><td>2.244957e-05</td></tr>
	<tr><th scope=row>RNA splicing, via transesterification reactions with bulged adenosine as nucleophile</th><td>19</td><td>Down</td><td>1.270731e-06</td><td>2.244957e-05</td></tr>
	<tr><th scope=row>mRNA splicing, via spliceosome</th><td>19</td><td>Down</td><td>1.270731e-06</td><td>2.244957e-05</td></tr>
	<tr><th scope=row>RNA splicing, via transesterification reactions</th><td>19</td><td>Down</td><td>1.270731e-06</td><td>2.244957e-05</td></tr>
	<tr><th scope=row>rRNA metabolic process</th><td>20</td><td>Up  </td><td>2.741145e-06</td><td>4.358420e-05</td></tr>
	<tr><th scope=row>non-membrane-bounded organelle assembly</th><td>12</td><td>Up  </td><td>5.155473e-06</td><td>7.452002e-05</td></tr>
	<tr><th scope=row>endoplasmic reticulum unfolded protein response</th><td> 6</td><td>Up  </td><td>6.423580e-06</td><td>7.856533e-05</td></tr>
	<tr><th scope=row>cellular response to unfolded protein</th><td> 6</td><td>Up  </td><td>6.423580e-06</td><td>7.856533e-05</td></tr>
	<tr><th scope=row>embryonic placenta development</th><td> 8</td><td>Up  </td><td>7.631357e-06</td><td>8.089238e-05</td></tr>
	<tr><th scope=row>placenta development</th><td> 8</td><td>Up  </td><td>7.631357e-06</td><td>8.089238e-05</td></tr>
	<tr><th scope=row>ribosome assembly</th><td>13</td><td>Up  </td><td>1.268233e-05</td><td>1.260307e-04</td></tr>
	<tr><th scope=row>protein-containing complex disassembly</th><td> 6</td><td>Up  </td><td>3.076370e-05</td><td>2.574436e-04</td></tr>
	<tr><th scope=row>integrated stress response signaling</th><td> 6</td><td>Up  </td><td>3.076370e-05</td><td>2.574436e-04</td></tr>
	<tr><th scope=row>response to starvation</th><td> 6</td><td>Up  </td><td>3.076370e-05</td><td>2.574436e-04</td></tr>
	<tr><th scope=row>regulation of mRNA splicing, via spliceosome</th><td> 6</td><td>Down</td><td>4.711705e-05</td><td>3.121504e-04</td></tr>
	<tr><th scope=row>regulation of RNA splicing</th><td> 6</td><td>Down</td><td>4.711705e-05</td><td>3.121504e-04</td></tr>
	<tr><th scope=row>protein localization to plasma membrane</th><td> 6</td><td>Down</td><td>4.711705e-05</td><td>3.121504e-04</td></tr>
	<tr><th scope=row>protein localization to cell periphery</th><td> 6</td><td>Down</td><td>4.711705e-05</td><td>3.121504e-04</td></tr>
	<tr><th scope=row>chromosome segregation</th><td> 6</td><td>Down</td><td>4.711705e-05</td><td>3.121504e-04</td></tr>
	<tr><th scope=row>in utero embryonic development</th><td>10</td><td>Up  </td><td>6.936306e-05</td><td>4.411491e-04</td></tr>
	<tr><th scope=row>regulation of signal transduction by p53 class mediator</th><td> 6</td><td>Up  </td><td>1.314391e-04</td><td>7.463860e-04</td></tr>
	<tr><th scope=row>signal transduction by p53 class mediator</th><td> 6</td><td>Up  </td><td>1.314391e-04</td><td>7.463860e-04</td></tr>
	<tr><th scope=row>regulation of transcription from RNA polymerase II promoter in response to stress</th><td> 6</td><td>Up  </td><td>1.314391e-04</td><td>7.463860e-04</td></tr>
	<tr><th scope=row>alternative mRNA splicing, via spliceosome</th><td> 5</td><td>Down</td><td>1.904121e-04</td><td>9.461102e-04</td></tr>
	<tr><th scope=row>negative regulation of mRNA splicing, via spliceosome</th><td> 5</td><td>Down</td><td>1.904121e-04</td><td>9.461102e-04</td></tr>
	<tr><th scope=row>negative regulation of mRNA processing</th><td> 5</td><td>Down</td><td>1.904121e-04</td><td>9.461102e-04</td></tr>
	<tr><th scope=row>negative regulation of RNA splicing</th><td> 5</td><td>Down</td><td>1.904121e-04</td><td>9.461102e-04</td></tr>
	<tr><th scope=row>intrinsic apoptotic signaling pathway</th><td> 9</td><td>Up  </td><td>2.600470e-04</td><td>1.252954e-03</td></tr>
	<tr><th scope=row>response to unfolded protein</th><td> 8</td><td>Up  </td><td>3.313066e-04</td><td>1.498628e-03</td></tr>
	<tr><th scope=row>response to topologically incorrect protein</th><td> 8</td><td>Up  </td><td>3.313066e-04</td><td>1.498628e-03</td></tr>
	<tr><th scope=row>protein-RNA complex assembly</th><td>17</td><td>Up  </td><td>3.487374e-04</td><td>1.498628e-03</td></tr>
	<tr><th scope=row>protein-RNA complex organization</th><td>17</td><td>Up  </td><td>3.487374e-04</td><td>1.498628e-03</td></tr>
	<tr><th scope=row>regulation of intrinsic apoptotic signaling pathway</th><td> 7</td><td>Up  </td><td>4.141291e-04</td><td>1.688373e-03</td></tr>
	<tr><th scope=row>regulation of DNA-templated transcription in response to stress</th><td> 7</td><td>Up  </td><td>4.141291e-04</td><td>1.688373e-03</td></tr>
	<tr><th scope=row>negative regulation of proteolysis</th><td> 6</td><td>Up  </td><td>5.015032e-04</td><td>1.944854e-03</td></tr>
	<tr><th scope=row>regulation of hemopoiesis</th><td> 6</td><td>Up  </td><td>5.015032e-04</td><td>1.944854e-03</td></tr>
	<tr><th scope=row>positive regulation of signal transduction by p53 class mediator</th><td> 5</td><td>Up  </td><td>5.750197e-04</td><td>2.031736e-03</td></tr>
	<tr><th scope=row>p38MAPK cascade</th><td> 5</td><td>Up  </td><td>5.750197e-04</td><td>2.031736e-03</td></tr>
	<tr><th scope=row>lymphocyte differentiation</th><td> 5</td><td>Up  </td><td>5.750197e-04</td><td>2.031736e-03</td></tr>
	<tr><th scope=row>cellular response to chemical stress</th><td> 5</td><td>Up  </td><td>5.750197e-04</td><td>2.031736e-03</td></tr>
	<tr><th scope=row>regulation of mRNA processing</th><td> 9</td><td>Down</td><td>1.127130e-03</td><td>3.895948e-03</td></tr>
	<tr><th scope=row>ribosomal large subunit biogenesis</th><td> 7</td><td>Up  </td><td>1.308083e-03</td><td>4.333026e-03</td></tr>
	<tr><th scope=row>fat cell differentiation</th><td> 7</td><td>Up  </td><td>1.308083e-03</td><td>4.333026e-03</td></tr>
	<tr><th scope=row>positive regulation of DNA metabolic process</th><td> 7</td><td>Down</td><td>1.879009e-03</td><td>6.097192e-03</td></tr>
	<tr><th scope=row>negative regulation of intrinsic apoptotic signaling pathway</th><td> 5</td><td>Up  </td><td>2.176594e-03</td><td>6.529781e-03</td></tr>
	<tr><th scope=row>cellular response to external stimulus</th><td> 5</td><td>Up  </td><td>2.176594e-03</td><td>6.529781e-03</td></tr>
	<tr><th scope=row>regulation of fat cell differentiation</th><td> 5</td><td>Up  </td><td>2.176594e-03</td><td>6.529781e-03</td></tr>
	<tr><th scope=row>regulation of leukocyte differentiation</th><td> 5</td><td>Up  </td><td>2.176594e-03</td><td>6.529781e-03</td></tr>
	<tr><th scope=row>positive regulation of translation</th><td> 6</td><td>Down</td><td>2.379612e-03</td><td>6.879241e-03</td></tr>
	<tr><th scope=row>positive regulation of amide metabolic process</th><td> 6</td><td>Down</td><td>2.379612e-03</td><td>6.879241e-03</td></tr>
</tbody>
</table>



The Gene Ontology Biological Function Profiling for the contrasts that informed on the Acetaminophen Experimental Factor is strongly supported by existing peer reviewed literature. The publications in the table below show that the associations generated through this computational approach performed well in establishing a biological profile consistent with known information. The functional groups that the specific biological function annotations discovered through cluster profiling are supported by the table below



| Functional Group | Your Identified GO Terms | Validating Publication Evidence | Biological Context in Hepatocytes |
|---|---|---|---|
| p53 Signal Transduction | p53 class mediator signaling, positive regulation of p53 signal transduction | [PMC5396540](https://pmc.ncbi.nlm.nih.gov/articles/PMC5396540/) | p53 acts as a protective factor by inhibiting sustained JNK activation and regulating cell survival. |
| Ribosome & Translation | ribosome biogenesis, rRNA processing, cytoplasmic translation | [Frontiers in Cell Dev Biol (2023)](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2023.1186638/full) | Stress-induced "translational reprogramming" occurs to conserve energy and manage NAPQI-induced protein damage. |
| ER Stress & UPR | response to unfolded protein, endoplasmic reticulum UPR | [PMC4913076](https://pmc.ncbi.nlm.nih.gov/articles/PMC4913076/) | Reactive NAPQI metabolites form protein adducts, triggering the Unfolded Protein Response (UPR) and ER stress. |
| Cell Death Signaling | intrinsic apoptotic signaling pathway, regulation of apoptosis | [PMC10281617](https://pmc.ncbi.nlm.nih.gov/articles/PMC10281617/) | APAP activates apoptotic pathways and senescence (p21/Cdkn1a). |
| Metabolic Adaptive Response | cellular response to chemical stress, response to starvation | [PMC5687251](https://pmc.ncbi.nlm.nih.gov/articles/PMC5687251/) | Reflects the exhaustion of glutathione (GSH) and the subsequent metabolic shift to oxidative stress management. |
| Developmental Reprogramming | embryonic organ development, placenta development | [PMC10499594](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499594/) | Hepatocytes transiently upregulate fetal-specific genes as part of a pro-regenerative response after acute injury. |



The Placenta/Embryonic terms are a common artifact in the hepatic regeneration studies where damaged hepatocytes often re-express fetal-specific genes [doi: 10.3390/biology14101361](https://doi.org/10.3390/biology14101361). 

## 19. Tables of Biological Functions Associated with the Phytochemical Treatment

#### Includes the results of the hypothesis testing. The results shown are above the chosen $p = 0.01$ threshold. It is important not to over interpret the p-values; the 99% significance threshold is meant to ensure a clean association table. There are additional subsets of this table showing where biological functions associated with the phytochemical factor are also shared with the Acetaminophen Factor and where they are distinct from the Acetaminophen Factor.


<table class="dataframe">
<caption>A data.frame: 49 × 4</caption>
<thead>
	<tr><th></th><th scope=col>NGenes</th><th scope=col>Direction</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>peptidyl-lysine modification</th><td> 5</td><td>Down</td><td>3.178580e-10</td><td>3.041272e-08</td></tr>
	<tr><th scope=row>lymphocyte differentiation</th><td> 5</td><td>Up  </td><td>3.825500e-10</td><td>3.041272e-08</td></tr>
	<tr><th scope=row>cellular oxidant detoxification</th><td> 9</td><td>Up  </td><td>6.453402e-08</td><td>2.052182e-06</td></tr>
	<tr><th scope=row>cellular detoxification</th><td> 9</td><td>Up  </td><td>6.453402e-08</td><td>2.052182e-06</td></tr>
	<tr><th scope=row>cellular response to toxic substance</th><td> 9</td><td>Up  </td><td>6.453402e-08</td><td>2.052182e-06</td></tr>
	<tr><th scope=row>B cell activation</th><td> 7</td><td>Up  </td><td>1.984353e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>cell redox homeostasis</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to hypoxia</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to decreased oxygen levels</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>homeostasis of number of cells</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to oxygen levels</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to oxidative stress</th><td>19</td><td>Up  </td><td>1.216394e-05</td><td>1.611722e-04</td></tr>
	<tr><th scope=row>cellular respiration</th><td> 7</td><td>Up  </td><td>1.898846e-05</td><td>2.156546e-04</td></tr>
	<tr><th scope=row>energy derivation by oxidation of organic compounds</th><td> 7</td><td>Up  </td><td>1.898846e-05</td><td>2.156546e-04</td></tr>
	<tr><th scope=row>regulation of DNA-binding transcription factor activity</th><td> 6</td><td>Down</td><td>3.236542e-05</td><td>3.430735e-04</td></tr>
	<tr><th scope=row>regulation of apoptotic signaling pathway</th><td> 6</td><td>Up  </td><td>3.709422e-05</td><td>3.469401e-04</td></tr>
	<tr><th scope=row>positive regulation of cellular catabolic process</th><td> 6</td><td>Up  </td><td>3.709422e-05</td><td>3.469401e-04</td></tr>
	<tr><th scope=row>myeloid leukocyte differentiation</th><td> 5</td><td>Up  </td><td>7.114817e-05</td><td>5.953979e-04</td></tr>
	<tr><th scope=row>regulation of leukocyte differentiation</th><td> 5</td><td>Up  </td><td>7.114817e-05</td><td>5.953979e-04</td></tr>
	<tr><th scope=row>proton transmembrane transport</th><td>10</td><td>Up  </td><td>8.959243e-05</td><td>7.122599e-04</td></tr>
	<tr><th scope=row>axonogenesis</th><td> 7</td><td>Down</td><td>1.271228e-04</td><td>9.187513e-04</td></tr>
	<tr><th scope=row>regulation of neuron projection development</th><td> 7</td><td>Down</td><td>1.271228e-04</td><td>9.187513e-04</td></tr>
	<tr><th scope=row>regulation of axonogenesis</th><td> 6</td><td>Down</td><td>2.752830e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>positive regulation of translation</th><td> 6</td><td>Up  </td><td>3.106405e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>positive regulation of amide metabolic process</th><td> 6</td><td>Up  </td><td>3.106405e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>regulation of hemopoiesis</th><td> 6</td><td>Up  </td><td>3.106405e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>histone modification</th><td>10</td><td>Down</td><td>3.676196e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>aerobic electron transport chain</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>ATP synthesis coupled electron transport</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>mitochondrial ATP synthesis coupled electron transport</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>respiratory electron transport chain</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>electron transport chain</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>regulation of fat cell differentiation</th><td> 5</td><td>Up  </td><td>6.692888e-04</td><td>3.040483e-03</td></tr>
	<tr><th scope=row>negative regulation of translation</th><td> 5</td><td>Up  </td><td>6.692888e-04</td><td>3.040483e-03</td></tr>
	<tr><th scope=row>negative regulation of amide metabolic process</th><td> 5</td><td>Up  </td><td>6.692888e-04</td><td>3.040483e-03</td></tr>
	<tr><th scope=row>oxidative phosphorylation</th><td>11</td><td>Up  </td><td>8.323541e-04</td><td>3.576873e-03</td></tr>
	<tr><th scope=row>aerobic respiration</th><td>11</td><td>Up  </td><td>8.323541e-04</td><td>3.576873e-03</td></tr>
	<tr><th scope=row>RNA localization</th><td> 7</td><td>Up  </td><td>8.981559e-04</td><td>3.661712e-03</td></tr>
	<tr><th scope=row>fat cell differentiation</th><td> 7</td><td>Up  </td><td>8.981559e-04</td><td>3.661712e-03</td></tr>
	<tr><th scope=row>ATP biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>purine ribonucleoside triphosphate biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>purine nucleoside triphosphate biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>ribonucleoside triphosphate biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>peptide hormone secretion</th><td> 6</td><td>Down</td><td>1.812798e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>peptide secretion</th><td> 6</td><td>Down</td><td>1.812798e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>regulation of transcription from RNA polymerase II promoter in response to stress</th><td> 6</td><td>Up  </td><td>2.014582e-03</td><td>6.963447e-03</td></tr>
	<tr><th scope=row>purine ribonucleotide metabolic process</th><td>16</td><td>Up  </td><td>2.634669e-03</td><td>8.549232e-03</td></tr>
	<tr><th scope=row>ribonucleotide metabolic process</th><td>16</td><td>Up  </td><td>2.634669e-03</td><td>8.549232e-03</td></tr>
	<tr><th scope=row>ribose phosphate metabolic process</th><td>16</td><td>Up  </td><td>2.634669e-03</td><td>8.549232e-03</td></tr>
</tbody>
</table>



#### Phytochemical associated functions also associated with Acetaminophen


<table class="dataframe">
<caption>A data.frame: 8 × 4</caption>
<thead>
	<tr><th></th><th scope=col>NGenes</th><th scope=col>Direction</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>lymphocyte differentiation</th><td>5</td><td>Up</td><td>3.825500e-10</td><td>3.041272e-08</td></tr>
	<tr><th scope=row>regulation of leukocyte differentiation</th><td>5</td><td>Up</td><td>7.114817e-05</td><td>5.953979e-04</td></tr>
	<tr><th scope=row>positive regulation of translation</th><td>6</td><td>Up</td><td>3.106405e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>positive regulation of amide metabolic process</th><td>6</td><td>Up</td><td>3.106405e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>regulation of hemopoiesis</th><td>6</td><td>Up</td><td>3.106405e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>regulation of fat cell differentiation</th><td>5</td><td>Up</td><td>6.692888e-04</td><td>3.040483e-03</td></tr>
	<tr><th scope=row>fat cell differentiation</th><td>7</td><td>Up</td><td>8.981559e-04</td><td>3.661712e-03</td></tr>
	<tr><th scope=row>regulation of transcription from RNA polymerase II promoter in response to stress</th><td>6</td><td>Up</td><td>2.014582e-03</td><td>6.963447e-03</td></tr>
</tbody>
</table>



#### Phytochemical associated functions not associated with Acetaminophen


<table class="dataframe">
<caption>A data.frame: 41 × 4</caption>
<thead>
	<tr><th></th><th scope=col>NGenes</th><th scope=col>Direction</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>peptidyl-lysine modification</th><td> 5</td><td>Down</td><td>3.178580e-10</td><td>3.041272e-08</td></tr>
	<tr><th scope=row>cellular oxidant detoxification</th><td> 9</td><td>Up  </td><td>6.453402e-08</td><td>2.052182e-06</td></tr>
	<tr><th scope=row>cellular detoxification</th><td> 9</td><td>Up  </td><td>6.453402e-08</td><td>2.052182e-06</td></tr>
	<tr><th scope=row>cellular response to toxic substance</th><td> 9</td><td>Up  </td><td>6.453402e-08</td><td>2.052182e-06</td></tr>
	<tr><th scope=row>B cell activation</th><td> 7</td><td>Up  </td><td>1.984353e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>cell redox homeostasis</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to hypoxia</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to decreased oxygen levels</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>homeostasis of number of cells</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to oxygen levels</th><td> 6</td><td>Up  </td><td>3.417840e-06</td><td>4.940332e-05</td></tr>
	<tr><th scope=row>response to oxidative stress</th><td>19</td><td>Up  </td><td>1.216394e-05</td><td>1.611722e-04</td></tr>
	<tr><th scope=row>cellular respiration</th><td> 7</td><td>Up  </td><td>1.898846e-05</td><td>2.156546e-04</td></tr>
	<tr><th scope=row>energy derivation by oxidation of organic compounds</th><td> 7</td><td>Up  </td><td>1.898846e-05</td><td>2.156546e-04</td></tr>
	<tr><th scope=row>regulation of DNA-binding transcription factor activity</th><td> 6</td><td>Down</td><td>3.236542e-05</td><td>3.430735e-04</td></tr>
	<tr><th scope=row>regulation of apoptotic signaling pathway</th><td> 6</td><td>Up  </td><td>3.709422e-05</td><td>3.469401e-04</td></tr>
	<tr><th scope=row>positive regulation of cellular catabolic process</th><td> 6</td><td>Up  </td><td>3.709422e-05</td><td>3.469401e-04</td></tr>
	<tr><th scope=row>myeloid leukocyte differentiation</th><td> 5</td><td>Up  </td><td>7.114817e-05</td><td>5.953979e-04</td></tr>
	<tr><th scope=row>proton transmembrane transport</th><td>10</td><td>Up  </td><td>8.959243e-05</td><td>7.122599e-04</td></tr>
	<tr><th scope=row>axonogenesis</th><td> 7</td><td>Down</td><td>1.271228e-04</td><td>9.187513e-04</td></tr>
	<tr><th scope=row>regulation of neuron projection development</th><td> 7</td><td>Down</td><td>1.271228e-04</td><td>9.187513e-04</td></tr>
	<tr><th scope=row>regulation of axonogenesis</th><td> 6</td><td>Down</td><td>2.752830e-04</td><td>1.899686e-03</td></tr>
	<tr><th scope=row>histone modification</th><td>10</td><td>Down</td><td>3.676196e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>aerobic electron transport chain</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>ATP synthesis coupled electron transport</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>mitochondrial ATP synthesis coupled electron transport</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>respiratory electron transport chain</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>electron transport chain</th><td> 8</td><td>Up  </td><td>4.089487e-04</td><td>2.031964e-03</td></tr>
	<tr><th scope=row>negative regulation of translation</th><td> 5</td><td>Up  </td><td>6.692888e-04</td><td>3.040483e-03</td></tr>
	<tr><th scope=row>negative regulation of amide metabolic process</th><td> 5</td><td>Up  </td><td>6.692888e-04</td><td>3.040483e-03</td></tr>
	<tr><th scope=row>oxidative phosphorylation</th><td>11</td><td>Up  </td><td>8.323541e-04</td><td>3.576873e-03</td></tr>
	<tr><th scope=row>aerobic respiration</th><td>11</td><td>Up  </td><td>8.323541e-04</td><td>3.576873e-03</td></tr>
	<tr><th scope=row>RNA localization</th><td> 7</td><td>Up  </td><td>8.981559e-04</td><td>3.661712e-03</td></tr>
	<tr><th scope=row>ATP biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>purine ribonucleoside triphosphate biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>purine nucleoside triphosphate biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>ribonucleoside triphosphate biosynthetic process</th><td>15</td><td>Up  </td><td>1.803210e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>peptide hormone secretion</th><td> 6</td><td>Down</td><td>1.812798e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>peptide secretion</th><td> 6</td><td>Down</td><td>1.812798e-03</td><td>6.405218e-03</td></tr>
	<tr><th scope=row>purine ribonucleotide metabolic process</th><td>16</td><td>Up  </td><td>2.634669e-03</td><td>8.549232e-03</td></tr>
	<tr><th scope=row>ribonucleotide metabolic process</th><td>16</td><td>Up  </td><td>2.634669e-03</td><td>8.549232e-03</td></tr>
	<tr><th scope=row>ribose phosphate metabolic process</th><td>16</td><td>Up  </td><td>2.634669e-03</td><td>8.549232e-03</td></tr>
</tbody>
</table>



#### The table below is a functional grouping assessment of the annotations associated with the phytochemical factor that were not also associated with the acetaminophen factor. Overall it profiles cellular response to inhibition of energy production machinery and reduction/oxidation balance.

 Grouping Profile | Key GO Terms | Biological Interpretation |
| :--- | :--- | :--- |
| **Mitochondrial Bioenergetics** | oxidative phosphorylation, ATP biosynthetic process, aerobic respiration, electron transport chain | This is the most dominant signal. The treatment significantly impacts how the hepatocyte generates energy, specifically via the mitochondria. |
| **Redox & Detoxification** | cellular oxidant detoxification, cell redox homeostasis, cellular response to toxic substance | Indicates a response to chemical-driven oxidative species or a direct modulation of the cell's antioxidant capacity. |
| **Hypoxia & Oxygen Sensing** | response to hypoxia, response to oxygen levels, response to decreased oxygen levels | This suggests the treatment might be mimicking or inducing a low-oxygen state, potentially through the stabilization of HIF (Hypoxia-Inducible Factor) pathways. |
| **Immune/Hematopoietic Signaling** | lymphocyte differentiation, B cell activation, regulation of hemopoiesis | In a pure hepatocyte line, this "leaky" signaling often indicates the activation of inflammatory cytokines or pathways (like NF-κB or IL-6) that are shared with immune cell development. |
| **Epigenetic/Post-Translational** | histone modification, peptidyl-lysine modification, RNA localization | This indicates the treatment is inducing structural changes at the chromatin level or modifying protein function after translation. |







#### The collective biological grouping profile, based on the associated Gene Ontology terms, suggests that the phytochemical treatment most significantly provokes mitochondrial activity which may or may not be related to potential low oxygen state induced by the phytochemical. 

## 20. Tables of Biological Functions Associated with the Compositional Treatment

#### Includes the results of the hypothesis testing. The results shown are above the chosen $p = 0.01$ threshold. It is important not to over interpret the p-values, but the 99% significance threshold is meant to ensure a clean association table. This table provides the means for peripheral level assessment of interactivity based on whether or not functional associations are shown to either be more significantly associated with the compositional treatment or where new functional association is introduced. 


<table class="dataframe">
<caption>A data.frame: 59 × 4</caption>
<thead>
	<tr><th></th><th scope=col>NGenes</th><th scope=col>Direction</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>cytoplasmic translation</th><td>61</td><td>Up  </td><td>2.875679e-14</td><td>4.572330e-12</td></tr>
	<tr><th scope=row>ribosome biogenesis</th><td>31</td><td>Up  </td><td>4.103466e-10</td><td>3.262255e-08</td></tr>
	<tr><th scope=row>ribosomal small subunit biogenesis</th><td>15</td><td>Up  </td><td>1.731157e-08</td><td>9.175133e-07</td></tr>
	<tr><th scope=row>establishment of organelle localization</th><td>11</td><td>Down</td><td>5.665943e-08</td><td>2.252212e-06</td></tr>
	<tr><th scope=row>response to unfolded protein</th><td> 8</td><td>Up  </td><td>9.168406e-08</td><td>2.429628e-06</td></tr>
	<tr><th scope=row>response to topologically incorrect protein</th><td> 8</td><td>Up  </td><td>9.168406e-08</td><td>2.429628e-06</td></tr>
	<tr><th scope=row>endoplasmic reticulum unfolded protein response</th><td> 6</td><td>Up  </td><td>6.766992e-07</td><td>1.276235e-05</td></tr>
	<tr><th scope=row>cellular response to unfolded protein</th><td> 6</td><td>Up  </td><td>6.766992e-07</td><td>1.276235e-05</td></tr>
	<tr><th scope=row>RNA splicing, via transesterification reactions with bulged adenosine as nucleophile</th><td>19</td><td>Down</td><td>8.829300e-07</td><td>1.276235e-05</td></tr>
	<tr><th scope=row>mRNA splicing, via spliceosome</th><td>19</td><td>Down</td><td>8.829300e-07</td><td>1.276235e-05</td></tr>
	<tr><th scope=row>RNA splicing, via transesterification reactions</th><td>19</td><td>Down</td><td>8.829300e-07</td><td>1.276235e-05</td></tr>
	<tr><th scope=row>ribosome assembly</th><td>13</td><td>Up  </td><td>1.061625e-06</td><td>1.406653e-05</td></tr>
	<tr><th scope=row>response to starvation</th><td> 6</td><td>Up  </td><td>3.062973e-06</td><td>3.746251e-05</td></tr>
	<tr><th scope=row>in utero embryonic development</th><td>10</td><td>Up  </td><td>5.472817e-06</td><td>6.215556e-05</td></tr>
	<tr><th scope=row>embryonic placenta development</th><td> 8</td><td>Up  </td><td>9.124330e-06</td><td>9.067303e-05</td></tr>
	<tr><th scope=row>placenta development</th><td> 8</td><td>Up  </td><td>9.124330e-06</td><td>9.067303e-05</td></tr>
	<tr><th scope=row>rRNA processing</th><td>16</td><td>Up  </td><td>1.138387e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>rRNA metabolic process</th><td>20</td><td>Up  </td><td>1.174015e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>regulation of transcription from RNA polymerase II promoter in response to stress</th><td> 6</td><td>Up  </td><td>1.269175e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>embryonic organ development</th><td> 6</td><td>Up  </td><td>1.269175e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>chromosome segregation</th><td> 6</td><td>Down</td><td>1.408170e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>protein-RNA complex assembly</th><td>17</td><td>Up  </td><td>1.421221e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>protein-RNA complex organization</th><td>17</td><td>Up  </td><td>1.421221e-05</td><td>9.824959e-05</td></tr>
	<tr><th scope=row>regulation of DNA-templated transcription in response to stress</th><td> 7</td><td>Up  </td><td>1.528496e-05</td><td>1.012628e-04</td></tr>
	<tr><th scope=row>regulation of mRNA splicing, via spliceosome</th><td> 6</td><td>Down</td><td>2.764375e-05</td><td>1.690521e-04</td></tr>
	<tr><th scope=row>regulation of RNA splicing</th><td> 6</td><td>Down</td><td>2.764375e-05</td><td>1.690521e-04</td></tr>
	<tr><th scope=row>protein-containing complex disassembly</th><td> 6</td><td>Up  </td><td>4.817108e-05</td><td>2.836742e-04</td></tr>
	<tr><th scope=row>protein localization to plasma membrane</th><td> 6</td><td>Down</td><td>5.309440e-05</td><td>2.911038e-04</td></tr>
	<tr><th scope=row>protein localization to cell periphery</th><td> 6</td><td>Down</td><td>5.309440e-05</td><td>2.911038e-04</td></tr>
	<tr><th scope=row>cellular response to external stimulus</th><td> 5</td><td>Up  </td><td>7.901046e-05</td><td>4.008904e-04</td></tr>
	<tr><th scope=row>intrinsic apoptotic signaling pathway</th><td> 9</td><td>Up  </td><td>8.044543e-05</td><td>4.008904e-04</td></tr>
	<tr><th scope=row>non-membrane-bounded organelle assembly</th><td>12</td><td>Up  </td><td>8.068235e-05</td><td>4.008904e-04</td></tr>
	<tr><th scope=row>regulation of signal transduction by p53 class mediator</th><td> 6</td><td>Up  </td><td>9.082624e-05</td><td>4.247462e-04</td></tr>
	<tr><th scope=row>signal transduction by p53 class mediator</th><td> 6</td><td>Up  </td><td>9.082624e-05</td><td>4.247462e-04</td></tr>
	<tr><th scope=row>histone modification</th><td>10</td><td>Down</td><td>2.357940e-04</td><td>1.071178e-03</td></tr>
	<tr><th scope=row>negative regulation of mRNA splicing, via spliceosome</th><td> 5</td><td>Down</td><td>2.615651e-04</td><td>1.094443e-03</td></tr>
	<tr><th scope=row>negative regulation of mRNA processing</th><td> 5</td><td>Down</td><td>2.615651e-04</td><td>1.094443e-03</td></tr>
	<tr><th scope=row>negative regulation of RNA splicing</th><td> 5</td><td>Down</td><td>2.615651e-04</td><td>1.094443e-03</td></tr>
	<tr><th scope=row>regulation of intrinsic apoptotic signaling pathway</th><td> 7</td><td>Up  </td><td>4.782021e-04</td><td>1.932178e-03</td></tr>
	<tr><th scope=row>ribosomal large subunit assembly</th><td> 6</td><td>Up  </td><td>5.346909e-04</td><td>1.932178e-03</td></tr>
	<tr><th scope=row>response to hypoxia</th><td> 6</td><td>Up  </td><td>5.346909e-04</td><td>1.932178e-03</td></tr>
	<tr><th scope=row>response to decreased oxygen levels</th><td> 6</td><td>Up  </td><td>5.346909e-04</td><td>1.932178e-03</td></tr>
	<tr><th scope=row>homeostasis of number of cells</th><td> 6</td><td>Up  </td><td>5.346909e-04</td><td>1.932178e-03</td></tr>
	<tr><th scope=row>response to oxygen levels</th><td> 6</td><td>Up  </td><td>5.346909e-04</td><td>1.932178e-03</td></tr>
	<tr><th scope=row>positive regulation of signal transduction by p53 class mediator</th><td> 5</td><td>Up  </td><td>5.586911e-04</td><td>1.974042e-03</td></tr>
	<tr><th scope=row>ribosomal large subunit biogenesis</th><td> 7</td><td>Up  </td><td>7.961040e-04</td><td>2.751751e-03</td></tr>
	<tr><th scope=row>alternative mRNA splicing, via spliceosome</th><td> 5</td><td>Down</td><td>9.008405e-04</td><td>3.047524e-03</td></tr>
	<tr><th scope=row>integrated stress response signaling</th><td> 6</td><td>Up  </td><td>9.248565e-04</td><td>3.063587e-03</td></tr>
	<tr><th scope=row>response to oxidative stress</th><td>19</td><td>Up  </td><td>9.638239e-04</td><td>3.127510e-03</td></tr>
	<tr><th scope=row>lymphocyte differentiation</th><td> 5</td><td>Up  </td><td>1.017985e-03</td><td>3.237191e-03</td></tr>
	<tr><th scope=row>fat cell differentiation</th><td> 7</td><td>Up  </td><td>1.301532e-03</td><td>4.057716e-03</td></tr>
	<tr><th scope=row>muscle cell development</th><td> 8</td><td>Down</td><td>1.561110e-03</td><td>4.697983e-03</td></tr>
	<tr><th scope=row>negative regulation of proteolysis</th><td> 6</td><td>Up  </td><td>1.565994e-03</td><td>4.697983e-03</td></tr>
	<tr><th scope=row>negative regulation of intrinsic apoptotic signaling pathway</th><td> 5</td><td>Up  </td><td>1.807643e-03</td><td>5.322503e-03</td></tr>
	<tr><th scope=row>positive regulation of cellular catabolic process</th><td> 6</td><td>Up  </td><td>2.595965e-03</td><td>7.369857e-03</td></tr>
	<tr><th scope=row>regulation of hemopoiesis</th><td> 6</td><td>Up  </td><td>2.595965e-03</td><td>7.369857e-03</td></tr>
	<tr><th scope=row>positive regulation of DNA metabolic process</th><td> 7</td><td>Down</td><td>2.642024e-03</td><td>7.369857e-03</td></tr>
	<tr><th scope=row>regulation of ubiquitin protein ligase activity</th><td> 5</td><td>Up  </td><td>3.128609e-03</td><td>8.431336e-03</td></tr>
	<tr><th scope=row>myeloid leukocyte differentiation</th><td> 5</td><td>Up  </td><td>3.128609e-03</td><td>8.431336e-03</td></tr>
</tbody>
</table>



### Table Comparing Functional Profiles in the Conditioned Factor and the Independent Treatments
| Functional Group | Status in This Profile | Comparison to Previous Profiles |
| :--- | :--- | :--- |
| **Proteotoxicity (UPR/Ribosome)** | **High** | Identical to the Tylenol profile |
| **p53/Apoptosis** | **High** | Identical to the Tylenol profile |
| **Oxygen/Hypoxia Sensing** | **Present** | Matches the phytochemical profile |
| **Mitochondrial ATP/OXPHOS** | **Absent** | Distinctly missing compared to phytochemical |
| **Ubiquitin/Muscle** | **New Signaling** | Unique to this specific list |


### Table of contrasts that inform on each of the atomic factors


<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th></th><th scope=col>factor</th><th scope=col>contrasts</th><th scope=col>length</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>CO         </td><td>CO_24hr_vs_Control, CO_4hr_vs_Control, CO_24hr_vs_Control_24hr, CO_24hr_vs_CO_4hr, CO_4hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_APAP_4hr                                                                                                                                                                                                                                                                                                                                              </td><td> 6</td></tr>
	<tr><th scope=row>2</th><td>APAP       </td><td>APAP_4hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_4hr_APAP_4hr, CO_4hr_APAP_24hr_vs_CO_4hr, APAP_4hr_vs_Control_24hr, CO_4hr_APAP_4hr_vs_CO_4hr                                                                                                                                                                                                                                                                                                                                           </td><td> 5</td></tr>
	<tr><th scope=row>3</th><td>CONDITIONED</td><td>CO_4hr_APAP_24hr_vs_Control, CO_4hr_APAP_4hr_vs_Control, CO_8hr_APAP_8hr_vs_Control, CO_4hr_APAP_24hr_vs_CO_24hr, CO_4hr_APAP_24hr_vs_Control_24hr, CO_4hr_APAP_24hr_vs_CO_8hr_APAP_8hr, CO_4hr_APAP_24hr_vs_APAP_4hr, CO_24hr_vs_CO_4hr_APAP_4hr, CO_24hr_vs_CO_8hr_APAP_8hr, CO_24hr_vs_APAP_4hr, CO_4hr_APAP_4hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_Control_24hr, CO_8hr_APAP_8hr_vs_CO_4hr_APAP_4hr, CO_8hr_APAP_8hr_vs_CO_4hr, CO_8hr_APAP_8hr_vs_APAP_4hr, APAP_4hr_vs_CO_4hr</td><td>16</td></tr>
</tbody>
</table>



#### Confounding in the available samples included the effect that time played in the expression profiles generated for each experimental group. Time was not separable as a factor because intra-experimental group time changes where obfuscated by identifying where the phytochemical treatment was primed by acetaminophen exposure, and inter-experimental group changes were largely implicit in the separation of how contrasts informed on the effects of the phytochemical and acetaminophen independently. The absence of the Mitochondrial ATP activity and Oxidative Phosphorylation in the "conditioned" factor profile could be emergent in the interactive treatments or it could represent the overall time exposure differences to the phytochemical in the contrasts that inform on the "conditioned" factor and the phytochemical factor. If not time dependent, this suggests that the induced mitochondrial response associated with the phytochemical treatment may be mitigating known negative effects of tylenol toxicity in the co treatments primed by acetaminophen. Similarly, Ubiquitin can be the result of longer exposures to acetaminophen that are found within the "Conditioned" profile or it can indicate cellular response to the two treatments where protein folding intervention is no longer manageable due to the treatments rather than time. A merged table with full results of the association testing is available as a [csv](../de_analysis/treatment_factor_annotation.csv).


<table class="dataframe">
<caption>A data.frame: 159 × 13</caption>
<thead>
	<tr><th scope=col>BF_annotation</th><th scope=col>apap_NGenes</th><th scope=col>co_NGenes</th><th scope=col>conditioned_NGenes</th><th scope=col>apap_Direction</th><th scope=col>co_Direction</th><th scope=col>conditioned_Direction</th><th scope=col>apap_PValue</th><th scope=col>co_PValue</th><th scope=col>conditioned_PValue</th><th scope=col>apap_FDR</th><th scope=col>co_FDR</th><th scope=col>conditioned_FDR</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>actin filament bundle assembly            </td><td> 5</td><td> 5</td><td> 5</td><td>Down</td><td>Down</td><td>Down</td><td>6.549671e-02</td><td>2.198029e-02</td><td>4.746132e-03</td><td>1.285676e-01</td><td>5.065023e-02</td><td>1.199021e-02</td></tr>
	<tr><td>actin filament bundle organization        </td><td> 5</td><td> 5</td><td> 5</td><td>Down</td><td>Down</td><td>Down</td><td>6.549671e-02</td><td>2.198029e-02</td><td>4.746132e-03</td><td>1.285676e-01</td><td>5.065023e-02</td><td>1.199021e-02</td></tr>
	<tr><td>actomyosin structure organization         </td><td> 5</td><td> 5</td><td> 5</td><td>Down</td><td>Down</td><td>Down</td><td>6.549671e-02</td><td>4.230925e-03</td><td>4.746132e-03</td><td>1.285676e-01</td><td>1.226672e-02</td><td>1.199021e-02</td></tr>
	<tr><td>adaptive thermogenesis                    </td><td> 9</td><td> 9</td><td> 9</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.499075e-01</td><td>6.634341e-01</td><td>6.625744e-01</td><td>2.628573e-01</td><td>8.177211e-01</td><td>7.216301e-01</td></tr>
	<tr><td>aerobic electron transport chain          </td><td> 8</td><td> 8</td><td> 8</td><td>Up  </td><td>Up  </td><td>Up  </td><td>5.193426e-01</td><td>4.089487e-04</td><td>8.536034e-01</td><td>5.955055e-01</td><td>2.031964e-03</td><td>8.870780e-01</td></tr>
	<tr><td>aerobic respiration                       </td><td>11</td><td>11</td><td>11</td><td>Up  </td><td>Up  </td><td>Up  </td><td>2.914080e-01</td><td>8.323541e-04</td><td>6.671674e-01</td><td>4.371120e-01</td><td>3.576873e-03</td><td>7.216301e-01</td></tr>
	<tr><td>alternative mRNA splicing, via spliceosome</td><td> 5</td><td> 5</td><td> 5</td><td>Down</td><td>Down</td><td>Down</td><td>1.904121e-04</td><td>9.885987e-01</td><td>9.008405e-04</td><td>9.461102e-04</td><td>9.885987e-01</td><td>3.047524e-03</td></tr>
	<tr><td>ATP biosynthetic process                  </td><td>15</td><td>15</td><td>15</td><td>Up  </td><td>Up  </td><td>Up  </td><td>7.773790e-01</td><td>1.803210e-03</td><td>8.776679e-01</td><td>8.078644e-01</td><td>6.405218e-03</td><td>8.888484e-01</td></tr>
	<tr><td>ATP metabolic process                     </td><td>12</td><td>12</td><td>12</td><td>Up  </td><td>Up  </td><td>Up  </td><td>4.381517e-01</td><td>7.872743e-02</td><td>5.255168e-01</td><td>5.574774e-01</td><td>1.345985e-01</td><td>6.235610e-01</td></tr>
	<tr><td>ATP synthesis coupled electron transport  </td><td> 8</td><td> 8</td><td> 8</td><td>Up  </td><td>Up  </td><td>Up  </td><td>5.193426e-01</td><td>4.089487e-04</td><td>8.536034e-01</td><td>5.955055e-01</td><td>2.031964e-03</td><td>8.870780e-01</td></tr>
	<tr><td>axonogenesis                              </td><td> 7</td><td> 7</td><td> 7</td><td>Up  </td><td>Down</td><td>Down</td><td>3.171006e-01</td><td>1.271228e-04</td><td>4.272818e-01</td><td>4.422719e-01</td><td>9.187513e-04</td><td>5.856708e-01</td></tr>
	<tr><td>B cell activation                         </td><td> 7</td><td> 7</td><td> 7</td><td>Up  </td><td>Up  </td><td>Up  </td><td>9.845977e-03</td><td>1.984353e-06</td><td>9.033262e-02</td><td>2.525017e-02</td><td>4.940332e-05</td><td>1.650906e-01</td></tr>
	<tr><td>blood coagulation                         </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Up  </td><td>Down</td><td>6.952477e-01</td><td>3.085194e-01</td><td>2.347450e-01</td><td>7.469215e-01</td><td>4.341114e-01</td><td>3.393132e-01</td></tr>
	<tr><td>canonical NF-kappaB signal transduction   </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Up  </td><td>Up  </td><td>9.641966e-01</td><td>9.050754e-02</td><td>1.925449e-01</td><td>9.641966e-01</td><td>1.483577e-01</td><td>3.092389e-01</td></tr>
	<tr><td>cell redox homeostasis                    </td><td> 6</td><td> 6</td><td> 6</td><td>Down</td><td>Up  </td><td>Up  </td><td>5.243445e-01</td><td>3.417840e-06</td><td>3.476864e-02</td><td>5.955055e-01</td><td>4.940332e-05</td><td>6.503781e-02</td></tr>
	<tr><td>cellular detoxification                   </td><td> 9</td><td> 9</td><td> 9</td><td>Down</td><td>Up  </td><td>Up  </td><td>4.417746e-01</td><td>6.453402e-08</td><td>2.453699e-02</td><td>5.574774e-01</td><td>2.052182e-06</td><td>4.700460e-02</td></tr>
	<tr><td>cellular oxidant detoxification           </td><td> 9</td><td> 9</td><td> 9</td><td>Down</td><td>Up  </td><td>Up  </td><td>4.417746e-01</td><td>6.453402e-08</td><td>2.453699e-02</td><td>5.574774e-01</td><td>2.052182e-06</td><td>4.700460e-02</td></tr>
	<tr><td>cellular respiration                      </td><td> 7</td><td> 7</td><td> 7</td><td>Up  </td><td>Up  </td><td>Up  </td><td>5.133805e-02</td><td>1.898846e-05</td><td>2.009238e-01</td><td>1.074046e-01</td><td>2.156546e-04</td><td>3.134984e-01</td></tr>
	<tr><td>cellular response to chemical stress      </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Down</td><td>Up  </td><td>5.750197e-04</td><td>9.885987e-01</td><td>5.278709e-03</td><td>2.031736e-03</td><td>9.885987e-01</td><td>1.199021e-02</td></tr>
	<tr><td>cellular response to external stimulus    </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Up  </td><td>Up  </td><td>2.176594e-03</td><td>2.611864e-01</td><td>7.901046e-05</td><td>6.529781e-03</td><td>3.917796e-01</td><td>4.008904e-04</td></tr>
	<tr><td>cellular response to oxidative stress     </td><td> 9</td><td> 9</td><td> 9</td><td>Up  </td><td>Up  </td><td>Up  </td><td>2.445655e-01</td><td>3.761972e-03</td><td>6.288544e-03</td><td>3.812344e-01</td><td>1.196307e-02</td><td>1.408280e-02</td></tr>
	<tr><td>cellular response to toxic substance      </td><td> 9</td><td> 9</td><td> 9</td><td>Down</td><td>Up  </td><td>Up  </td><td>4.417746e-01</td><td>6.453402e-08</td><td>2.453699e-02</td><td>5.574774e-01</td><td>2.052182e-06</td><td>4.700460e-02</td></tr>
	<tr><td>cellular response to tumor necrosis factor</td><td> 8</td><td> 8</td><td> 8</td><td>Up  </td><td>Up  </td><td>Down</td><td>3.479518e-01</td><td>8.166096e-03</td><td>8.403352e-01</td><td>4.810812e-01</td><td>1.997553e-02</td><td>8.870780e-01</td></tr>
	<tr><td>cellular response to unfolded protein     </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Up  </td><td>Up  </td><td>6.423580e-06</td><td>6.161270e-01</td><td>6.766992e-07</td><td>7.856533e-05</td><td>7.837136e-01</td><td>1.276235e-05</td></tr>
	<tr><td>chromosome segregation                    </td><td> 6</td><td> 6</td><td> 6</td><td>Down</td><td>Down</td><td>Down</td><td>4.711705e-05</td><td>9.875684e-01</td><td>1.408170e-05</td><td>3.121504e-04</td><td>9.885987e-01</td><td>9.824959e-05</td></tr>
	<tr><td>circadian regulation of gene expression   </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Down</td><td>Down</td><td>4.626061e-01</td><td>3.716782e-02</td><td>2.347450e-01</td><td>5.746435e-01</td><td>7.674913e-02</td><td>3.393132e-01</td></tr>
	<tr><td>coagulation                               </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Up  </td><td>Down</td><td>6.952477e-01</td><td>3.085194e-01</td><td>2.347450e-01</td><td>7.469215e-01</td><td>4.341114e-01</td><td>3.393132e-01</td></tr>
	<tr><td>cold-induced thermogenesis                </td><td> 9</td><td> 9</td><td> 9</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.499075e-01</td><td>6.634341e-01</td><td>6.625744e-01</td><td>2.628573e-01</td><td>8.177211e-01</td><td>7.216301e-01</td></tr>
	<tr><td>cytoplasmic translation                   </td><td>61</td><td>61</td><td>61</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.348321e-13</td><td>9.675632e-01</td><td>2.875679e-14</td><td>2.143830e-11</td><td>9.885987e-01</td><td>4.572330e-12</td></tr>
	<tr><td>detoxification                            </td><td>12</td><td>12</td><td>12</td><td>Down</td><td>Up  </td><td>Down</td><td>6.863268e-02</td><td>1.608926e-01</td><td>9.301972e-01</td><td>1.330804e-01</td><td>2.532863e-01</td><td>9.360845e-01</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>response to lipopolysaccharide                                                      </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Up  </td><td>Up  </td><td>7.207851e-03</td><td>2.369143e-02</td><td>8.683965e-03</td><td>1.878768e-02</td><td>5.305546e-02</td><td>1.865879e-02</td></tr>
	<tr><td>response to molecule of bacterial origin                                            </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Up  </td><td>Up  </td><td>7.207851e-03</td><td>2.369143e-02</td><td>8.683965e-03</td><td>1.878768e-02</td><td>5.305546e-02</td><td>1.865879e-02</td></tr>
	<tr><td>response to oxidative stress                                                        </td><td>19</td><td>19</td><td>19</td><td>Up  </td><td>Up  </td><td>Up  </td><td>4.175169e-01</td><td>1.216394e-05</td><td>9.638239e-04</td><td>5.574774e-01</td><td>1.611722e-04</td><td>3.127510e-03</td></tr>
	<tr><td>response to oxygen levels                                                           </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Up  </td><td>Up  </td><td>1.555114e-01</td><td>3.417840e-06</td><td>5.346909e-04</td><td>2.628573e-01</td><td>4.940332e-05</td><td>1.932178e-03</td></tr>
	<tr><td>response to starvation                                                              </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Up  </td><td>Up  </td><td>3.076370e-05</td><td>6.161270e-01</td><td>3.062973e-06</td><td>2.574436e-04</td><td>7.837136e-01</td><td>3.746251e-05</td></tr>
	<tr><td>response to topologically incorrect protein                                         </td><td> 8</td><td> 8</td><td> 8</td><td>Up  </td><td>Up  </td><td>Up  </td><td>3.313066e-04</td><td>3.846244e-01</td><td>9.168406e-08</td><td>1.498628e-03</td><td>5.272007e-01</td><td>2.429628e-06</td></tr>
	<tr><td>response to toxic substance                                                         </td><td> 8</td><td> 8</td><td> 8</td><td>Down</td><td>Down</td><td>Down</td><td>1.570531e-01</td><td>2.531013e-02</td><td>2.183596e-01</td><td>2.628573e-01</td><td>5.589320e-02</td><td>3.370794e-01</td></tr>
	<tr><td>response to tumor necrosis factor                                                   </td><td>15</td><td>15</td><td>15</td><td>Up  </td><td>Up  </td><td>Down</td><td>4.842505e-01</td><td>2.951225e-02</td><td>5.651952e-01</td><td>5.955055e-01</td><td>6.341145e-02</td><td>6.548072e-01</td></tr>
	<tr><td>response to unfolded protein                                                        </td><td> 8</td><td> 8</td><td> 8</td><td>Up  </td><td>Up  </td><td>Up  </td><td>3.313066e-04</td><td>3.846244e-01</td><td>9.168406e-08</td><td>1.498628e-03</td><td>5.272007e-01</td><td>2.429628e-06</td></tr>
	<tr><td>ribonucleoside triphosphate biosynthetic process                                    </td><td>15</td><td>15</td><td>15</td><td>Up  </td><td>Up  </td><td>Up  </td><td>7.773790e-01</td><td>1.803210e-03</td><td>8.776679e-01</td><td>8.078644e-01</td><td>6.405218e-03</td><td>8.888484e-01</td></tr>
	<tr><td>ribonucleoside triphosphate metabolic process                                       </td><td>12</td><td>12</td><td>12</td><td>Up  </td><td>Up  </td><td>Up  </td><td>4.381517e-01</td><td>7.872743e-02</td><td>5.255168e-01</td><td>5.574774e-01</td><td>1.345985e-01</td><td>6.235610e-01</td></tr>
	<tr><td>ribonucleotide biosynthetic process                                                 </td><td>12</td><td>12</td><td>12</td><td>Up  </td><td>Up  </td><td>Up  </td><td>3.118586e-01</td><td>7.872743e-02</td><td>5.255168e-01</td><td>4.388099e-01</td><td>1.345985e-01</td><td>6.235610e-01</td></tr>
	<tr><td>ribonucleotide metabolic process                                                    </td><td>16</td><td>16</td><td>16</td><td>Up  </td><td>Up  </td><td>Up  </td><td>1.998225e-01</td><td>2.634669e-03</td><td>2.969761e-01</td><td>3.177178e-01</td><td>8.549232e-03</td><td>4.178691e-01</td></tr>
	<tr><td>ribose phosphate biosynthetic process                                               </td><td>12</td><td>12</td><td>12</td><td>Up  </td><td>Up  </td><td>Up  </td><td>3.118586e-01</td><td>7.872743e-02</td><td>5.255168e-01</td><td>4.388099e-01</td><td>1.345985e-01</td><td>6.235610e-01</td></tr>
	<tr><td>ribose phosphate metabolic process                                                  </td><td>16</td><td>16</td><td>16</td><td>Up  </td><td>Up  </td><td>Up  </td><td>1.998225e-01</td><td>2.634669e-03</td><td>2.969761e-01</td><td>3.177178e-01</td><td>8.549232e-03</td><td>4.178691e-01</td></tr>
	<tr><td>ribosomal large subunit assembly                                                    </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Down</td><td>Up  </td><td>5.224041e-03</td><td>9.875684e-01</td><td>5.346909e-04</td><td>1.457232e-02</td><td>9.885987e-01</td><td>1.932178e-03</td></tr>
	<tr><td>ribosomal large subunit biogenesis                                                  </td><td> 7</td><td> 7</td><td> 7</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.308083e-03</td><td>9.866339e-01</td><td>7.961040e-04</td><td>4.333026e-03</td><td>9.885987e-01</td><td>2.751751e-03</td></tr>
	<tr><td>ribosomal small subunit biogenesis                                                  </td><td>15</td><td>15</td><td>15</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.091673e-08</td><td>5.143400e-01</td><td>1.731157e-08</td><td>5.785868e-07</td><td>6.930513e-01</td><td>9.175133e-07</td></tr>
	<tr><td>ribosome assembly                                                                   </td><td>13</td><td>13</td><td>13</td><td>Up  </td><td>Up  </td><td>Up  </td><td>1.268233e-05</td><td>7.503610e-01</td><td>1.061625e-06</td><td>1.260307e-04</td><td>9.107435e-01</td><td>1.406653e-05</td></tr>
	<tr><td>ribosome biogenesis                                                                 </td><td>31</td><td>31</td><td>31</td><td>Up  </td><td>Up  </td><td>Up  </td><td>2.293456e-09</td><td>7.050327e-01</td><td>4.103466e-10</td><td>1.823297e-07</td><td>8.623093e-01</td><td>3.262255e-08</td></tr>
	<tr><td>RNA catabolic process                                                               </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Up  </td><td>Up  </td><td>2.092579e-02</td><td>4.628949e-03</td><td>5.278709e-03</td><td>4.753144e-02</td><td>1.226672e-02</td><td>1.199021e-02</td></tr>
	<tr><td>RNA localization                                                                    </td><td> 7</td><td> 7</td><td> 7</td><td>Up  </td><td>Up  </td><td>Up  </td><td>4.936774e-01</td><td>8.981559e-04</td><td>1.707207e-02</td><td>5.955055e-01</td><td>3.661712e-03</td><td>3.525272e-02</td></tr>
	<tr><td>RNA splicing, via transesterification reactions                                     </td><td>19</td><td>19</td><td>19</td><td>Down</td><td>Down</td><td>Down</td><td>1.270731e-06</td><td>9.790902e-01</td><td>8.829300e-07</td><td>2.244957e-05</td><td>9.885987e-01</td><td>1.276235e-05</td></tr>
	<tr><td>RNA splicing, via transesterification reactions with bulged adenosine as nucleophile</td><td>19</td><td>19</td><td>19</td><td>Down</td><td>Down</td><td>Down</td><td>1.270731e-06</td><td>9.790902e-01</td><td>8.829300e-07</td><td>2.244957e-05</td><td>9.885987e-01</td><td>1.276235e-05</td></tr>
	<tr><td>rRNA metabolic process                                                              </td><td>20</td><td>20</td><td>20</td><td>Up  </td><td>Down</td><td>Up  </td><td>2.741145e-06</td><td>5.755497e-01</td><td>1.174015e-05</td><td>4.358420e-05</td><td>7.690118e-01</td><td>9.824959e-05</td></tr>
	<tr><td>rRNA processing                                                                     </td><td>16</td><td>16</td><td>16</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.080292e-07</td><td>2.163502e-01</td><td>1.138387e-05</td><td>4.294160e-06</td><td>3.339775e-01</td><td>9.824959e-05</td></tr>
	<tr><td>signal transduction by p53 class mediator                                           </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.314391e-04</td><td>9.875684e-01</td><td>9.082624e-05</td><td>7.463860e-04</td><td>9.885987e-01</td><td>4.247462e-04</td></tr>
	<tr><td>temperature homeostasis                                                             </td><td> 9</td><td> 9</td><td> 9</td><td>Up  </td><td>Down</td><td>Up  </td><td>1.499075e-01</td><td>6.634341e-01</td><td>6.625744e-01</td><td>2.628573e-01</td><td>8.177211e-01</td><td>7.216301e-01</td></tr>
	<tr><td>toll-like receptor signaling pathway                                                </td><td> 5</td><td> 5</td><td> 5</td><td>Up  </td><td>Up  </td><td>Up  </td><td>5.339689e-02</td><td>9.050754e-02</td><td>5.202680e-01</td><td>1.102611e-01</td><td>1.483577e-01</td><td>6.235610e-01</td></tr>
	<tr><td>wound healing                                                                       </td><td> 6</td><td> 6</td><td> 6</td><td>Up  </td><td>Up  </td><td>Down</td><td>6.952477e-01</td><td>3.085194e-01</td><td>2.347450e-01</td><td>7.469215e-01</td><td>4.341114e-01</td><td>3.393132e-01</td></tr>
</tbody>
</table>



_____
# IX. Conclusion
_____

#### While this analysis proved effective for creating a strong profile of phytochemical and acetaminophen treatments, clarifications for the interactive effects of these treatments remain undetermined. There was clear evidence of the independent effects each treatment had on hepatic expression and this evidence was validated in the conditioned profile that captures the same functions from the two independent treatments as well as the consistency of the acetaminophen treatment profile with known effects on hepatocyte gene expression. From this, the efficacy of the analysis pipeline is supported by the results and it is very likely that the profile for the phytochemical is characterized by strong associations with stimulation of mitochondrial activity. Based on this association, there may be potential for therapeutic resources related to contexts where conditions interfere with mitochondrial activity. 

#### **Future steps could include:**
1. additional samples for current experiment:
    - clarification of suggested interactive differential expression by establishing additional rna sequencing samples
      - Multiple samples were filtered out of this analysis because they were shown to have count profiles not indicative of a common treatment. Establishing additional samples would allow for testing designs in DESeq2 that meet full rank requirements for the design matrix upon which the differential analysis is dependent
    - Control groups at each time exposure intervals would allow for clearer separation of the impact the vehicle itself has on expression levels

    - Additional samples can be generated independently for analysis along with the viable samples that informed on this analysis provided the same environment conditions can be established for the new samples. 
  
1. Future experimentation can be devised to explicitly focus on the potential mitigating effects that the induced mitochondrial activity associated with phytochemical treatment might have on conditions where this activity is insufficient with analysis performed to establish whether or not the environment that induced the response reflects a supportive environment for cell/tissue viability.

1. More comprehensive clustering and enrichment analysis can be performed using a gene expression set not filtered by differential expression analysis results. Since the results of analysis restricted to the DEG set picked up biological signal, it can serve as validation of clustering results for fine tuning parameters such as those required for WGCNA. Additional information can likely be determined for the cell functions associated with each aspect of the treatment protocol.
