1. clean_counts.R:
    - supports clean_counts snakemake rule
    - creates nTPMs and updates column names in data extracted from counts.txt so that that DESeq metadata ('coldata') can be parsed from them
    - outputs: 
        - "updated_sampleNames_map.csv" 
            - contains parallel columns with pairings indicating what new name replaced the old name so that it can be demonstrated the identified of samples were not swapped
        - "de_analysis/counts_matrix.csv" 
            - for DESeq2 analysis
        - "de_analysis/tpm_counts.csv" 
            - for presenting counts in normalized convention

2. DEanalysis.R:
    - generates test results for experimental group contrasts and exports them for secondary analysis
    - visuals
        - per test sample dist matrix heatmap and sample pca plot
            - check group replicate proximity and separation
        - volcano plots for each contrast 
        -