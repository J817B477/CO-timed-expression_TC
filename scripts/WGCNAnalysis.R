suppressPackageStartupMessages({
 library(WGCNA)
 library(DESeq2)
 library(gridExtra)
 library(parallel)
 library(foreach)
 library(doParallel)
})

# allows multi-threading
allowWGCNAThreads()

# gets isolated counts data and augmented data set:
## wgcna performed on independent counts matrix
## relevant results added to augmented data
matrix_df = read.csv(snakemake@input[[1]])
full_df = read.csv(snakemake@input[[2]])

## assigns gene identifiers to rows so only samples are attributes
rownames(matrix_df) = matrix_df$Geneid
matrix_df$Geneid = NULL 


# check for sample outliers and initial gene filtering
## read in sample metadata (generated in DE analysis)
gsg = goodSamplesGenes(t(matrix_df))

print("Sample Outlier Analysis:")
summary(gsg)

cat("\nAll Samples are OK:", gsg$allOK)

cat("\nNumber of Good Genes:") 
print(table(gsg$goodGenes))
cat("\nNumber of Good Samples:")
print(table(gsg$goodSamples))

htree = hclust(dist(t(matrix_df)), method = "average")

pdf("WGCNAnalysis/sample_htree.pdf", width = 7, height = 5)
plot (htree)
dev.off()
# remove flagged genes
matrix_df = matrix_df[gsg$goodGenes == TRUE,]
cat("remaining genes:\n", nrow(matrix_df), "\n")

# removes genes with < 15 counts in >= 75% of samples
matrix_df75 = matrix_df[rowSums(matrix_df >= 15) >= floor(.75 * ncol(matrix_df)),]

cat("Number of Genes with at least 15 counts in 75% of samples: ", nrow(matrix_df75),"\n")

# Create Variance Stabilized (Normalized) counts
## read in sample metadata (generated in DE analysis)
coldata = read.csv(snakemake@input[[3]])
coldata$'Sample' = coldata$'X'
coldata$'X' = NULL

print(coldata)
dds = DESeqDataSetFromMatrix(countData = matrix_df75,
                             colData = coldata,
                             design = ~1) # 1 = no model; keeps analysis experimentally agnostic
dds_norm = vst(dds)

norm_matrix = assay(dds_norm)
print(head(norm_matrix))

# choose soft threshold
find_best_power = function(norm_counts){
    if (dim(norm_counts)[1] > dim(norm_counts)[2]){norm_counts = t(norm_counts)}

    # circumstantial filtering of df
    row_vars = apply(norm_counts, 1, var)
    threshold = quantile(row_vars, .30)

    low_var_indices = which(row_vars <= threshold)
    power_seeking_df = norm_counts[-low_var_indices,]

    # establishes range of test powers to choose from    
    powers = c(c(1:10), seq(from = 12, to = 50, by = 2))

    # deploys WGNCA package function for characterizing power candidates
    sft = pickSoftThreshold(power_seeking_df,
                      powerVector = powers,
                      networkType = "signed",
                      verbose = 5)

    # extracts sft index and power index as parallel vectors
    sft_indices = sft$fitIndices$SFT.R.sq
    powers = sft$fitIndices$Power

    # finds location of minimum sft_index that is greater than the threshold
    power_index = which(sft_indices >= 0.8)[1]
    sft_value = sft_indices[power_index]
    power = powers[power_index]

    # provides user with information regarding the best power
    cat("\nsft index value: ", sft_value, "\n")
    cat("power: ", power, "\n")
    
    return(power)
}

# parameter fine tuning

if(!file.exists("WGCNAnalysis/wgcna_parameter_grid_metrics.csv")){    
    # establish test values for parameters
    sft_power = find_best_power(norm_matrix)
    sft_power_range = seq(sft_power-6,sft_power+2,by = 2)
    dp_split_range = seq(0,4,by = 1)
    min_mod_range = seq(20,100,by=10)
    merge_ht_range = seq(.15,.30,by =.05)


    results <- data.frame(
    # param args
    beta = numeric(),
    deepSplit = integer(),
    minModuleSize = integer(),
    mergeCutHeight = numeric(),
    # metrics
    nMods = integer(),
    effD = numeric(),
    delta_effD_nMods = numeric(),
    globalMedian_mod_cohesion = numeric(),
    grey_proportion = numeric(),
    stringsAsFactors = FALSE
    )

    temp_cor = cor
    cor = WGCNA::cor

    # limit cor expansion inside each worker
    Sys.setenv(OMP_NUM_THREADS = 1)
    Sys.setenv(OPENBLAS_NUM_THREADS = 1)
    Sys.setenv(MKL_NUM_THREADS = 1)

    # parallel processing of parameter test value combinations: 1 per sft_pwr candidate
    n_cores = length(sft_power_range)
    cl = makeCluster(n_cores)
    registerDoParallel(cl)  


    results_df = foreach(beta = sft_power_range, .combine = 'rbind', .packages='WGCNA') %dopar% {
        express_df = t(norm_matrix)
        # same TOMdiss per sft_power
        adj = adjacency(express_df, power = beta, type = "signed")
        TOMdiss = 1- TOMsimilarity(adj)
        geneTree = hclust(as.dist(TOMdiss), method = "average")


        inner_results <- data.frame()

        for(ds in dp_split_range){
            for(ms in min_mod_range){
                for(mh in merge_ht_range){
                    mods = cutreeDynamic(dendro = geneTree,
                                    distM = TOMdiss,
                                    deepSplit = ds,
                                    minClusterSize = ms)
                    
                    # mods and their eigengenes
                    merged = mergeCloseModules(express_df, colors = mods, cutHeight = mh)
                    mergedMods = merged$colors
                    mergedMEs = merged$newMEs

                    # generate metrics
                    n_mods = length(unique(mergedMods)) - ifelse(any(mergedMods == "grey"),1,0)

                    eig_vals <- eigen(cor(mergedMEs))$values
                    effD = sum(eig_vals)^2 / sum(eig_vals^2)

                    delta_effD_nMods = n_mods - effD

                    mod_member_table = table(merged$colors)


                    grey_proportion = if("grey" %in% names(mod_member_table)){
                        mod_member_table["grey"]/sum(mod_member_table)
                    } else {0}

                    # df of gene's correlation with each eigen gene in results
                    kME_all = signedKME(express_df, mergedMEs)
                    print(kME_all)

                    #creates named list
                    mod_names = setdiff(unique(mergedMods),"grey")
                    mod_k_list = numeric(length(mod_names))
                    names(mod_k_list) = mod_names

                    # populates list with mod-level cohesion
                    for(name in mod_names){
                        mod_genes_idx = which(mergedMods == name)
                        kME_colname = paste0("kME", name)
                        mod_k_list[name] = median(abs(kME_all[mod_genes_idx, kME_colname]))
                    }

                    globalMedian_mod_cohesion = median(mod_k_list)

                    inner_results = rbind(inner_results, data.frame(
                        beta = beta,
                        deepSplit = ds,
                        minModuleSize = ms,
                        mergeCutHeight = mh,
                        n_mods = n_mods,
                        effD = effD,
                        delta_effD_nMods = delta_effD_nMods,
                        globalMedian_mod_cohesion = globalMedian_mod_cohesion,
                        grey_proportion = grey_proportion
                    ))
                }
            }
        }
        inner_results
    }

    stopCluster(cl)

    print(head(results_df))

    write.csv(results_df, "WGCNAnalysis/wgcna_parameter_grid_metrics.csv", row.names = FALSE)
} else {
   results_df = read.csv("WGCNAnalysis/wgcna_parameter_grid_metrics.csv")
}


#~~~~~~~~~~~~~~~~~~will update to run with best params ~~~~~~~~~~~~~~~#

# # build modules network: produces only 3 module eigengenes
# modules_net = blockwiseModules(t(norm_matrix),
#                                maxBlockSize = 14000,
#                                TOMType = "signed",
#                                power = sft_power,
#                                mergeCutHeight = 0.25,
#                                numericLabels = FALSE,
#                                randomSeed = 54,
#                                verbose = 3)

# eigengenes_df = modules_net$MEs

# cor = temp_cor
# print(head(eigengenes_df))
# print(table(modules_net$colors))

# plotDendroAndColors(dendro = modules_net$dendrograms[[1]],
#                     colors = modules_net$colors,
#                     main = "Gene dendrogram and module colors")