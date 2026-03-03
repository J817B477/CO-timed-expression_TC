suppressPackageStartupMessages({

  library(DESeq2)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
})

##############################################################
#---- read in counts df and formats counts matrix object ----#
##############################################################

counts_df = read.csv(snakemake@input[[1]])

# create raw counts matrix
counts_matrix = counts_df[, !grepl(".nTPM", colnames(counts_df))]

rownames(counts_matrix) = counts_matrix$Geneid
counts_matrix$Geneid = NULL

# subsets to only count-value columns 
counts_matrix = counts_matrix[,-c(1:5)]



##################################################
#---- generate sample metadata (full samples)----#
##################################################

samples = colnames(counts_matrix)
# Extract CO treatment time if present
co_time = str_extract(samples, "(?<=CO_)[0-9]+hr")
co_time[is.na(co_time)] = "0hr"

# Extract APAP treatment time if present
apap_time = str_extract(samples, "(?<=APAP_)[0-9]+hr")
apap_time[is.na(apap_time)] = "0hr"


treatment = sapply(samples, function(x){sub("_.$","",x)})
# Extract replicate (assume last number in name)
replicate = str_extract(samples, "[0-9]+$")

# construct colData object containing metadata
time_factor_levels = c("0hr","4hr","8hr","24hr")
coldata = data.frame(
  row.names = samples,
  treatment = factor(treatment,
                     levels = c("Control", 
                                setdiff(unique(treatment),"Control"))),
  # as.numeric treats duration as a factor rather than category
  co_time = factor(co_time, levels = time_factor_levels), 
  apap_time = factor(apap_time, levels = time_factor_levels),
  replicate = factor(replicate)
)

#~~~~~~~~~~~~~~~~~~~~~~~Experimental Group Analysis~~~~~~~~~~~~~~~~~~~~~~~#

#######################################################
#---- Sample comparison by grouping (all samples) ----#
#######################################################

# construct dds object: design assumes treatments are dependent
dds_treatment = DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~ treatment)


# filter for low expression levels
## filters dds to remove sparse genes, which are defined here to be
## the set of genes that do not have 10 or more counts in at least n
## samples, where n = the typical minimum group size.
filtered_dds_treatment = dds_treatment[rowSums(counts(dds_treatment)>=10)>=3,]

## normalize data for heteroskedasticity
rld = rlog(filtered_dds_treatment, blind = TRUE) # doesn't account for treatment groups

## create sample distance matrix
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)

## create hierarchical heatmap of sample distances
pdf("de_analysis/full-sample_dist_heatmap.pdf", width = 7, height = 5)  # width/height in inches

colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         main = "Sample Heatmap: All Samples",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

dev.off()

## create pca plot for sample groupings
pcaData = plotPCA(rld, intgroup = "treatment", returnData = TRUE)

percentVar = round(100 * attr(pcaData, "percentVar"))

pdf("de_analysis/full-sample_dist_PCAplot.pdf", width = 7, height = 5)

ggplot(pcaData, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = rownames(pcaData)), size = 3) + 
  xlim(min(pcaData$PC1) - 5, max(pcaData$PC1) + 5) +  # adjust these margins
  ylim(min(pcaData$PC2) - 2, max(pcaData$PC2) + 2) +
  labs(title = "PCA sample plot: all samples",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       caption = "<add caption here>") +
  coord_fixed(ratio = 2) +  # Optional: ensures 1 unit on x == 1 unit on y
  theme_minimal() 

dev.off()

###################################################################
#---- Differential Analysis: Experimental Group (All Samples) ----#
###################################################################
de_dds_treatment = DESeq(dds_treatment)

contrast_groups = resultsNames(de_dds_treatment)[-1]
cat("\nTreatment Contrast Names:\n",contrast_groups, '\n')

res_df_list = list()
for(contrast in contrast_groups){
  res_df = as.data.frame(results(de_dds_treatment, name = contrast))

  res_df$significance = with(res_df, 
                              ifelse(padj < 0.05 & log2FoldChange > 0, "upregulated",
                              ifelse(padj < 0.05 & log2FoldChange < 0, "downregulated",
                              "not significant")))

  res_df$significance[is.na(res_df$significance)] = "not significant" 

  plot_name = paste0("de_analysis/treatment_", contrast, "_volcano_plot.pdf")

  pdf(plot_name, width = 7, height = 5)
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    labs(title = "Volcano Plot: Differentially Expressed Genes",
        subtitle = paste0("Experimental Contrast: ", contrast),
        x = expression(log[2]*(fold~change)),
        y = expression(-log[10]*(Adj.~P-value))
    ) +
    geom_point(aes(color = significance), alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("upregulated" = "red", 
                                  "downregulated" = "#145291", 
                                  "not significant" = "grey")) +
    theme_minimal()
  dev.off()

  contrast_str = sub("treatment_","",contrast)

  res_df_list[[contrast_str]] = res_df
}

res_df_treatment = do.call(cbind, res_df_list)
colnames(res_df_treatment) = paste0("groupTest_", colnames(res_df_treatment))
res_df_treatment$Geneid <- rownames(res_df_list[[1]])
#~~~~~~~~~~~~~~~~~~~~~Treatment Interaction Analysis~~~~~~~~~~~~~~~~~~~~~#


######################################################
#---- generate sample metadata (filtered samples)----#
######################################################

# subsets counts_matrix to exclude samples with 8hr exposure intervals
counts_matrix_interaction = counts_matrix[, !grepl(".8.|.24.", colnames(counts_matrix))]
# subsets coldata to exclude samples with 8hr exposure intervals
coldata_interaction = coldata[!grepl(".8.|.24.", rownames(coldata)),]

# creates dds object omitting 8hr and 24hr exposure samples
dds_interaction = DESeqDataSetFromMatrix(countData = counts_matrix_interaction,
                              colData = coldata_interaction,
                              design = ~ co_time * apap_time)



######################################################################
#---- treatment interaction sample comparison (filtered samples) ----#
######################################################################

# filter for low expression levels
## filters dds to remove sparse genes, which are defined here to be 
## the set of genes that do not have 10 or more counts in at least n 
## samples, where n = the typical minimum group size.
filtered_dds_interaction = dds_interaction[rowSums(counts(dds_interaction)>=10)>=3,]

## normalize data for heteroskedasticity
rld_interaction = rlog(filtered_dds_interaction, blind = TRUE) # doesn't account for treatment groups

## create sample distance matrix
sampleDists_interaction = dist(t(assay(rld_interaction)))
sampleDistMatrix_interaction = as.matrix(sampleDists_interaction)

## create hierarchical heatmap of sample distances
pdf("de_analysis/filtered-sample_dist_heatmap.pdf", width = 7, height = 5)  # width/height in inches

colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         main = "Sample Heatmap: Filtered Samples",
         clustering_distance_rows = sampleDists_interaction,
         clustering_distance_cols = sampleDists_interaction,
         col = colors)

dev.off()

## create pca plot for sample groupings
pcaData_interaction = plotPCA(rld_interaction, intgroup = "treatment", returnData = TRUE)

percentVar_interaction = round(100 * attr(pcaData_interaction, "percentVar"))

pdf("de_analysis/filtered-sample_dist_PCAplot.pdf", width = 7, height = 5)

ggplot(pcaData_interaction, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = rownames(pcaData_interaction)), size = 3) + 
  xlim(min(pcaData_interaction$PC1) - 5, max(pcaData_interaction$PC1) + 5) +  
  ylim(min(pcaData_interaction$PC2) - 2, max(pcaData_interaction$PC2) + 2) +
  labs(title = "PCA sample plot: filtered samples",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       caption = "<add caption here>") +
  coord_fixed(ratio = 2) +  # Optional: ensures 1 unit on x == 1 unit only
  theme_minimal() 

dev.off()

##########################################################################
#---- Differential Analysis: CO/APAP interaction (omits 8/24hr ints) ----#
##########################################################################

de_dds_interaction = DESeq(dds_interaction)

contrast_groups = resultsNames(de_dds_interaction)[-1]
cat("\nTreatment Contrast Names:\n",contrast_groups, '\n')

res_df_list = list()
for(contrast in contrast_groups){
  res_df = as.data.frame(results(de_dds_interaction, name = contrast))

  res_df$significance = with(res_df, 
                              ifelse(padj < 0.05 & log2FoldChange > 0, "upregulated",
                              ifelse(padj < 0.05 & log2FoldChange < 0, "downregulated",
                              "not significant")))

  res_df$significance[is.na(res_df$significance)] = "not significant" 

  plot_name = paste0("de_analysis/interaction_", contrast, "_volcano_plot.pdf")

  pdf(plot_name, width = 7, height = 5)
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    labs(title = "Volcano Plot: Differentially Expressed Genes",
        subtitle = paste0("Experimental Contrast: ", contrast),
        x = expression(log[2]*(fold~change)),
        y = expression(-log[10]*(Adj.~P-value))
    ) +
    geom_point(aes(color = significance), alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("upregulated" = "red", 
                                  "downregulated" = "#145291", 
                                  "not significant" = "grey")) +
    theme_minimal()
  dev.off()

  res_df_list[[contrast]] = res_df
}

res_df_interaction = do.call(cbind, res_df_list)
colnames(res_df_interaction) = paste0("interactionTest_", colnames(res_df_interaction))
res_df_interaction$Geneid <- rownames(res_df_list[[1]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Combining Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# augments counts_df with treatment group DE results data

de_counts_df =  merge(counts_df, res_df_treatment,
                   by = "Geneid",
                   all = TRUE)

de_counts_df =  merge(de_counts_df, res_df_interaction,
                   by = "Geneid",
                   all = TRUE)

print("Final Col Names:") 
print(colnames(de_counts_df))


## sort rows by padj values
padj_cols <- grep("padj", colnames(de_counts_df), value = TRUE)
de_counts_df <- de_counts_df[do.call(order, de_counts_df[padj_cols]), ]

## assigns gene symbols 
ens_ids = de_counts_df$Geneid
## strips off version extensions of ensemble ids
ens_ids_clean = sub("\\..*", "", ens_ids)

## gets gene symbols from ensemble ids
gene_symbols = AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

de_counts_df$gene_symbol = gene_symbols

# subsets to requested attributes
de_sub_df = de_counts_df[, c(
  "Geneid",
  "gene_symbol",
  grep("nTPM|padj|significance",
       colnames(de_counts_df),
       value = TRUE)
)]


# subsets to only DGEs in at least one contrast
de_sub_df = de_sub_df[apply(de_sub_df, 1, function(row) any(grepl("regulated", row))), ]

write.csv(de_counts_df,snakemake@output[[1]])
write.csv(de_sub_df,snakemake@output[[2]])
