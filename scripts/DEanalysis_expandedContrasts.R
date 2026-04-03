#' script level docstring:
#' Script executes differential expression analysis for cross product
#' of contrasts. Test design uses glm with experimental groups as
#' predictors. Interactive test model was attempted since experimental
#' groups are sometimes exposed to both phytochemical and acetaminophen.
#' Full rank deficiency restricted use of majority of samples in this
#' restricting time exposure to 4 hours which generated very limited results.
#' (Generated coldata object reflects this attempt to perform interactive test).
#' To address this limitation, secondary analysis was performed on comprehensive
#' set of contrasts in:
#' - coding_notebooks/deg_secondaryAnlysis_expandedContrasts.ipynb

# prevent unnecessary terminal output
suppressPackageStartupMessages({

  library(DESeq2)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
})

# create directory for storing volcano plots
dir.create("de_analysis/expandedContrasts_plots")

##############################################################
#---- read in counts df and formats counts matrix object ----#
##############################################################

counts_df <- read.csv(snakemake@input[[1]])

# create raw counts matrix
counts_matrix <- counts_df[, !grepl(".nTPM", colnames(counts_df))]
tpm_df <- counts_df[, c("Geneid", grep(".nTPM", colnames(counts_df), value = TRUE))]

rownames(counts_matrix) <- counts_matrix$Geneid
counts_matrix$Geneid <- NULL

# subset to only count-value columns
counts_matrix <- counts_matrix[, -c(1:5)]


####################################
#---- generate sample metadata ----#
####################################

# generates vector of sample names for parsing
samples <- colnames(counts_matrix)

# extract CO treatment time if present
## ?<= means matches pattern only if preceded by 'CO' pattern (lookbehind)
co_time <- str_extract(samples, "(?<=CO_)[0-9]+hr")
co_time[is.na(co_time)] <- "0hr"

# extract APAP treatment time if present
apap_time <- str_extract(samples, "(?<=APAP_)[0-9]+hr")
apap_time[is.na(apap_time)] <- "0hr"

# establish treatment of each sample (removes replicate id)
treatment <- sapply(samples, function(x) {
  sub("_.$", "", x)
})

# extract replicate (assumes last number in name)
replicate <- str_extract(samples, "[0-9]+$")

# construct colData object containing metadata
time_factor_levels <- c("0hr", "4hr", "8hr", "24hr")
coldata <- data.frame(
  row.names = samples,
  treatment = factor(treatment,
                     levels = c("Control",
                                setdiff(unique(treatment), "Control"))),
  co_time = factor(co_time, levels = time_factor_levels), 
  apap_time = factor(apap_time, levels = time_factor_levels),
  replicate = factor(replicate)
)

# write coldata to csv (for documenting)
write.csv(coldata, "de_analysis/coldata.csv")


#####################################################################
#---- Differential Analysis: Experimental Group (All Contrasts) ----#
#####################################################################

# create dds object: required object for facilitating DESeq testing
dds_treatment <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                        colData = coldata,
                                        design = ~ treatment)

#(start)~~~~~~~~~~~~~~~~~~~~visuals for replicate cohesion~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# filter for low expression levels
## filters dds to remove sparse genes, which are defined here to be
## the set of genes that do not have 10 or more counts in at least n
## samples, where n = the typical minimum group size.
filtered_dds_treatment <- dds_treatment[rowSums(counts(dds_treatment) >= 10) >= 3, ]

## normalize data for heteroskedasticity
rld <- rlog(filtered_dds_treatment, blind = TRUE) # doesn't account for treatment groups

## create sample distance matrix
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

## create hierarchical heatmap of sample distances
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
ht <- pheatmap(sampleDistMatrix,
         main = "Sample Heatmap: All Samples",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

pdf("de_analysis/sample_dist_heatmap.pdf", width = 7, height = 5)
print(ht)
dev.off()

svg("de_analysis/sample_dist_heatmap.svg", width = 7, height = 5)
print(ht)
dev.off()

## create pca plot for sample groupings
pcaData <- plotPCA(rld, intgroup = "treatment", returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = treatment)) +
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

pdf("de_analysis/sample_dist_PCAplot.pdf", width = 7, height = 5)
print(pca_plot)
dev.off()

svg("de_analysis/sample_dist_PCAplot.svg", width = 7, height = 5)
print(pca_plot)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~visuals for replicate cohesion~~~~~~~~~~~~~~~~~~~~~~~~(end)#

# create test object
de_dds_treatment <- DESeq(dds_treatment)

# get experimental groups
treatments <- levels(coldata$treatment)

# generate all possible pairwise combinations
contrast_pairs <- as.data.frame(combn(treatments, 2))[c(2, 1),] # must be this order

# initialize storage
res_df_list <- list()
padj_list <- list()
l2fc_list <- list()
dir_list  <- list()

# generate test results
for (i in seq_len(ncol(contrast_pairs))){

  # df is wide form
  group1 <- contrast_pairs[1, i]
  group2 <- contrast_pairs[2, i]

  # ensures shorter exposure is base of current contrast
  if (grepl("[0-9]+", group1) && grepl("[0-9]+", group2)){
    g1 <- max(as.numeric(unlist(str_extract_all(group1, "[0-9]+"))))
    g2 <- max(as.numeric(unlist(str_extract_all(group2, "[0-9]+"))))

    if (g2 > g1) {
      temp_group <- group1
      group1 <- group2
      group2 <- temp_group
    } else if (g1 == g2) {
      if (nchar(group2) > nchar(group1)){
        temp_group <- group1
        group1 <- group2
        group2 <- temp_group
      }
    }
  }

  # ensures that 24_control is only first if control is second
  if (group1 == 'Control_24hr' && group2 != 'Control'){
    temp_group <- group1
    group1 <- group2
    group2 <- temp_group
  }



  # create name for current contrast
  contrast_name <- paste0(group1, "_vs_", group2)

  # extract results for contrast
  res <- results(de_dds_treatment, contrast = c("treatment", group1, group2))
  res_df <- as.data.frame(res)

  # add significance labels
  res_df$significance <- with(
    res_df, 
    ifelse(padj < 0.05 & log2FoldChange > 0, "upregulated",
      ifelse(padj < 0.05 & log2FoldChange < 0, "downregulated",
        "not significant"
      )
    )
  )

  # NAs are implicitly not significant test results
  res_df$significance[is.na(res_df$significance)] <- "not significant"

  # append to "mothership" list (full results per contrast)
  ## add the contrast name to column headers to distinguish
  res_labeled <- res_df
  colnames(res_labeled) <- paste0(contrast_name, "_", colnames(res_labeled))
  res_df_list[[contrast_name]] <- res_labeled

  # store specific attributes for the attribute DFs
  padj_list[[contrast_name]] <- res_df$padj
  l2fc_list[[contrast_name]] <- res_df$log2FoldChange
  dir_list[[contrast_name]]  <- res_df$significance

  # make volcano plot for contrast
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    labs(title = "Volcano Plot: Differentially Expressed Genes",
      subtitle = paste0("Experimental Contrast: ", contrast_name),
      x = expression(log[2]*(fold~change)),
      y = expression(-log[10]*(Adj.~P-value))
    ) +
    geom_point(aes(color = significance), alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("upregulated" = "red", 
                                  "downregulated" = "#145291", 
                                  "not significant" = "grey")) +
    theme_minimal()

  # write plot to pdf
  pdf_plot_name <- paste0("de_analysis/expandedContrasts_plots/",
                    contrast_name,
                    "_volcano_plot.pdf")
  pdf(pdf_plot_name, width = 7, height = 5)
  print(p)
  dev.off()

  # write plot to csv
  svg_plot_name <- paste0("de_analysis/expandedContrasts_plots/",
                    contrast_name,
                    "_volcano_plot.svg")
  svg(svg_plot_name, width = 7, height = 5)
  print(p)
  dev.off()
}



################################
#---- Create Attribute DFs ----#
################################

# create ens_ids vector
ens_ids <- rownames(res_df_list[[1]])
ens_ids_clean <- sub("\\..*", "", ens_ids)

# get gene symbols from ensemble ids
gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# create adjusted-p df
df_padj <- as.data.frame(do.call(cbind, padj_list), row.names = ens_ids)
df_padj$Geneid <- rownames(df_padj)
rownames(df_padj) <- NULL
df_padj$gene_symbol <- gene_symbols

# reorder padj_df columns
n <- ncol(df_padj)
df_padj <- df_padj[, c((n - 1):n, 1:(n - 2))]

# add tpm counts
df_padj <- merge(df_padj, tpm_df, by = "Geneid", all.x = TRUE)

# create log2(fold change) df
df_l2fc <- as.data.frame(do.call(cbind, l2fc_list), row.names = ens_ids)
df_l2fc$Geneid <- rownames(df_l2fc)
rownames(df_l2fc) <- NULL
df_l2fc$gene_symbol <- gene_symbols

# reorder l2fc_df columns
n <- ncol(df_l2fc)
df_l2fc <- df_l2fc[, c((n - 1):n, 1:(n - 2))]

# add tpm counts
df_l2fc <- merge(df_l2fc, tpm_df, by = "Geneid", all.x = TRUE)

# create direction df
df_direction <- as.data.frame(do.call(cbind, dir_list), row.names = ens_ids)
df_direction$Geneid <- rownames(df_direction)
rownames(df_direction) <- NULL
df_direction$gene_symbol <- gene_symbols

# reorder direction_df columns
n <- ncol(df_direction)
df_direction <- df_direction[, c((n - 1):n, 1:(n - 2))]

# add tpm counts
df_direction <- merge(df_direction, tpm_df, by = "Geneid", all.x = TRUE)

# created MotherShip df: all wald test data results
df_mothership <- do.call(cbind, res_df_list)
df_mothership$Geneid <- ens_ids
df_mothership$gene_symbol <- gene_symbols

# reorder mothership_df columns
n <- ncol(df_mothership)
df_mothership <- df_mothership[, c((n - 1):n, 1:(n - 2))]

# add tpm counts
df_mothership <- merge(df_mothership, tpm_df, by = "Geneid", all.x = TRUE)


# export dfs to .csv
write.csv(df_padj, snakemake@output[[1]])
write.csv(df_l2fc, snakemake@output[[2]])
write.csv(df_direction, snakemake@output[[3]])
write.csv(df_mothership, snakemake@output[[4]])
