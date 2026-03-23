suppressPackageStartupMessages({

  library(DESeq2)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
})

dir.create("de_analysis/expandedContrasts_plots-tables")
##############################################################
#---- read in counts df and formats counts matrix object ----#
##############################################################

counts_df <- read.csv(snakemake@input[[1]])

# create raw counts matrix
counts_matrix <- counts_df[, !grepl(".nTPM", colnames(counts_df))]
tpm_df <- counts_df[, c("Geneid", grep(".nTPM", colnames(counts_df), value = TRUE))]

rownames(counts_matrix) <- counts_matrix$Geneid
counts_matrix$Geneid <- NULL

# subsets to only count-value columns
counts_matrix <- counts_matrix[, -c(1:5)]


##################################################
#---- generate sample metadata (full samples)----#
##################################################

samples <- colnames(counts_matrix)
# Extract CO treatment time if present
co_time <- str_extract(samples, "(?<=CO_)[0-9]+hr")
co_time[is.na(co_time)] <- "0hr"

# Extract APAP treatment time if present
apap_time <- str_extract(samples, "(?<=APAP_)[0-9]+hr")
apap_time[is.na(apap_time)] <- "0hr"


treatment <- sapply(samples, function(x) {sub("_.$", "", x)})
# Extract replicate (assume last number in name)
replicate <- str_extract(samples, "[0-9]+$")

# construct colData object containing metadata
time_factor_levels <- c("0hr", "4hr", "8hr", "24hr")
coldata <- data.frame(
  row.names = samples,
  treatment = factor(treatment,
                     levels = c("Control",
                                setdiff(unique(treatment), "Control"))),
  # as.numeric treats duration as a factor rather than category
  co_time = factor(co_time, levels = time_factor_levels), 
  apap_time = factor(apap_time, levels = time_factor_levels),
  replicate = factor(replicate)
)

#####################################################################
#---- Differential Analysis: Experimental Group (All Contrasts) ----#
#####################################################################
dds_treatment = DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~ treatment)

de_dds_treatment <- DESeq(dds_treatment)

# get experimental groups
treatments <- levels(coldata$treatment)
# Generate all possible pairwise combinations
contrast_pairs <- as.data.frame(combn(treatments, 2))[c(2, 1),] #must be this order

# initialize storage
res_df_list <- list()
padj_list <- list()
l2fc_list <- list()
dir_list  <- list()

# generate test results
for(i in seq_len(ncol(contrast_pairs))){
  group1 <- contrast_pairs[1, i]
  group2 <- contrast_pairs[2, i]

  # ensures longer exposures relative to lower exposures
  if (grepl("[0-9]+", group1) & grepl("[0-9]+", group2)){
    g1 <- max(as.numeric(unlist(str_extract_all(group1, "[0-9]+"))))
    g2 <- max(as.numeric(unlist(str_extract_all(group2, "[0-9]+"))))

    if (g2 > g1){
      temp_group <- group1
      group1 <- group2
      group2 <- temp_group
    } else if (g1 == g2){
      if (nchar(group2) > nchar(group1)){
        temp_group <- group1
        group1 <- group2
        group2 <- temp_group
      }
    }
  }

  # ensures that 24_control is only first if control is second
  if (group1 == 'Control_24hr' & group2 != 'Control'){
    temp_group <- group1
    group1 <- group2
    group2 <- temp_group
  }



  # Create a clean name for this contrast
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

  res_df$significance[is.na(res_df$significance)] <- "not significant"

  # append to mothership list (full results per contrast)
  # add the contrast name to column headers to distinguish
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

  pdf_plot_name <- paste0("de_analysis/expandedContrasts_plots-tables/",
                    contrast_name,
                    "_volcano_plot.pdf")
  pdf(pdf_plot_name, width = 7, height = 5)
  print(p)
  dev.off()

  svg_plot_name <- paste0("de_analysis/expandedContrasts_plots-tables/",
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
