#' script doc string:
#' Script optimizes counts data for differential analysis
#' preprocessing and execution. New column names encode experimental
#' design in consistent parsable format. TPM count normalization
#' reflects request for format from collaborators.

# import counts data generated with feature counts
counts_df <- read.delim(snakemake@input[[1]],
                        comment.char = "#",
                        check.names = FALSE)

###############################################
#---- renames columns to parsable strings ----#
###############################################

# save sequence of old columns names
old_colnames <- colnames(counts_df)

# create vector of new column names
new_colnames <- c(
"Geneid","Chr",
"Start","End",
"Strand","Length",
"CO_4hr_APAP_24hr_1","CO_4hr_APAP_24hr_2",
"CO_24hr_3","CO_24hr_2",
"CO_24hr_1","Control_24hr_3",
"Control_24hr_2","Control_24hr_1",
"CO_4hr_APAP_4hr_3","CO_4hr_APAP_4hr_2",
"CO_4hr_APAP_4hr_1","CO_8hr_APAP_8hr_3",
"CO_8hr_APAP_8hr_2","CO_8hr_APAP_8hr_1",
"CO_4hr_3","CO_4hr_2",
"CO_4hr_1","APAP_4hr_3",
"APAP_4hr_2","APAP_4hr_1",
"Control_3","Control_2",
"Control_1"
)

# write map of new and old colname pairings to csv (for verification)
write.csv(cbind(old_colnames, new_colnames), "updated_sampleNames_map.csv")

# update new colnames
colnames(counts_df) <- new_colnames

# remove metadata from counts matrix: keeps only sample raw counts and geneids
counts_matrix <- counts_df[, c(1, 7:ncol(counts_df))]

# write counts only matrix with updated column identifiers to csv
write.csv(counts_matrix, snakemake@output[[1]], row.names = FALSE)

#############################################################################
#---- create csv with tpm: keeps metadata for future functional mapping ----#
#############################################################################

# generate length vector from annotations provided by featureCounts
gene_lengths <- counts_df$"Length"

# print number of transcripts without length data
print(paste0("Number of transcripts without length data:",
             sum(is.na(gene_lengths))))

# gene lengths as kilobases
kb_lengths <- gene_lengths / 1000

# isolate sample columns to mutate into added tpm sample counts
conv_cols <- colnames(counts_matrix)[-1]

# create new column with tpm counts for each sample
for (col in conv_cols){
  # create unit-id'd new column name
  newName <- paste0(col,'.nTPM', collapse = "")

  # vectorized transformation of raw counts to tpm normalization
  rpk <- counts_df[[col]] / kb_lengths
  tpm <- (rpk/sum(rpk)) * 1e6
  counts_df[[newName]] <- tpm
}
# write expanded dataframe to csv
write.csv(counts_df,snakemake@output[[2]], row.names = FALSE)
