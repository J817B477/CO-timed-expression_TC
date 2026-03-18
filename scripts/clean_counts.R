
#imports counts data generated with feature counts
counts_df = read.delim(snakemake@input[[1]], comment.char = "#", check.names = FALSE)


# rename columns: removes files extensions, labeling be sample name only
old_colnames = colnames(counts_df)


new_colnames= c(
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
"Control_1", 'sample_22',
'sample_23', 'sample_24',
'sample_25'
)
write.csv(cbind(old_colnames, new_colnames), "updated_sampleNames_map.csv")
colnames(counts_df) = new_colnames

# remove metatadata from counts matrix: keeps only sample raw counts and geneids
counts_matrix = counts_df[, c(1, 7:ncol(counts_df))]

write.csv(counts_matrix, snakemake@output[[1]], row.names = F)


# create csv with tpm: keeps all metadata for future functional mapping
gene_lengths = counts_df$"Length"

# indicates where nulls exist
print(paste0("Number of transcripts without length data:", sum(is.na(gene_lengths))))

# lengths as kilobases
kb_lengths = gene_lengths/1000

# isolates sample columns to mutate into added tpm sample counts
conv_cols = colnames(counts_matrix)[-1]

# creates new column with tpm counts for each sample
for (col in conv_cols){
    newName = paste0(col,'.nTPM', collapse = "")
    rpk = counts_df[[col]]/kb_lengths
    tpm = (rpk/sum(rpk)) * 1e6
    counts_df[[newName]] = tpm
}
# writes expanded dataframe to csv
write.csv(counts_df,snakemake@output[[2]],row.names = F)
