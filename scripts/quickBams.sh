#!/usr/bin/env bash
set -euo pipefail

REF=~/genomeRefs/hg38/hg38.fa
THREADS=8

# Folder containing FASTQ files is passed as the first argument
FASTQ_DIR="$1"
ALIGNMENTS_DIR="$2"

# Loop over all .fastq or .fastq.gz files in the folder
for READS in "$FASTQ_DIR"/*.fastqsanger.gz; do
    # Skip if no files match
    [ -e "$READS" ] || continue

    BASENAME=$(basename "$READS") # strips directory
    SAMPLE=${BASENAME%%.*}   # strips .fastq, .fastq.gz, etc.
    OUT="${ALIGNMENTS_DIR}/${SAMPLE}.sorted.bam"

    echo "Processing $SAMPLE …"

    bwa mem -t "$THREADS" "$REF" "$READS" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$OUT"

    samtools index "$OUT"
done
