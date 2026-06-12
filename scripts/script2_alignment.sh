#!/bin/bash

# Set input and output directory to the current location
INPUT_DIR=$(dirname "${BASH_SOURCE[0]}")
OUTPUT_DIR="${INPUT_DIR}"

# OPTIMIZED HARDWARE TARGET: Balanced perfectly for your 4 WSL cores
THREADS=4
SAMTOOLS_THREADS=2

# Loop through all trimmed fastq files in the directory
for FASTQ_FILE in "${INPUT_DIR}"/*_trimmed.fastq
do
    # Safety check: skip the loop if no matching files are found
    [ -e "$FASTQ_FILE" ] || continue

    # Get the filename without the extension suffix
    FILENAME=$(basename "${FASTQ_FILE}" _trimmed.fastq)

    echo "===================================================="
    echo "Processing Sample Alignment: ${FILENAME}"
    echo "===================================================="

    # Align the reads and stream directly into a multi-threaded samtools sort
    hisat2 -q -x "${INPUT_DIR}/nctc8325_index" -U "${FASTQ_FILE}" -p "${THREADS}" | \
        samtools sort -@ "${SAMTOOLS_THREADS}" -o "${OUTPUT_DIR}/${FILENAME}_sorted.bam" -

    echo "Indexing BAM for downstream analysis..."
    samtools index "${OUTPUT_DIR}/${FILENAME}_sorted.bam"

done

echo "===================================================="
echo "Pipeline completed successfully with indices generated!"
echo "===================================================="