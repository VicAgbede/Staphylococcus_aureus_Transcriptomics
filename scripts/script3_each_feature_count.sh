#!/bin/bash

# ==============================================================================
# Pipeline Step 3: High-Throughput Feature Quantification (Unified Data Folder)
# ==============================================================================

# Strict error handling: exit immediately if any command fails
set -e

# --- Configuration & Paths ---
# Everything is now consolidated into your single Data folder
DATA_DIR="/mnt/e/Staph_Research/Data"
THREADS=4

# Official RefSeq GFF matched annotation located inside the Data folder
ANNOTATION="${DATA_DIR}/GCF_000013425.1_ASM1342v1_genomic.gff"

echo "========================================================================="
echo " Starting featureCounts quantification loop..."
echo " Working Directory: ${DATA_DIR}"
echo " Using annotation: $(basename "${ANNOTATION}")"
echo "========================================================================="

# --- Execution Loop ---
# Iterate through all sorted BAM alignment files in the data directory
for BAM_FILE in "${DATA_DIR}"/*_sorted.bam; do
    
    # Safeguard: Check if any BAM files actually exist matching the pattern
    [ -e "${BAM_FILE}" ] || { echo "No sorted BAM files found in ${DATA_DIR}. Exiting."; exit 1; }

    # Extract base filename without path or extensions for clean output tracking
    BASENAME=$(basename "${BAM_FILE}")
    FILENAME="${BASENAME%_sorted.bam}"

    echo ""
    echo "Processing Sample: ${FILENAME}"
    echo "-------------------------------------------------------------------------"

    # Run featureCounts with our verified robust bacterial parameters:
    # -F GFF       : Explicitly parses GFF3 format structural lines
    # -t gene      : Counts across whole gene intervals (captures UTR/leader regions)
    # -g Name      : Extracts standard gene symbols/IDs from the GFF attribute column
    # -s 0         : Enforces unstranded library mode
    # -M --fraction: Rescues multi-mapping reads (e.g., rRNAs) using fractional splits
    featureCounts -T "${THREADS}" \
                  -F GFF \
                  -a "${ANNOTATION}" \
                  -t gene \
                  -g Name \
                  -s 0 \
                  -M \
                  --fraction \
                  -o "${DATA_DIR}/${FILENAME}_counts.txt" \
                  "${BAM_FILE}"

    echo "Finished quantification for ${FILENAME}."
done

echo ""
echo "========================================================================="
echo " Pipeline Complete! All count tables saved directly to: ${DATA_DIR}"
echo "========================================================================="