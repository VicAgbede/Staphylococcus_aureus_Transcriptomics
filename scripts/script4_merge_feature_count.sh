#!/bin/bash

# ==============================================================================
# Pipeline Step 4: Safe Count Matrix Merging (featureCounts Delimited Output)
# ==============================================================================

# Strict error handling: exit immediately if any command fails
set -e

# Define your data folder where all *_counts.txt files are stored
DATA_DIR="/mnt/e/Staph_Research/Data"
OUTPUT="${DATA_DIR}/merged_matrix.csv"

# Navigate to your data workspace
cd "${DATA_DIR}"

echo "========================================================================="
echo " Starting matrix merge in: ${DATA_DIR}"
echo "========================================================================="

# 1. Grab the first true sample file to initialize the master matrix 
# (This avoids picking up any 'test_' remnants if they are floating around)
FIRST_FILE=$(ls DRR*_counts.txt | head -n 1)

if [ -z "${FIRST_FILE}" ]; then
    echo "Error: No matching sample count files (DRR*_counts.txt) found!"
    exit 1
fi

echo "Initializing master matrix with baseline file: ${FIRST_FILE}"

# Extract Geneid (Col 1) and the sample's count file name header (Col 7)
# grep -v "^#" completely strips away the software version commentary line
grep -v "^#" "${FIRST_FILE}" | awk 'NR==1 {print "Geneid," $7} NR>1 {print $1 "," $7}' > "${OUTPUT}"

# 2. Loop and merge the remaining sample files sequentially
for file in DRR*_counts.txt; do
    if [ "$file" != "${FIRST_FILE}" ]; then
        echo "Appending column from sample: ${file}"
        
        # Extract just Column 7 (Count values) along with its sample header name
        grep -v "^#" "$file" | awk 'NR==1 {print $7} NR>1 {print $7}' | paste -d, "${OUTPUT}" - > temp_merge.csv
        
        # Atomically update our master matrix file
        mv temp_merge.csv "${OUTPUT}"
    fi
done

echo ""
echo "========================================================================="
echo " Succession! Your pristine matrix is ready at: ${OUTPUT}"
echo "========================================================================="

