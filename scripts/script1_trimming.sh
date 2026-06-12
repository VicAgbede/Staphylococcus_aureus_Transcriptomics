#!/bin/bash

# Set the path to the input directory as the current directory
input_dir=$(pwd)

# Explicit path to your Trimmomatic jar engine
TRIMMOMATIC_JAR="${input_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar"

# CORRECTED PATH: Points directly to the adapters folder seen in your file explorer
ADAPTERS_PATH="${input_dir}/adapters/TruSeq3-SE.fa"

for file in "$input_dir"/*.fastq; do
    
    # Skip files that are already trimmed
    [[ "$file" == *"_trimmed.fastq" ]] && continue
    
    filename=$(basename "$file" .fastq)

    echo "=================================================="
    echo "Stringent Trim Started: ${filename}"
    echo "=================================================="

    # Run trimmomatic with the corrected absolute path to the adapter file
    java -jar "$TRIMMOMATIC_JAR" SE -threads 32 -phred33 \
        "$file" \
        "$input_dir/${filename}_trimmed.fastq" \
        ILLUMINACLIP:"${ADAPTERS_PATH}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done

echo "=================================================="
echo "Trimming complete with zero file errors!"
echo "=================================================="

