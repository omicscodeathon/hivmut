#!/bin/bash

# Define the output directory
output_dir="illumina"

# Create the directory if it does not exist
mkdir -p "$output_dir"

# Iterate over each accession in the file and download the corresponding FASTQ files
while read -r accession; do
    fastq-dump --split-files --outdir "$output_dir" "$accession"
done < samples.txt
