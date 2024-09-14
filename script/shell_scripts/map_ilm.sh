#!/bin/bash

# Set the directory containing the processed FASTQ files
fastq_dir="$HOME/HIV-Project/ilm_processed_output"
## Reference genome
reference_genome="$HOME/HIV-Project/references/reference_genome.fa"
## Result directory
results_dir="$HOME/HIV-Project/processed_ilm_bam_files"
### Capture mapping statistics
stats_dir="$results_dir/statistics"

# Number of threads to use
threads=60

# Create the results and statistics directories if they don't exist
mkdir -p "$results_dir"
mkdir -p "$stats_dir"

# Extract the directory containing the reference genome
reference_dir=$(dirname "$reference_genome")

# Check if BWA-MEM2 index files exist, if not, generate them
if [ ! -f "${reference_genome}.bwt.2bit.64" ]; then
    echo "Index files for the reference genome not found. Generating index..."
    bwa-mem2 index -p "$reference_dir/$(basename "$reference_genome")" "$reference_genome"
    echo "Index generation completed."
else
    echo "Index files for the reference genome found. Skipping index generation."
fi

# Loop through each paired FASTQ file
for fq1 in "$fastq_dir"/*_1.fastq; do
    # Extract the sample name by removing the '.fastp_1.fastq' suffix
    base_name=$(basename "$fq1" .fastp_1.fastq)
    fq2="${fastq_dir}/${base_name}.fastp_2.fastq"
    
    # Check if the reverse read file exists
    if [ ! -f "$fq2" ]; then
        echo "Error: $fq2 does not exist."
        exit 1
    fi

    # Define file paths using the sample name
    sam_file="${results_dir}/${base_name}.sam"
    bam_file="${results_dir}/${base_name}.bam"
    sorted_bam_file="${results_dir}/${base_name}.sorted.bam"
    fixmate_bam_file="${results_dir}/${base_name}.fixmate.bam"
    marked_bam_file="${results_dir}/${base_name}.marked.bam"
    stats_file="${stats_dir}/${base_name}_mapping_stats.txt"
    depth_file="${stats_dir}/${base_name}_depth.txt"

    # Define the read group information
    read_group="@RG\tID:${base_name}\tSM:${base_name}\tPL:ILLUMINA\tLB:${base_name}_lib\tPU:${base_name}_unit"

    # Align with BWA-MEM2 using settings for short paired-end reads and 60 threads
    bwa-mem2 mem -t $threads -R "$read_group" "$reference_genome" "$fq1" "$fq2" > "$sam_file"
    
    # Convert SAM to BAM
    samtools view -@ $threads -Sb "$sam_file" > "$bam_file"
    
    # Sort the BAM file by queryname before running fixmate
    samtools sort -@ $threads -n -o "$sorted_bam_file" "$bam_file"
    
    # Add or fix mate-pair information in the BAM file
    samtools fixmate -m -@ $threads "$sorted_bam_file" "$fixmate_bam_file"
    
    # Sort the BAM file by coordinates after fixmate (needed by markdup)
    samtools sort -@ $threads -o "$sorted_bam_file" "$fixmate_bam_file"
    
    # Index the sorted BAM file before marking duplicates
    samtools index "$sorted_bam_file"

    # Check if the sorted BAM file exists and is non-empty before running markdup
    if [ -s "$sorted_bam_file" ]; then
        # Mark duplicates
        samtools markdup -@ $threads "$sorted_bam_file" "$marked_bam_file"
        
        # Index the final BAM file
        samtools index "$marked_bam_file"

        # Capture mapping statistics
        samtools flagstat "$marked_bam_file" > "$stats_file"
        
        # Capture depth/coverage information
        samtools depth -a "$marked_bam_file" > "$depth_file"
        
        # Optional: Capture summary coverage statistics (requires samtools >= 1.10)
        samtools coverage "$marked_bam_file" >> "$stats_file"
        
        # Clean up intermediate files
        rm "$sam_file" "$bam_file" "$sorted_bam_file" "$fixmate_bam_file"
    else
        echo "Error: Sorted BAM file $sorted_bam_file does not exist or is empty. Skipping markdup for $base_name."
    fi
done

echo "Short-read mapping pipeline completed. Processed BAM files and statistics are stored in $results_dir."
