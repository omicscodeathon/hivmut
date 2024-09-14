#!/bin/bash

# Set the directory containing the processed FASTQ files
fastq_dir="$HOME/HIV-Project/ont_processed_output"
## Reference genome
reference_genome="$HOME/HIV-Project/references/reference_genome.fa"
## Result directory
results_dir="$HOME/HIV-Project/processed_ont_bam_files"
### Capture mapping statistics
stats_dir="$results_dir/statistics"

# Number of threads to use
threads=60

# Create the results and statistics directories if they don't exist
mkdir -p "$results_dir"
mkdir -p "$stats_dir"

# Loop through each FASTQ file
for fq in "$fastq_dir"/*_preprocessed.fastq.gz; do
    # Extract the sample name by removing the '.fastq.gz' suffix
    base_name=$(basename "$fq" _preprocessed.fastq.gz)
    
    # Define file paths using the sample name
    sam_file="${results_dir}/${base_name}.sam"
    bam_file="${results_dir}/${base_name}.bam"
    sorted_bam_file="${results_dir}/${base_name}.sorted.bam"
    marked_bam_file="${results_dir}/${base_name}.marked.bam"
    stats_file="${stats_dir}/${base_name}_mapping_stats.txt"
    depth_file="${stats_dir}/${base_name}_depth.txt"

    # Define the read group information
    read_group="@RG\tID:${base_name}\tSM:${base_name}\tPL:ONT\tLB:${base_name}_lib\tPU:${base_name}_unit"

    # Align with minimap2 using settings for ONT reads and 60 threads
    minimap2 -ax map-ont -t $threads -R "$read_group" "$reference_genome" "$fq" > "$sam_file"
    
    # Convert SAM to BAM
    samtools view -@ $threads -Sb "$sam_file" > "$bam_file"
    
    # Sort the BAM file by coordinates
    samtools sort -@ $threads -o "$sorted_bam_file" "$bam_file"
    
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
        rm "$sam_file" "$bam_file"
    else
        echo "Error: Sorted BAM file $sorted_bam_file does not exist or is empty. Skipping markdup for $base_name."
    fi
done

echo "ONT mapping pipeline completed. Processed BAM files and statistics are stored in $results_dir."
