#!/bin/bash

# Set paths to directories and files
bam_dir="$HOME/HIV-Project/BAMS"  # Path to BAM files
reference_genome="$HOME/HIV-Project/references/reference_genome.fa"  # Path to HIV-1 reference genome
gff_file="$HOME/HIV-Project/references/reference_genome.gff"  # Path to GFF annotation file
vcf_dir="$HOME/HIV-Project/ivar_VCFiles"  # Output directory for individual VCF files
merged_vcf="$HOME/HIV-Project/merged_genotyping.vcf.gz"  # Merged output VCF file
filtered_vcf="$HOME/HIV-Project/filtered_genotyping.vcf.gz"  # Filtered output VCF file

# Number of CPU threads to use
threads=60  # Adjust based on your system's core count

# Create the VCF output directory if it doesn't exist
mkdir -p "$vcf_dir"

# Step 1: Sort BAM files and call variants with ivar
for bam in "$bam_dir"/*.bam; do
    base_name=$(basename "$bam" .bam)
    sorted_bam="$vcf_dir/${base_name}_sorted.bam"
    vcf_file="$vcf_dir/${base_name}_variants.vcf"

    echo "==================================================="
    echo "Processing sample: $base_name"
    echo "Input BAM file: $bam"
    echo "Output VCF file: $vcf_file"
    echo "Using $threads CPU threads"
    echo "==================================================="

    # Step 1: Sort the BAM file (if not already sorted)
    samtools sort -o "$sorted_bam" "$bam"

    # Step 2: Generate VCF file using samtools and ivar with GFF annotation and reference genome
    samtools mpileup -A -d 0 -Q 0 -B -f "$reference_genome" "$sorted_bam" | \
    ivar variants -p "$vcf_file" -q 20 -t 0.01 -m 10 -r "$reference_genome" -g "$gff_file"

    echo "---------------------------------------------------"
    echo "Finished processing sample: $base_name"
    echo "---------------------------------------------------"
    echo ""
done

# Step 2: Merging Individual VCF Files
echo "==================================================="
echo "Merging individual VCF files into one VCF"
echo "==================================================="

bcftools merge -m all -O z -o "$merged_vcf" "$vcf_dir"/*.vcf

echo "---------------------------------------------------"
echo "Finished merging VCF files"
echo "---------------------------------------------------"

# Step 3: Filtering the Variants Based on Quality and Depth
echo "==================================================="
echo "Filtering variants based on quality and depth"
echo "==================================================="

bcftools filter -s LowQual -e '%QUAL<20 || DP<10' -Oz -o "$filtered_vcf" "$merged_vcf"

echo "---------------------------------------------------"
echo "Finished filtering variants"
echo "---------------------------------------------------"

# Step 4: Indexing the Final VCF File
echo "==================================================="
echo "Indexing the final filtered VCF file"
echo "==================================================="

bcftools index "$filtered_vcf"

echo "---------------------------------------------------"
echo "Variant calling pipeline completed. Final VCF is indexed and ready for analysis."
echo "==================================================="
