#!/bin/bash

# Genomic Data Processing Script
# Processes BED files to extract peak and control sequences
# Usage: ./getFasta.sh <sample_name> [reference_genome] [input_dir] [output_dir]

# Check if required argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_name> [reference_genome] [input_dir] [output_dir]"
    echo "  sample_name: Required - name of the sample to process"
    echo "  reference_genome: Optional - path to reference genome (default: ./annotation/GRCh38.p14.genome.fa)"
    echo "  input_dir: Optional - input directory (default: new_result)"
    echo "  output_dir: Optional - output directory (default: new_fasta)"
    exit 1
fi

# Command line arguments
SAMPLE_NAME="$1"
REFERENCE_GENOME="${2:-./annotation/GRCh38.p14.genome.fa}"
INPUT_DIR="${3:-new_result}"
OUTPUT_DIR="${4:-new_fasta}"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "Processing sample: $SAMPLE_NAME"
echo "Reference genome: $REFERENCE_GENOME"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "---"

# Process peak regions
# Extract columns 1-4, add score column (0), and strand column
echo "Creating peak BED file..."
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, 0, $5}' \
    "${INPUT_DIR}/${SAMPLE_NAME}_control.bed" > "${INPUT_DIR}/${SAMPLE_NAME}_peak.bed"

# Extract FASTA sequences for peak regions
echo "Extracting peak sequences..."
bedtools getfasta \
    -s \
    -fi "$REFERENCE_GENOME" \
    -bed "${INPUT_DIR}/${SAMPLE_NAME}_peak.bed" \
    -name > "${OUTPUT_DIR}/${SAMPLE_NAME}_peak.fasta"

# Process control/pair regions  
# Extract columns 1, 6-7, 4, add score column (0), and strand column
echo "Creating control BED file..."
awk 'BEGIN{OFS="\t"} {print $1, $6, $7, $4, 0, $5}' \
    "${INPUT_DIR}/${SAMPLE_NAME}_control.bed" > "${INPUT_DIR}/${SAMPLE_NAME}_pair.bed"

# Extract FASTA sequences for control regions
echo "Extracting control sequences..."
bedtools getfasta \
    -s \
    -fi "$REFERENCE_GENOME" \
    -bed "${INPUT_DIR}/${SAMPLE_NAME}_pair.bed" \
    -name > "${OUTPUT_DIR}/${SAMPLE_NAME}_control.fasta"

echo "Processing complete for sample: $SAMPLE_NAME"