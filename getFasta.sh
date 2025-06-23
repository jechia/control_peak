#!/bin/bash

# Genomic Data Processing Script
# Processes BED files to extract peak and control sequences
# Usage: ./getFasta.sh -s <sample_name> -r [reference_genome] -i [input_dir] -o [output_dir]

REFERENCE_GENOME="./annotation/GRCh38.p14.genome.fa"
INPUT_DIR="result"
OUTPUT_DIR="fasta"

# Help message
usage() {
    echo "Usage: $0 -s <sample_name> [-r reference_genome] [-i input_dir] [-o output_dir]"
    echo "  -s SAMPLE_NAME       Required: name of the sample to process"
    echo "  -r REFERENCE_GENOME  Optional: path to reference genome (default: $REFERENCE_GENOME)"
    echo "  -i INPUT_DIR         Optional: input directory (default: $INPUT_DIR)"
    echo "  -o OUTPUT_DIR        Optional: output directory (default: $OUTPUT_DIR)"
    exit 1
}

# Parse options
while getopts ":s:r:i:o:" opt; do
  case ${opt} in
    s ) SAMPLE_NAME="$OPTARG"
      ;;
    r ) REFERENCE_GENOME="$OPTARG"
      ;;
    i ) INPUT_DIR="$OPTARG"
      ;;
    o ) OUTPUT_DIR="$OPTARG"
      ;;
    \? )
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    : )
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if required argument is provided
if [ -z "$SAMPLE_NAME" ]; then
    echo "Error: sample name is required."
    usage
fi

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory '$OUTPUT_DIR' does not exist. Creating it..."
    mkdir -p "$OUTPUT_DIR"
fi

echo "Processing sample: ${INPUT_DIR}/$SAMPLE_NAME"
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

echo "Processing complete for sample: ${SAMPLE_NAME}"
echo "FASTA files saved to: ${OUTPUT_DIR}"