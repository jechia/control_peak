# Control Peak Pipeline for eCLIP-seq

**Author:** Yue Hu  
**Date:** June 10, 2025

## Introduction

A pipeline to generate control peaks based on eCLIP-seq peaks. Each control peak matches the transcript region and gene with the original eCLIP-seq peak.

## Installation

### Quick Install

```bash
git clone https://github.com/jechia/control_peak.git
cd control_peak
```

### Requirements

The following software and tools are required:

- Bash 4.0 or higher
- Standard Unix tools (awk, sed, grep)
- bedtools
- python-intervals

## Basic Usage

### Prepare Annotation Files

```bash
# Download annotation from GENCODE database (GENCODE_V46)
sh downAnno.sh

# Generate customized annotation based on GENCODE
sh genAnno.sh --gtf annotation/gencode.v46.primary_assembly.annotation.gtf \
              --rna-list RNA.list
```
#### Note: This Annotation generation part is currently only available for linux system. For other systems, please use the generated annotation for human based on GENCODE V46 in "annotation.zip" or "annotation.tar.gz".

### Generate Control Peaks
We're using this with our example dataset K562_RBFOX2_peaks.bed

```bash
# Annotate the peaks
mkdir anno
bedtools intersect -a K562_RBFOX2_peaks.bed \
                   -b annotation/gencode_v46_transcripts.bed \
                   -f 1 -wa -wb -s > anno/K562_RBFOX2.bed

# Generate control peaks
python eCLIP_control_v4.py -i anno/K562_RBFOX2.bed \
                           -a annotation/gencode_v46_anno.bed \
                           -g annotation/genes.bed \
                           -p 10

# Get fasta for the peak and control peaks
sh getFasta.sh -s K562_RBFOX2 -a ./annotation/GRCh38.p14.genome.fa -i result -o fasta
```
