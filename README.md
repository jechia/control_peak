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

### Generate Control Peaks

```bash
# Annotate the peaks
bedtools intersect -a <sample>_peaks.bed \
                   -b annotation/gencode_v46_transcripts.bed \
                   -f 1 -wa -wb -s > anno/<sample>.bed

# Generate control peaks
python eCLIP_control_v4.py -i anno/<sample>.bed \
                           -a annotation/gencode_v46_anno.bed \
                           -g annotation/genes.bed \
                           -p 10

# Get fasta for the peak and control peaks
sh getFasta.sh <sample>
```
