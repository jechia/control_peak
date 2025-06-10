#!/bin/bash

# Usage:
#   bash genAnno.sh --gtf <gtf_file> --rna-list <rna_list_file> [--out-dir <output_dir>]
#
# Description:
#   Generates a BED annotation from a GENCODE GTF file, distinguishing UTRs, CDS, and exons,
#   and filtering RNA exons based on biotype.


# ------------------ Usage ------------------
usage() {
    echo "Usage: $0 --gtf <gtf_file> --rna-list <rna_list_file> [--out-dir <output_dir>]"
    echo ""
    echo "Arguments:"
    echo "  --gtf        Path to GENCODE GTF file (required)"
    echo "  --rna-list   Path to RNA biotype list (required)"
    echo "  --out-dir    Output directory (default: annotation)"
    exit 1
}

# ------------------ Parse Arguments ------------------
GTF=""
RNA_LIST=""
OUT_DIR="annotation"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gtf) GTF="$2"; shift ;;
        --rna-list) RNA_LIST="$2"; shift ;;
        --out-dir) OUT_DIR="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# ------------------ Validate ------------------
if [[ -z "$GTF" || -z "$RNA_LIST" ]]; then
    echo "❌ Error: --gtf and --rna-list are required."
    usage
fi

if [[ ! -f "$GTF" ]]; then
    echo "❌ GTF file not found: $GTF"
    exit 1
fi

if [[ ! -f "$RNA_LIST" ]]; then
    echo "❌ RNA list file not found: $RNA_LIST"
    exit 1
fi

# Output directory check
if [[ -d "$OUT_DIR" && "$(ls -A "$OUT_DIR")" ]]; then
    echo "⚠️  Warning: Output directory '$OUT_DIR' already exists and is not empty."
    read -p "❓ Continue and overwrite existing files? [y/N]: " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo "❌ Exiting."
        exit 1
    fi
else
    mkdir -p "$OUT_DIR"
fi

# ------------------ Begin Processing ------------------
echo "✅ Starting annotation generation..."
echo "GTF file:       $GTF"
echo "RNA list:       $RNA_LIST"
echo "Output dir:     $OUT_DIR"


# ----------- Step 0: Extract transcript and gene BED from GTF -----------
echo "Step 0: Extracting transcripts.bed genes.bed ..."
awk -v OFS="\t" '$3=="transcript" { match($0, /transcript_name "([^"]+)"/, a); print $1, $4-1, $5, a[1], $6, $7 }' "$GTF" > "$OUT_DIR/transcripts.bed"
awk -v OFS="\t" '$3=="transcript" { match($0, /gene_name "([^"]+)"/, a); print $1, $4-1, $5, a[1], $6, $7 }' "$GTF" > "$OUT_DIR/genes.bed"

# ----------- Step 1: Get stop codons -----------
echo "Step 1: Extracting stop codons..."
awk -v OFS="\t" '$3=="stop_codon" {
    if ($7 == "-") { s = $4 - 1; e = $5 }
    else           { s = $4; e = $5 + 1 }
    print $1, s, e, $20, 0, $7
}' "$GTF" | tr -d '";' > "$OUT_DIR/gencode_v46_stop_codon.bed"

# ----------- Step 2: Get UTR regions -----------
echo "Step 2: Extracting UTRs..."
awk -v OFS="\t" '$3=="UTR"{print $1, $4, $5, $20, 0, $7, $16}' "$GTF" | tr -d '";' > "$OUT_DIR/gencode_v46_UTR.bed"

# ----------- Step 3: Annotate UTRs as 5' or 3' -----------
echo "Step 3: Classifying UTRs as 5'/3'..."
awk -v OFS="\t" '
BEGIN {
    while ((getline < "'"$OUT_DIR/gencode_v46_stop_codon.bed"'") > 0) {
        start[$4] = $2; end[$4] = $3
    }
}
{
    if ($2 >= start[$4]) {
        type = ($6 == "-") ? "UTR5" : "UTR3"
    } else {
        type = ($6 == "-") ? "UTR3" : "UTR5"
    }
    print $1, $2, $3, $4"_"type"_"$4, $5, $6, $7
}' "$OUT_DIR/gencode_v46_UTR.bed" > "$OUT_DIR/gencode_v46_UTR_anno.bed"

# ----------- Step 4: Extract CDS -----------
echo "Step 4: Extracting CDS..."
awk -v OFS="\t" '$3=="CDS"{print $1, $4, $5, $20"_CDS_"$4, 0, $7, $16}' "$GTF" | tr -d '";' > "$OUT_DIR/gencode_v46_CDS.bed"

# ----------- Step 5: Combine CDS and UTR -----------
echo "Step 5: Combining CDS and UTR annotations..."
cat "$OUT_DIR/gencode_v46_CDS.bed" "$OUT_DIR/gencode_v46_UTR_anno.bed" > "$OUT_DIR/gencode_v46_pc_trans_anno.bed"

# ----------- Step 6: Filter exons by RNA biotype -----------
echo "Step 6: Filtering exons for selected RNA biotypes..."

# Extract exons where transcript_type is in RNA.list
gawk -v OFS="\t" '$3=="exon"{print $1, $4, $5, $20"_exon_"$4, 0, $7, $18, $16}' "$GTF" | tr -d '";' | \
awk -v OFS="\t" '
BEGIN {
    while ((getline < "'"RNA.list"'") > 0) {
        biotype[$1] = 1
    }
}
$7 in biotype {
    print $1, $2, $3, $4, $5, $6, $8
}' > "$OUT_DIR/gencode_v46_exons_RNA.bed"

# ----------- Step 7: Extract protein coding transcript names -----------
echo "Step 7: Generating transcript name list..."
cat "$OUT_DIR/gencode_v46_exons_RNA.bed" "$OUT_DIR/gencode_v46_pc_trans_anno.bed" | awk '{split($4, a, "_"); print a[1]}'  | sort | uniq > "$OUT_DIR/trans_name.list"

# ----------- Step 8: Filter transcript BED -----------
awk 'BEGIN{while((getline<"'"$OUT_DIR/trans_name.list"'")>0){d[$1]=1}}$4 in d' "$OUT_DIR/transcripts.bed" > "$OUT_DIR/gencode_v46_transcripts.bed"

# ----------- Step 9: Get full exon annotations for selected transcripts -----------
echo "Step 8: Extracting exon annotations for final transcript list..."
gawk -v OFS="\t" '$3=="exon"{print $1, $4, $5, $20, 0, $7, $18, $16}' "$GTF" | tr -d '";' | \
awk -v OFS="\t" '
BEGIN {
    while ((getline < "'"$OUT_DIR/trans_name.list"'") > 0) {
        trans[$1] = 1
    }
}
$4 in trans {
    print $1, $2, $3, $4"_exon_"$2, $5, $6, $8
}' > "$OUT_DIR/gencode_v46_exons.bed"

# ----------- Step 10: Final annotation file -----------
echo "Step 9: Creating final annotation file..."
cat "$OUT_DIR/gencode_v46_exons.bed" "$OUT_DIR/gencode_v46_pc_trans_anno.bed" > "$OUT_DIR/gencode_v46_anno.bed"

echo "✅ Done. Final annotation is at:"
echo "$OUT_DIR/gencode_v46_anno.bed"
