mkdir -p annotation
cd annotation

# Download the annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz

# Uncompress
gunzip gencode.v46.primary_assembly.annotation.gtf.gz

# Download the genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz

# Uncompress
gunzip GRCh38.p14.genome.fa.gz

cd ..