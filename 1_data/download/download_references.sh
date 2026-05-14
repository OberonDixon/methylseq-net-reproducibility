#!/usr/bin/env bash
set -euo pipefail

GENOMES_DIR=$(methylseqnet-repro-path genomes)
HAPLOTYPES_DIR=$(methylseqnet-repro-path haplotypes)
ANNOTATIONS_DIR=$(methylseqnet-repro-path annotations)

mkdir -p "$GENOMES_DIR"
mkdir -p "$HAPLOTYPES_DIR"
mkdir -p "$ANNOTATIONS_DIR"

# Download reference genome and haplotype information

wget -O "$GENOMES_DIR/hg38.2bit" \
    https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
wget -O "$GENOMES_DIR/hg38.fa.gz" \
    https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip "$GENOMES_DIR/hg38.fa.gz"
samtools faidx "$GENOMES_DIR/hg38.fa"

wget -O "$HAPLOTYPES_DIR/NA12878.vcf.gz" \
    https://hgdownload.soe.ucsc.edu/gbdb/hg38/platinumGenomes/NA12878.vcf.gz

wget -O "$HAPLOTYPES_DIR/NA12878.vcf.gz.tbi" \
    https://hgdownload.soe.ucsc.edu/gbdb/hg38/platinumGenomes/NA12878.vcf.gz.tbi

# Download external annotations

GTF="$ANNOTATIONS_DIR/gencode.v49.basic.annotation.gtf"

wget -O "${GTF}.gz" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz

zcat "${GTF}.gz" \
    | grep -v "^#" \
    | sort -k1,1 -k4,4n \
    | bgzip > "${GTF}.sorted.gz"

tabix -p gff "${GTF}.sorted.gz"

rm "${GTF}.gz"

wget -O "$ANNOTATIONS_DIR/hg38.archetype_motifs.v1.0.bed.gz" \
    https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz

wget -O "$ANNOTATIONS_DIR/hg38.archetype_motifs.v1.0.bed.gz.tbi" \
    https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz.tbi

wget -O "$ANNOTATIONS_DIR/41586_2022_5580_MOESM4_ESM.xlsx" \
    https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05580-6/MediaObjects/41586_2022_5580_MOESM4_ESM.xlsx

wget -O "$ANNOTATIONS_DIR/tabula_sapiens_2.0_table4_TFs-by-celltype.xlsx" \
    "https://www.biorxiv.org/content/biorxiv/early/2025/08/27/2024.12.03.626516/DC21/embed/media-21.xlsx?download=true"

# Copy pre-packaged annotations

cp ../../reference/transcription_factor_metadata.csv "$ANNOTATIONS_DIR/transcription_factor_metadata.csv"

cp ../../reference/imprinted_genes.csv "$ANNOTATIONS_DIR/imprinted_genes.csv"

cp -r ../../reference/pwms "$ANNOTATIONS_DIR/pwms"

cp -r ../../borzoi-human-splits "$ANNOTATIONS_DIR/borzoi-human-splits"