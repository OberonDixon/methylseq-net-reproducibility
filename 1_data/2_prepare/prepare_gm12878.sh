#!/usr/bin/env bash
set -euo pipefail
RNA_DIR=$(methylseqnet-repro-path rna_tracks)
HAPLOTYPES_DIR=$(methylseqnet-repro-path haplotypes)
GENOMES_DIR=$(methylseqnet-repro-path genomes)
VCF_FILE="$HAPLOTYPES_DIR/NA12878.vcf.gz"

bcftools consensus -f "$GENOMES_DIR/hg38.fa" -H 1 \
    -i 'TYPE="snp"' ${VCF_FILE} \
    > "${HAPLOTYPES_DIR}/NA12878.hg38.hp1"
bcftools consensus -f "$GENOMES_DIR/hg38.fa" -H 2 \
    -i 'TYPE="snp"' ${VCF_FILE} \
    > "${HAPLOTYPES_DIR}/NA12878.hg38.hp2"

echo "haplotagging isoseq data"
echo "----------------------------------------------"

echo "indexing isoseq bam"
samtools index "$RNA_DIR/GM12878.kinnex.bam"

echo "running whatshap"
whatshap haplotag \
    --reference "$GENOMES_DIR/hg38.fa" \
    --ignore-read-groups \
    --skip-missing-contigs \
    "$HAPLOTYPES_DIR/NA12878.vcf.gz" \
    "$RNA_DIR/GM12878.kinnex.bam" \
    | samtools view -@ 4 -h \
    | tee >(samtools view -@ 2 -b -d HP:1 -o "$RNA_DIR/GM12878.kinnex.HP1.bam") \
    | samtools view -@ 2 -b -d HP:2 -o "$RNA_DIR/GM12878.kinnex.HP2.bam"

echo "indexing phased bams"
samtools index "$RNA_DIR/GM12878.kinnex.HP1.bam"
samtools index "$RNA_DIR/GM12878.kinnex.HP2.bam"

echo "isoseq collapse -> bedgz tss counts"
echo "----------------------------------------------"

SCRIPT_DIR="$(dirname "$0")"

for label in "" "HP1" "HP2"; do
    if [[ -z "$label" ]]; then
        bam="${RNA_DIR}/GM12878.kinnex.bam"
        gff="${RNA_DIR}/GM12878.kinnex.collapsed.no5exon.gff"
        abundance="${RNA_DIR}/GM12878.kinnex.collapsed.no5exon.abundance.txt"
        bed="${RNA_DIR}/GM12878.kinnex.no5exon.tss.counts.bed.gz"
        tag="unphased"
    else
        bam="${RNA_DIR}/GM12878.kinnex.${label}.bam"
        gff="${RNA_DIR}/GM12878.kinnex.${label}.collapsed.no5exon.gff"
        abundance="${RNA_DIR}/GM12878.kinnex.${label}.collapsed.no5exon.abundance.txt"
        bed="${RNA_DIR}/GM12878.kinnex.${label}.no5exon.tss.counts.bed.gz"
        tag="$label"
    fi

    echo "collapsing ${tag} isoseq data"
    isoseq collapse --min-aln-coverage 0.85 "$bam" "$gff"

    echo "${tag} metrics for min-aln-coverage 0.85"
    awk '{sum+=$3} END {print "Total reads:", sum}' "$abundance"
    wc -l "$abundance"
    cut -f3 "$abundance" | sort -n | uniq -c | tail -20

    python "${SCRIPT_DIR}/isoseq_collapsed_to_bedgz.py" "$gff" "$abundance" "$bed"
done