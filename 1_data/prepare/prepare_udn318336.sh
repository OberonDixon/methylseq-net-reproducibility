#!/usr/bin/env bash
set -euo pipefail
RNA_DIR=$(methylseqnet-repro-path rna_tracks)
HAPLOTYPES_DIR=$(methylseqnet-repro-path haplotypes)
GENOMES_DIR=$(methylseqnet-repro-path genomes)
OUT_DIR="$HAPLOTYPES_DIR"

TRIO_FASTA="$HAPLOTYPES_DIR/transloc-asm-v0.20.0-r639.trio.fa.gz"
REFERENCE="$GENOMES_DIR/hg38.fa"
SAMPLE="UDN318336"
THREADS=16
OUTPUT_STEM=${HAPLOTYPES_DIR}/${SAMPLE}.hg38
VCF_FILE=${OUTPUT_STEM}.dip.vcf.gz

echo "=== Trio Assembly to Phased VCF Pipeline (using dipcall) ==="
echo ""

# Create output directories
mkdir -p ${OUT_DIR}/haplotypes

echo "=== Step 1: Extracting haplotypes from trio assembly ==="
# Extract haplotype 1
zcat ${TRIO_FASTA} | \
awk '/^>h1/{p=1} /^>/ && !/^>h1/{p=0} p' | \
sed 's/^>h1/>/' > ${OUT_DIR}/haplotypes/hp1.fa

# Extract haplotype 2
zcat ${TRIO_FASTA} | \
awk '/^>h2/{p=1} /^>/ && !/^>h2/{p=0} p' | \
sed 's/^>h2/>/' > ${OUT_DIR}/haplotypes/hp2.fa

echo "Haplotype 1 contigs: $(grep -c '^>' ${OUT_DIR}/haplotypes/hp1.fa)"
echo "Haplotype 2 contigs: $(grep -c '^>' ${OUT_DIR}/haplotypes/hp2.fa)"

echo ""
echo "=== Step 2: Running dipcall (alignment + variant calling + phasing) ==="
echo "This may take 1-2 hours..."

# Run dipcall - does everything in one command!
run-dipcall \
    -t ${THREADS} \
    "$OUTPUT_STEM" \
    ${REFERENCE} \
    ${OUT_DIR}/haplotypes/hp1.fa \
    ${OUT_DIR}/haplotypes/hp2.fa | make -j 2 -f -

echo ""
echo "=== Pipeline Complete! ==="
echo ""
echo "Output files:"
echo "  Main VCF (all variants):     ${OUT_DIR}/${SAMPLE}.hg38.dip.vcf.gz"
echo "  Confident regions BED:       ${OUT_DIR}/${SAMPLE}.hg38.dip.bed"
echo "  PAF alignment (maternal):    ${OUT_DIR}/${SAMPLE}.hg38.hap1.paf.gz"
echo "  PAF alignment (paternal):    ${OUT_DIR}/${SAMPLE}.hg38.hap2.paf.gz"
echo "  BAM alignment (maternal):    ${OUT_DIR}/${SAMPLE}.hg38.hap1.bam"
echo "  BAM alignment (paternal):    ${OUT_DIR}/${SAMPLE}.hg38.hap2.bam"
echo ""
echo "VCF Statistics:"
bcftools stats ${OUT_DIR}/${SAMPLE}.hg38.dip.vcf.gz | grep "^SN"
echo ""
echo "Variant counts by type:"
bcftools view -H ${OUT_DIR}/${SAMPLE}.hg38.dip.vcf.gz | \
    awk '{print length($4), length($5)}' | \
    awk '{
        if ($1 == 1 && $2 == 1) snp++
        else if ($1 < $2) ins++
        else if ($1 > $2) del++
        else complex++
    }
    END {
        print "SNPs:", snp
        print "Insertions:", ins
        print "Deletions:", del
        print "Complex:", complex
    }'

tabix -f -p vcf ${VCF_FILE}

bcftools consensus -f "$GENOMES_DIR/hg38.fa" -H 1 \
    -i 'TYPE="snp"' ${VCF_FILE} \
    > "${HAPLOTYPES_DIR}/${SAMPLE}.hg38.HP1.fa"
bcftools consensus -f "$GENOMES_DIR/hg38.fa" -H 2 \
    -i 'TYPE="snp"' ${VCF_FILE} \
    > "${HAPLOTYPES_DIR}/${SAMPLE}.hg38.HP2.fa"

python generate_balanced_translocation.py

# REFERENCE="$GENOMES_DIR/hg38.fa"

# # Define VCF and BAM file pairs
# declare -a VCF_FILES=(
#     "${HAPLOTYPES_DIR}/${SAMPLE}.hg38.dip.vcf.gz"
#     "${HAPLOTYPES_DIR}/${SAMPLE}.hg38.dip.vcf.gz"
# )

# declare -a BAM_FILES=(
#     "${RNA_DIR}/PS00289.m84039_230408_220827_s4.RNA.bam"
#     "${RNA_DIR}/PS00449.m84046_240329_064709_s2.IsoSeqX_bc11_5p--IsoSeqX_3p.flnc.bam"
# )

# SCRIPT_DIR="$(dirname "$0")"

# # Loop through paired files
# for i in "${!VCF_FILES[@]}"; do
#     vcf_file="${VCF_FILES[$i]}"
#     bam_file="${BAM_FILES[$i]}"
    
#     # Check if files exist
#     if [[ ! -f "$vcf_file" ]]; then
#         echo "Error: VCF not found: $vcf_file"
#         continue
#     fi
#     if [[ ! -f "$bam_file" ]]; then
#         echo "Error: BAM not found: $bam_file"
#         continue
#     fi
    
#     # Extract output basename (remove .bam extension)
#     bam_basename="${bam_file%.bam}"
#     aligned_bam="${bam_basename}.hg38.bam"

#     echo "Indexing vcf: $vcf_file"
#     tabix -f -p vcf "$vcf_file"
    
#     # Check if BAM needs alignment
#     echo "Checking if BAM is aligned..."
#     if samtools view "$bam_file" | head -1 | cut -f3 | grep -q '^\*$'; then
#         echo "BAM is unaligned, aligning to genome with minimap2..."
#         pbmm2 align \
#             --preset ISOSEQ \
#             --sort \
#             -j 16 \
#             "$REFERENCE" \
#             "$bam_file" \
#             "$aligned_bam"
        
#         echo "Indexing aligned BAM..."
#         samtools index "$aligned_bam"
        
#         bam_file="$aligned_bam"
#         echo "Using aligned BAM: $aligned_bam"
#     else
#         echo "BAM appears to be aligned, using as-is"
#         # Index if not already indexed
#         if [[ ! -f "${bam_file}.bai" ]]; then
#             echo "Indexing BAM..."
#             samtools index "$bam_file"
#         fi
#     fi
    
#     echo "Haplotagging: $bam_file"
    
#     # Run whatshap haplotag and split
#     whatshap haplotag \
#         --reference "$REFERENCE" \
#         --ignore-read-groups \
#         --skip-missing-contigs \
#         "$vcf_file" \
#         "$bam_file" \
#     | samtools view -@ 4 -h \
#     | tee >(samtools view -@ 2 -b -d HP:1 -o "${bam_basename}.HP1.bam") \
#     | samtools view -@ 2 -b -d HP:2 -o "${bam_basename}.HP2.bam"

#     echo "Indexing haplotype BAMs..."
#     samtools index "${bam_basename}.HP1.bam"
#     samtools index "${bam_basename}.HP2.bam"

#     echo "Isoseq collapse"

#     isoseq collapse \
#         --min-aln-coverage 0.85 \
#         "${bam_file}" \
#         "${bam_basename}.collapsed.no5exon.gff"

#     echo "Unphased metrics for min-aln-coverage 0.85"
#     awk '{sum+=$3} END {print "Total reads:", sum}' "${bam_basename}.collapsed.no5exon.abundance.txt"
#     wc -l "${bam_basename}.collapsed.no5exon.abundance.txt"
#     cut -f3 "${bam_basename}.collapsed.no5exon.abundance.txt" | sort -n | uniq -c | tail -20
#     python "${SCRIPT_DIR}/isoseq_collapsed_to_bedgz.py" \
#         "${bam_basename}.collapsed.no5exon.gff" \
#         "${bam_basename}.collapsed.no5exon.abundance.txt" \
#         "${bam_basename}.no5exon.tss.counts.bed.gz"

#     isoseq collapse \
#         --min-aln-coverage 0.85 \
#         "${bam_basename}.HP1.bam" \
#         "${bam_basename}.HP1.collapsed.no5exon.gff"

#     echo "HP1 metrics for min-aln-coverage 0.85"
#     awk '{sum+=$3} END {print "Total reads:", sum}' "${bam_basename}.HP1.collapsed.no5exon.abundance.txt"
#     wc -l "${bam_basename}.HP1.collapsed.no5exon.abundance.txt"
#     cut -f3 "${bam_basename}.HP1.collapsed.no5exon.abundance.txt" | sort -n | uniq -c | tail -20
#     python "${SCRIPT_DIR}/isoseq_collapsed_to_bedgz.py" \
#         "${bam_basename}.HP1.collapsed.no5exon.gff" \
#         "${bam_basename}.HP1.collapsed.no5exon.abundance.txt" \
#         "${bam_basename}.HP1.no5exon.tss.counts.bed.gz"

#     isoseq collapse \
#         --min-aln-coverage 0.85 \
#         "${bam_basename}.HP2.bam" \
#         "${bam_basename}.HP2.collapsed.no5exon.gff"

#     echo "HP2 metrics for min-aln-coverage 0.85"
#     awk '{sum+=$3} END {print "Total reads:", sum}' "${bam_basename}.HP2.collapsed.no5exon.abundance.txt"
#     wc -l "${bam_basename}.HP2.collapsed.no5exon.abundance.txt"
#     cut -f3 "${bam_basename}.HP2.collapsed.no5exon.abundance.txt" | sort -n | uniq -c | tail -20
#     python "${SCRIPT_DIR}/isoseq_collapsed_to_bedgz.py" \
#         "${bam_basename}.HP2.collapsed.no5exon.gff" \
#         "${bam_basename}.HP2.collapsed.no5exon.abundance.txt" \
#         "${bam_basename}.HP2.no5exon.tss.counts.bed.gz"
    
#     echo "Completed: $(basename $bam_file)"
#     echo "---"
# done