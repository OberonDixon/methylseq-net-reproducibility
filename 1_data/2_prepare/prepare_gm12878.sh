#!/usr/bin/env bash
set -euo pipefail

# GM12878 (NA12878) Kinnex Iso-Seq. The unoriented BAM from the Vollger et al. download is haplotagged with whatshap and oriented
# by its segment-adapter tags (orient_bam_dldr.py explains why); bam_5prime_tss_beds.py then counts the 5' ends of its reads into
# the TSS count beds of all reads and of the HP1 and HP2 reads (see its docstring). The same final step as
# prepare_cell_lines_isoseq.sh for the other cell lines, whose BAMs come already haplotagged; here the haplotagged, oriented BAM
# has to be made first.
#
# Each step is skipped when its outputs exist, so the script can be re-run to resume. The haplotagged BAM (and so the phased VCF
# and whatshap) is only prepared when something needs it: nothing is done, and no VCF is needed, when the three beds are already
# complete. Environment:
#   ORIENTED_DIR     where the tagged, oriented BAM goes; default: the directory next to the rna_tracks directory called
#                    rna_bams_dldr_oriented
#   BED_DIR          where the three beds go, default ${ORIENTED_DIR}/beds_5prime
#   VCF_FILE         the phased NA12878 VCF, default ${HAPLOTYPES_DIR}/NA12878.vcf.gz. It needs its .tbi index
#   FORCE_TAG=1      redo the haplotagging and orientation, and then the beds, as the HP tags change
#   FORCE_BEDS=1     redo the beds, reusing the tagged BAM
#   PYTHON           the python to use, which needs pysam (default: python)
# The HP tags depend on the whatshap version (2.8 tags about 6% more alignments than 2.4 on chr22, with the same haplotype for
# every alignment that both tag), so the environment's whatshap version should be recorded with the results.
RNA_DIR=$(methylseqnet-repro-path rna_tracks)
HAPLOTYPES_DIR=$(methylseqnet-repro-path haplotypes)
GENOMES_DIR=$(methylseqnet-repro-path genomes)
VCF_FILE="${VCF_FILE:-$HAPLOTYPES_DIR/NA12878.vcf.gz}"
ORIENTED_DIR="${ORIENTED_DIR:-$(dirname "$RNA_DIR")/rna_bams_dldr_oriented}"
BED_DIR="${BED_DIR:-${ORIENTED_DIR}/beds_5prime}"
THREADS="${THREADS:-8}"
PYTHON="${PYTHON:-python}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$ORIENTED_DIR" "$BED_DIR"

# Exit with an error naming the tools that are not on the PATH
require_tools() {
    local missing=0 tool
    for tool in "$@"; do
        command -v "$tool" >/dev/null 2>&1 || { echo "Error: ${tool} is not on PATH" >&2; missing=1; }
    done
    return "$missing"
}
"$PYTHON" -c 'import pysam' 2>/dev/null || {
    echo "Error: ${PYTHON} cannot import pysam; activate the methylseqnet-reproducibility environment or set PYTHON to its python" >&2
    exit 1
}

require_vcf() {
    [[ -f "$VCF_FILE" ]] || { echo "Error: the phased NA12878 VCF is needed but not found at $VCF_FILE (set VCF_FILE)." >&2; exit 1; }
}

# The SNP-only haplotype genomes (used as the sequences of the haplotyped datasets), skipped when they already exist
for hp in 1 2; do
    haplotype_fasta="${HAPLOTYPES_DIR}/NA12878.hg38.hp${hp}.fa"
    if [[ -s "$haplotype_fasta" ]]; then
        echo "${haplotype_fasta} already exists, skipping"
    else
        require_vcf
        require_tools bcftools || exit 1
        bcftools consensus -f "$GENOMES_DIR/hg38.fa" -H "$hp" \
            -i 'TYPE="snp"' "${VCF_FILE}" \
            > "${haplotype_fasta}.tmp"
        mv "${haplotype_fasta}.tmp" "$haplotype_fasta"
    fi
done

TAGGED_BAM="${ORIENTED_DIR}/GM12878.kinnex.tagged.bam"
UNPHASED_BED="${BED_DIR}/GM12878.kinnex.tss.counts.bed.gz"
HP1_BED="${BED_DIR}/GM12878.kinnex.HP1.tss.counts.bed.gz"
HP2_BED="${BED_DIR}/GM12878.kinnex.HP2.tss.counts.bed.gz"

# New HP tags mean new beds
[[ -z "${FORCE_TAG:-}" ]] || export FORCE_BEDS=1

if [[ -s "${UNPHASED_BED}.tbi" && -s "${HP1_BED}.tbi" && -s "${HP2_BED}.tbi" && -z "${FORCE_BEDS:-}" ]]; then
    echo "the three beds are already complete, nothing to do (FORCE_BEDS=1 or FORCE_TAG=1 redo them)"
    exit 0
fi

echo "haplotagging and orienting isoseq data"
echo "----------------------------------------------"

# One BAM of all reads with the HP tags (reads whatshap cannot phase just have none), with the strand of the reads corrected.
# Haplotagging does not depend on the strand flag, so the order of the two steps does not matter.
if [[ -s "$TAGGED_BAM" && -s "${TAGGED_BAM}.bai" && -z "${FORCE_TAG:-}" ]]; then
    echo "${TAGGED_BAM} already exists, skipping"
else
    require_vcf
    require_tools whatshap samtools || exit 1
    HAPLOTAGGED_TMP="${TAGGED_BAM%.bam}.haplotagged.tmp.bam"
    ORIENTED_TMP="${TAGGED_BAM%.bam}.oriented.tmp.bam"
    [[ -f "$RNA_DIR/GM12878.kinnex.bam.bai" ]] || samtools index -@ "$THREADS" "$RNA_DIR/GM12878.kinnex.bam"

    echo "running whatshap $(whatshap --version)"
    whatshap haplotag \
        --reference "$GENOMES_DIR/hg38.fa" \
        --ignore-read-groups \
        --skip-missing-contigs \
        -o "$HAPLOTAGGED_TMP" \
        "$VCF_FILE" \
        "$RNA_DIR/GM12878.kinnex.bam"

    echo "orienting reads by their segment adapters"
    "$PYTHON" "${SCRIPT_DIR}/orient_bam_dldr.py" "$HAPLOTAGGED_TMP" "$ORIENTED_TMP"   # also writes the index
    mv "$ORIENTED_TMP" "$TAGGED_BAM"
    mv "${ORIENTED_TMP}.bai" "${TAGGED_BAM}.bai"
    rm -f "$HAPLOTAGGED_TMP" "${HAPLOTAGGED_TMP}.bai"
fi

echo "5' end counts -> bedgz tss counts"
echo "----------------------------------------------"

"$PYTHON" "${SCRIPT_DIR}/bam_5prime_tss_beds.py" "$TAGGED_BAM" \
    --unphased-output "$UNPHASED_BED" --hp1-output "$HP1_BED" --hp2-output "$HP2_BED" --threads "$THREADS"
