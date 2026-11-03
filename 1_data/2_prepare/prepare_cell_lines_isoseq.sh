#!/usr/bin/env bash
set -euo pipefail

# Pre-haplotagged Vollger Iso-Seq BAMs for the cell lines in vollger_cell_lines (GM12878 is
# handled at higher depth by prepare_gm12878.sh). The BAMs already carry whatshap HP tags, so
# there is no haplotagging here: split by HP tag -> isoseq collapse -> tss counts bedgz.
ISOSEQ_DIR="${ISOSEQ_DIR:-/global/scratch/projects/vector_streetslab/oberon/datasets/vollger_isoseq}"
BAM_DIR="${ISOSEQ_DIR}/bams"
INTERMEDIATE_DIR="${ISOSEQ_DIR}/intermediate"
BED_DIR="${ISOSEQ_DIR}/beds"
LOG_DIR="${INTERMEDIATE_DIR}/logs"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Samples processed concurrently, and threads per sample (JOBS * THREADS ~ cores on the node)
JOBS="${JOBS:-4}"
THREADS="${THREADS:-8}"

mkdir -p "$INTERMEDIATE_DIR" "$BED_DIR" "$LOG_DIR"

# Iso-Seq sample IDs, one per cell line (fiber-seq ID -> iso-seq ID, from links.txt):
#   HepG2 PS01307 -> PS01309, Panc1 PS01316 -> PS01318, Caco-2 PS01319 -> PS01321,
#   THP-1 PS01326 -> PS01327, Jurkat PS01329 -> PS01330, Hek293 PS01332 -> PS01333,
#   Hap1 PS01367 -> PS01368, K562 PS01370 -> PS01371, H1 PS01373 -> PS01374,
#   H9 PS01376 -> PS01377
SAMPLES=(
    PS01309 PS01318 PS01321 PS01327 PS01330
    PS01333 PS01368 PS01371 PS01374 PS01377
)

process_sample() {
    local sample="$1"
    local name="ENCODE-${sample}_untreated_${sample}"
    local full_bam="${BAM_DIR}/${name}.tagged.bam"
    [[ -f "$full_bam" ]] || { echo "Error: $full_bam not found." >&2; return 1; }

    echo "splitting ${sample} by HP tag"
    # reads without an HP tag are in neither HP bam
    local hp
    for hp in 1 2; do
        local hp_bam="${INTERMEDIATE_DIR}/${name}.HP${hp}.bam"
        if [[ ! -s "$hp_bam" ]]; then
            samtools view -@ "$THREADS" -b -d "HP:${hp}" -o "${hp_bam}.part" "$full_bam"
            mv "${hp_bam}.part" "$hp_bam"
        fi
        samtools index -@ "$THREADS" "$hp_bam"
    done
    samtools index -@ "$THREADS" "$full_bam"

    echo "isoseq collapse -> bedgz tss counts"
    local label bam stem bed tag gff abundance
    for label in "" "HP1" "HP2"; do
        if [[ -z "$label" ]]; then
            bam="$full_bam"
            stem="${INTERMEDIATE_DIR}/${name}.collapsed.no5exon"
            bed="${BED_DIR}/${name}.no5exon.tss.counts.bed.gz"
            tag="unphased"
        else
            bam="${INTERMEDIATE_DIR}/${name}.${label}.bam"
            stem="${INTERMEDIATE_DIR}/${name}.${label}.collapsed.no5exon"
            bed="${BED_DIR}/${name}.${label}.no5exon.tss.counts.bed.gz"
            tag="$label"
        fi
        gff="${stem}.gff"
        abundance="${stem}.abundance.txt"

        echo "collapsing ${sample} ${tag} isoseq data"
        isoseq collapse -j "$THREADS" --min-aln-coverage 0.85 "$bam" "$gff"

        echo "${sample} ${tag} metrics for min-aln-coverage 0.85"
        awk '{sum+=$3} END {print "Total reads:", sum}' "$abundance"
        wc -l "$abundance"

        python "${SCRIPT_DIR}/isoseq_collapsed_to_bedgz.py" "$gff" "$abundance" "$bed"
    done
}

# Job pool: up to $JOBS samples at once, each logging to its own file. Each job records its
# exit status in a .status file (bash 4.4 has no `wait -p`), and the pool is throttled with
# `wait -n`.
run_sample() {
    local sample="$1" rc
    rm -f "${LOG_DIR}/${sample}.status"
    # plain subshell (not in a ||/if list) so `set -e` stops the sample at its first failure
    set +e
    ( set -e; process_sample "$sample" ) > "${LOG_DIR}/${sample}.log" 2>&1
    rc=$?
    set -e
    echo "$rc" > "${LOG_DIR}/${sample}.status"
    return "$rc"
}

for sample in "${SAMPLES[@]}"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || true; done
    echo "START ${sample} (log: ${LOG_DIR}/${sample}.log)"
    run_sample "$sample" &
done
wait || true

fail=0
for sample in "${SAMPLES[@]}"; do
    if [[ "$(cat "${LOG_DIR}/${sample}.status" 2>/dev/null)" == "0" ]]; then
        echo "DONE ${sample}"
    else
        echo "FAIL ${sample} (see ${LOG_DIR}/${sample}.log)" >&2
        fail=1
    fi
done

if (( fail )); then
    echo "Some samples failed; re-run to resume." >&2
    exit 1
fi
echo "All samples complete."
