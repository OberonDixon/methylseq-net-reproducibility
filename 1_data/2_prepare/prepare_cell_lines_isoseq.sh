#!/usr/bin/env bash
set -euo pipefail

# TSS count beds of the Vollger Iso-Seq BAMs for the cell lines in vollger_cell_lines (GM12878 is handled at higher depth by
# prepare_gm12878.sh). The BAMs already carry whatshap HP tags, so there is no haplotagging here: bam_5prime_tss_beds.py counts
# the 5' ends of the reads of each BAM into the beds of all reads and of the HP1 and HP2 reads (see its docstring).
#
# A sample whose three beds exist is skipped, so the script can be re-run to resume. Environment:
#   ISOSEQ_DIR      the Iso-Seq dataset directory
#   BED_DIR         where the beds go, default ${ISOSEQ_DIR}/beds_5prime
#   FORCE_BEDS=1    redo the beds even if they exist
#   PYTHON          the python to use, which needs pysam (default: python)
ISOSEQ_DIR="${ISOSEQ_DIR:-/global/scratch/projects/vector_streetslab/oberon/datasets/vollger_isoseq}"
BAM_DIR="${ISOSEQ_DIR}/bams"
BED_DIR="${BED_DIR:-${ISOSEQ_DIR}/beds_5prime}"
LOG_DIR="${ISOSEQ_DIR}/logs"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-python}"

# Samples processed concurrently, and BAM decompression threads per sample
JOBS="${JOBS:-4}"
THREADS="${THREADS:-4}"

mkdir -p "$BED_DIR" "$LOG_DIR"

"$PYTHON" -c 'import pysam' 2>/dev/null || {
    echo "Error: ${PYTHON} cannot import pysam; activate the methylseqnet-reproducibility environment or set PYTHON to its python" >&2
    exit 1
}

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
    local bam="${BAM_DIR}/${name}.tagged.bam"
    local unphased_bed="${BED_DIR}/${name}.tss.counts.bed.gz"
    local hp1_bed="${BED_DIR}/${name}.HP1.tss.counts.bed.gz"
    local hp2_bed="${BED_DIR}/${name}.HP2.tss.counts.bed.gz"
    [[ -f "$bam" ]] || { echo "Error: $bam not found." >&2; return 1; }

    if [[ -s "${unphased_bed}.tbi" && -s "${hp1_bed}.tbi" && -s "${hp2_bed}.tbi" && -z "${FORCE_BEDS:-}" ]]; then
        echo "the beds of ${name} already exist, skipping (FORCE_BEDS=1 redoes them)"
        return 0
    fi
    "$PYTHON" "${SCRIPT_DIR}/bam_5prime_tss_beds.py" "$bam" \
        --unphased-output "$unphased_bed" --hp1-output "$hp1_bed" --hp2-output "$hp2_bed" --threads "$THREADS"
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
