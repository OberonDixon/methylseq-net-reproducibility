#!/usr/bin/env bash
set -euo pipefail

# Pre-haplotagged Vollger Iso-Seq BAMs for the cell lines in vollger_cell_lines (GM12878 is
# handled at higher depth by prepare_gm12878.sh). The BAMs already carry whatshap HP tags, so
# there is no haplotagging here: ONE isoseq collapse of all reads -> tss counts bedgz of all reads
# and of the HP1 and HP2 reads. The reads are attributed to haplotypes after the collapse, using
# the collapse's read_stat.txt and the BAM's HP tags (see isoseq_tss_beds.sh), so that the unphased
# and haplotype beds share their TSS positions. (Collapsing the HP1 and HP2 reads separately, as
# this script used to, moves the TSS of the same reads between the beds.)
#
# Each step is skipped when its outputs exist, so the script can be re-run to resume. Environment:
#   ISOSEQ_DIR      the Iso-Seq dataset directory
#   BED_DIR         where all three beds of a sample (unphased, HP1, HP2) go, default ${ISOSEQ_DIR}/beds_single_collapse.
#                   The older ${ISOSEQ_DIR}/beds holds beds from separate collapses of the HP1, HP2 and all reads, whose TSS
#                   positions do not match; it is left alone, and must not be used as BED_DIR because its existing files would
#                   make the steps skip
#   FORCE_COLLAPSE=1  redo isoseq collapse even if complete;  FORCE_BEDS=1  redo the beds, reusing the collapse
ISOSEQ_DIR="${ISOSEQ_DIR:-/global/scratch/projects/vector_streetslab/oberon/datasets/vollger_isoseq}"
BAM_DIR="${ISOSEQ_DIR}/bams"
INTERMEDIATE_DIR="${ISOSEQ_DIR}/intermediate"
BED_DIR="${BED_DIR:-${ISOSEQ_DIR}/beds_single_collapse}"
LOG_DIR="${INTERMEDIATE_DIR}/logs"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Samples processed concurrently, and threads per sample (JOBS * THREADS ~ cores on the node)
JOBS="${JOBS:-4}"
THREADS="${THREADS:-8}"

mkdir -p "$INTERMEDIATE_DIR" "$BED_DIR" "$LOG_DIR"

source "${SCRIPT_DIR}/isoseq_tss_beds.sh"
check_isoseq_requirements || exit 1

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
    [[ -f "$bam" ]] || { echo "Error: $bam not found." >&2; return 1; }

    make_tss_beds "$bam" \
        "${INTERMEDIATE_DIR}/${name}.collapsed.no5exon" \
        "${BED_DIR}/${name}.no5exon.tss.counts.bed.gz" \
        "${BED_DIR}/${name}.HP1.no5exon.tss.counts.bed.gz" \
        "${BED_DIR}/${name}.HP2.no5exon.tss.counts.bed.gz"
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
