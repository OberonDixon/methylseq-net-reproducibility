#!/usr/bin/env bash
# Shared steps for turning a haplotagged Iso-Seq BAM into TSS count beds. Sourced by the prepare_* scripts; not run directly.
#
# The caller must set SCRIPT_DIR (the directory of the helper python scripts, i.e. this one) and THREADS. The helper scripts
# run with $PYTHON (default: python), which needs pysam: activate the methylseqnet-reproducibility environment, or set PYTHON.
#
#   make_tss_beds BAM STEM UNPHASED_BED HP1_BED HP2_BED [HP_TAG_BAMS]
#
#     BAM           aligned, sorted BAM of all reads, with the whatshap HP tags (reads without a tag are fine)
#     STEM          prefix of the `isoseq collapse` outputs: STEM.gff, STEM.abundance.txt, STEM.read_stat.txt, STEM.report.json, ...
#     UNPHASED_BED  TSS counts of all the reads of the collapse
#     HP1_BED, HP2_BED  the same isoforms with only the reads of that haplotype, so that all three beds have the same TSS positions
#     HP_TAG_BAMS   optional, space-separated BAMs to take the HP tags from instead of BAM, e.g. the HP1 and HP2 BAMs of an
#                   earlier haplotagging, so that it need not be redone. BAM is then only used if the collapse has to be (re)done.
#
# The three steps are independent and each is skipped when its outputs already exist, so completed work is never redone:
#   1. one isoseq collapse of ALL reads                                  (FORCE_COLLAPSE=1 redoes it)
#   2. UNPHASED_BED from the collapse                                    (FORCE_BEDS=1 redoes it)
#   3. HP1_BED and HP2_BED from the collapse, read_stat.txt and the BAM's HP tags (FORCE_BEDS=1 redoes them)
# To redo only the beds after changing the bed code, set FORCE_BEDS=1: the existing collapse is reused.
# Whenever the collapse is (re)done in a call, the beds are redone too, since existing beds would belong to the old collapse.
# Note that redoing the beds overwrites UNPHASED_BED in place, so do not do it while anything is reading that file.

PYTHON="${PYTHON:-python}"

# Fail early, before any long step, if tools the steps need are missing.
require_tools() {
    local missing=0 tool
    for tool in "$@"; do
        command -v "$tool" >/dev/null 2>&1 || { echo "Error: ${tool} is not on PATH" >&2; missing=1; }
    done
    return "$missing"
}

# The tools of the collapse and bed steps; extra tools can be passed as arguments.
check_isoseq_requirements() {
    local missing=0
    require_tools isoseq bgzip tabix "$@" || missing=1
    "$PYTHON" -c 'import pysam' 2>/dev/null || {
        echo "Error: ${PYTHON} ($(command -v "$PYTHON" || echo 'not found')) cannot import pysam; activate the methylseqnet-reproducibility environment or set PYTHON to its python" >&2
        missing=1
    }
    return "$missing"
}

collapse_is_complete() {
    local stem="$1"
    [[ -s "${stem}.gff" && -s "${stem}.abundance.txt" && -s "${stem}.read_stat.txt" && -s "${stem}.report.json" ]]
}

# True when make_tss_beds would find nothing to do: the collapse is complete and the three beds exist, and neither FORCE_COLLAPSE
# nor FORCE_BEDS is set. A caller can use it to skip preparing the inputs (such as a haplotagged BAM) that only the steps need.
#   tss_beds_up_to_date STEM UNPHASED_BED HP1_BED HP2_BED
tss_beds_up_to_date() {
    collapse_is_complete "$1" && [[ -s "$2" && -s "$3" && -s "$4" ]] && [[ -z "${FORCE_COLLAPSE:-}" && -z "${FORCE_BEDS:-}" ]]
}

make_tss_beds() {
    local bam="$1" stem="$2" unphased_bed="$3" hp1_bed="$4" hp2_bed="$5" hp_tag_bams="${6:-$1}"
    local -a hp_tag_bam_list
    read -r -a hp_tag_bam_list <<< "$hp_tag_bams"
    local gff="${stem}.gff" abundance="${stem}.abundance.txt" read_stat="${stem}.read_stat.txt"
    local redo_beds="${FORCE_BEDS:-}"

    if collapse_is_complete "$stem" && [[ -z "${FORCE_COLLAPSE:-}" ]]; then
        echo "isoseq collapse: ${stem} is already complete, skipping"
    else
        echo "isoseq collapse of all reads in ${bam}"
        [[ -f "${bam}.bai" ]] || samtools index -@ "$THREADS" "$bam"
        isoseq collapse -j "$THREADS" --min-aln-coverage 0.85 "$bam" "$gff"
        collapse_is_complete "$stem" || { echo "Error: isoseq collapse did not produce all of its outputs for ${stem}" >&2; return 1; }
        redo_beds=1
    fi
    echo "metrics for min-aln-coverage 0.85"
    awk '{sum+=$3} END {print "Total reads:", sum}' "$abundance"
    wc -l "$abundance"

    if [[ -s "$unphased_bed" && -z "$redo_beds" ]]; then
        echo "unphased TSS bed: ${unphased_bed} already exists, skipping"
    else
        mkdir -p "$(dirname "$unphased_bed")"
        "$PYTHON" "${SCRIPT_DIR}/isoseq_collapsed_to_bedgz.py" "$gff" "$abundance" "$unphased_bed"
    fi

    if [[ -s "$hp1_bed" && -s "$hp2_bed" && -z "$redo_beds" ]]; then
        echo "haplotype TSS beds: ${hp1_bed} and ${hp2_bed} already exist, skipping"
    else
        mkdir -p "$(dirname "$hp1_bed")" "$(dirname "$hp2_bed")"
        "$PYTHON" "${SCRIPT_DIR}/isoseq_hp_counts_to_bedgz.py" \
            --gff "$gff" --abundance "$abundance" --read-stat "$read_stat" --bam "${hp_tag_bam_list[@]}" \
            --hp1-output "$hp1_bed" --hp2-output "$hp2_bed" --threads "$THREADS"
    fi
}
