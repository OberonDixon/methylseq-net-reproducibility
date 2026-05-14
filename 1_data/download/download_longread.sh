#!/usr/bin/env bash
# Download data from https://www.biorxiv.org/content/10.1101/2023.09.26.559521v1
# Requires the methylseqnet-reproducibility conda environment for path resolution.
#
# Usage: ./download_longread.sh
 
set -euo pipefail
 
# Resolve output directories
PACBIO_5MC_DIR=$(methylseqnet-repro-path pacbio_5mC_tracks)
echo "Methylation atlas directory: $PACBIO_5MC_DIR"
mkdir -p "$PACBIO_5MC_DIR"

FIBERSEQ_DIR=$(methylseqnet-repro-path fiberseq_tracks)
FIBERSEQ_BW_DIR="$FIBERSEQ_DIR/GM12878_trackHub/bw"
echo "Fiber-seq bigwigs directory: $FIBERSEQ_BW_DIR"
mkdir -p "$FIBERSEQ_BW_DIR"

RNA_DIR=$(methylseqnet-repro-path rna_tracks)
echo "RNA-seq directory: $RNA_DIR"
mkdir -p "$RNA_DIR"

BASE_URL="https://s3-us-west-1.amazonaws.com/stergachis-manuscript-data/2023/Vollger_et_al_long-read_multi-ome"
# Download PacBio 5mC bigwigs (hap1 and hap2)
wget -P "$PACBIO_5MC_DIR" "${BASE_URL}/5mC/GM12878_WGS-pb-5mC.hap1.bw"
wget -P "$PACBIO_5MC_DIR" "${BASE_URL}/5mC/GM12878_WGS-pb-5mC.hap2.bw"
# Download Fiber-seq bigwigs (hap1 and hap2)
wget -P "$FIBERSEQ_BW_DIR" "${BASE_URL}/GM12878_pacbiome/trackHub/bw/hap1.acc.bw"
wget -P "$FIBERSEQ_BW_DIR" "${BASE_URL}/GM12878_pacbiome/trackHub/bw/hap2.acc.bw"
# Download RNA-seq bam
wget -O "$RNA_DIR/GM12878.kinnex.bam" "${BASE_URL}/isoseq/GM12878.isoseq.bam"