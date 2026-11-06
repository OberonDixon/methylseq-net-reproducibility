#!/usr/bin/env bash
set -euo pipefail

# Download ENCODE SCREEN (V4, GRCh38) candidate cis-regulatory elements:
#   GRCh38-PLS.bed : promoter-like signatures
#   GRCh38-ELS.bed : enhancer-like signatures (column 6 = pELS / dELS)

ANNOTATIONS_DIR=$(methylseqnet-repro-path annotations)
mkdir -p "$ANNOTATIONS_DIR"

for CLASS in PLS ELS; do
    wget -O "$ANNOTATIONS_DIR/GRCh38-${CLASS}.bed" \
        "https://downloads.wenglab.org/cCREs/GRCh38-${CLASS}.bed"
done
