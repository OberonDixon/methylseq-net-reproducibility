#!/usr/bin/env bash
# Download GEO accession GSE186458 (methylation atlas) and decompress
# Requires the methylseqnet-reproducibility conda environment for path resolution.
#
# Usage: ./download_methylation_atlas.sh
 
set -euo pipefail
 
# Resolve output directory
METHYLATION_DIR=$(methylseqnet-repro-path methylation_atlas)
echo "Methylation atlas directory: $METHYLATION_DIR"
mkdir -p "$METHYLATION_DIR"
 
# Download
TARBALL="$METHYLATION_DIR/GSE186458_RAW.tar"
echo "Downloading GSE186458 from GEO..."
wget --content-disposition \
     --output-document="$TARBALL" \
     "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE186458&format=file"
 
# Decompress
EXTRACT_DIR="$METHYLATION_DIR/GSE186458_RAW"
echo "Extracting to $EXTRACT_DIR ..."
mkdir -p "$EXTRACT_DIR"
tar -xf "$TARBALL" -C "$EXTRACT_DIR"
 
echo ""
echo "Done. Files in $EXTRACT_DIR:"
ls "$EXTRACT_DIR"