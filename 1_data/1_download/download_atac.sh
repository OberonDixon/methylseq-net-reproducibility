#!/usr/bin/env bash
# Download all CATLAS (Zhang et al. 2021 Cell) per-cell-type ATAC-seq bigwigs
# Source: https://decoder-genetics.wustl.edu/catlasv1/humanenhancer/data/bw/

# Usage: ./download_catlas_bigwigs.sh [output_dir]
#   output_dir defaults to ./bigwigs

ATAC_DIR=$(methylseqnet-repro-path atac_atlas)/bigwigs
mkdir -p "$ATAC_DIR"

set -euo pipefail

BASE_URL="https://decoder-genetics.wustl.edu/catlasv1/humanenhancer/data/bw"
OUTDIR="$ATAC_DIR"

mkdir -p "$OUTDIR"

echo "Fetching file list from $BASE_URL ..."
LISTING=$(curl -fsSL "$BASE_URL/")

mapfile -t SERVER_FILES < <(echo "$LISTING" | grep -oP 'href="\K[^"/]+\.bw(?=")' | sort -u)

N=${#SERVER_FILES[@]}
if [[ $N -eq 0 ]]; then
    echo "ERROR: No .bw files found in directory listing. Check the URL." >&2
    exit 1
fi

echo "Found $N bigwig files. Downloading to $OUTDIR/ ..."
echo ""

COUNT=0
FAILED=()

for server_fname in "${SERVER_FILES[@]}"; do
    # Convert underscores back to %20 to match original local filenames
    local_fname="${server_fname//_/%20}"
    outfile="$OUTDIR/$local_fname"

    if [[ -f "$outfile" ]]; then
        echo "  [skip] $local_fname"
        COUNT=$((COUNT + 1))
        continue
    fi

    echo "  $local_fname"
    if curl -fsSL "$BASE_URL/$server_fname" -o "$outfile"; then
        COUNT=$((COUNT + 1))
    else
        echo "  [FAILED] $local_fname"
        FAILED+=("$local_fname")
        rm -f "$outfile"
    fi
done

echo ""
echo "Done: $COUNT / $N files in $OUTDIR/"

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo "Failed (${#FAILED[@]}):"
    for f in "${FAILED[@]}"; do echo "  $f"; done
    exit 1
fi