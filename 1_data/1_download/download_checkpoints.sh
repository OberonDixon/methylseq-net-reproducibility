#!/usr/bin/env bash
set -euo pipefail

CHECKPOINTS_DIR=$(methylseqnet-repro-path model_checkpoints)

mkdir -p "$CHECKPOINTS_DIR"

ZENODO_DOI="https://zenodo.org/api/records/20275120/files"

# Download files
curl -L "${ZENODO_DOI}/checkpoints.zip/content" -o "${CHECKPOINTS_DIR}/checkpoints.zip"
curl -L "${ZENODO_DOI}/README.md/content" -o "${CHECKPOINTS_DIR}/README.md"

# Unzip contents in place
unzip -d "$CHECKPOINTS_DIR" "${CHECKPOINTS_DIR}/checkpoints.zip"
rm "${CHECKPOINTS_DIR}/checkpoints.zip"