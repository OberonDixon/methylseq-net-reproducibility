#!/bin/bash

# source activate methylseqnet-reproducibility

PREPROCESSED=$(methylseqnet-repro-path preprocessed_datasets)
CONFIG_DIR=$(methylseqnet-repro-path configs)
ANNOTATIONS_DIR=$(methylseqnet-repro-path annotations)

methylseqnet-peaks --dataset-paths \
    $PREPROCESSED/atlas/fold0.h5 \
    $PREPROCESSED/atlas/fold1.h5 \
    $PREPROCESSED/atlas/fold2.h5 \
    $PREPROCESSED/atlas/fold5.h5 \
    $PREPROCESSED/atlas/fold6.h5 \
    $PREPROCESSED/atlas/fold7.h5 \
    --output-directory $PREPROCESSED/motif_insertion/peaks/random_celltype_peaks_cpg.05-.15_10cts \
    --label-substrings "" \
    --data-type "ATAC-seq" \
    --peak-threshold 10 \
    --num-peaks 500 \
    --min-peak-distance 131072 \
    --cpg-density-range 0.05 0.15 \
    --cpg-density-window 2048 \
    --random-seeds 1

methylseqnet-insert-motifs \
    $ANNOTATIONS_DIR \
    $PREPROCESSED/motif_insertion/peaks/random_celltype_peaks_cpg.05-.15_10cts \
    $PREPROCESSED/motif_insertion/fastas/motif_inserted_random_celltype_cpg.05-.15_10cts_peaks_2048 \
    --input-len 16384 \
    --shuffle-len 2048 \
    --tfs-file "${CONFIG_DIR}/transcription_factors.txt" \
    --overwrite

methylseqnet-preprocess --config "$CONFIG_DIR/preprocess/preprocess_motif_insert_alltypes_cpg.05-.15_peaks_10cts.gin"