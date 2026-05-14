#!/bin/bash
#SBATCH --job-name=preprocess-longread
#SBATCH --account=fc_streets
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A.err

source activate methylseqnet-reproducibility

CONFIG_DIR=$(methylseqnet-repro-path configs)

methylseqnet-preprocess --config "$CONFIG_DIR/preprocess/preprocess_borzoi_longread.gin"
methylseqnet-preprocess --config "$CONFIG_DIR/preprocess/preprocess_borzoi_longread_chrX.gin"