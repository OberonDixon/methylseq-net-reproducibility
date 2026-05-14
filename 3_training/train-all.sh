#!/bin/bash
#SBATCH --job-name=train_methylseqnet_factorized
#SBATCH --account=fc_nilah
#SBATCH --partition=savio3_gpu
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --requeue
#SBATCH --time=72:00:00
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/train_methylseqnet_%A_%a.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/train_methylseqnet_%A_%a.err
#SBATCH --exclude=n0217.savio3,n0215.savio3,n0130.savio4,n0132.savio4,n0134.savio3,n0135.savio3,n0136.savio3,n0137.savio3,n0138.savio3,n0143.savio3,n0144.savio3,n0145.savio3,n0158.savio3,n0159.savio3,n0160.savio3,n0161.savio3,n0174.savio3,n0175.savio3,n0176.savio3
#SBATCH --array=36-36 # Specify the range of array jobs (e.g., 0-2 for 3 configs)
# Command(s) to run:
# Define an array of config files
source activate methylseqnet-reproducibility
CONFIG_DIR=$(methylseqnet-repro-path configs)
CONFIG_FILES=(
    "$CONFIG_DIR/train/borzoi-rep0_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized_atlas+longread_imputed.gin"
    "$CONFIG_DIR/train/borzoi-rep0_concat_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_FiLM_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_probe_atlas+longread.gin"
    "$CONFIG_DIR/train/basenji2-reinit_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/basenji2-reinit_factorized_atlas+longread_imputed.gin"
    "$CONFIG_DIR/train/borzoi-rep1_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized_atlas_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized_atlas_imputed.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized1k_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized1k_atlas+longread_imputed.gin"
    "$CONFIG_DIR/train/basenji2_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/basenji2_factorized_atlas+longread_imputed.gin"
    "$CONFIG_DIR/train/borzoi-rep2_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep3_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized128_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized256_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-transformer4x256_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-transformer8x128_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-dilate1_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-dilate3_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-dilate7_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-cond1to15_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-cond1to3_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-cond3to1_atlas+longread_true.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized-cond15to1_atlas+longread_true.gin"
    "$CONFIG_DIR/train/bassetlike_methylseq-only_atlas+longread_onehot.gin"
    "$CONFIG_DIR/train/bassetlike_methylseq-only_atlas+longread_interp.gin"
    "$CONFIG_DIR/train/bassetlike_methylseq-only_atlas+longread_smoothed.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized_longread_imputed.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized_longread_true.gin"
    "$CONFIG_DIR/train/bassetlike_methylseq-only_atlas+longread_smoothed1k.gin"
    "$CONFIG_DIR/train/bassetlike_methylseq-only_atlas+longread_smoothed8k.gin"
    "$CONFIG_DIR/train/bassetlike_seq-only_atlas+longread.gin"
    "$CONFIG_DIR/train/bassetlike_methyl-only_atlas+longread_smoothed.gin"
    "$CONFIG_DIR/train/borzoi-rep0_factorized_atlas+longread_unimputed.gin"
)
CONFIG_FILE=${CONFIG_FILES[$SLURM_ARRAY_TASK_ID]}
UNIQUE_IDENTIFIER="slurm${SLURM_ARRAY_JOB_ID}task${SLURM_ARRAY_TASK_ID}"
methylseqnet-train --config $CONFIG_FILE --unique_identifier $UNIQUE_IDENTIFIER --batch_size 1
