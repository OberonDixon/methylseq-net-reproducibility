#!/bin/bash
#SBATCH --job-name=preprocess-cell-lines
#SBATCH --account=fc_nilah
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A_%a.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A_%a.err
#SBATCH --array=0-15 # 8 folds x 2 configs (phased, unphased); tasks 0-7 are the phased folds and 8-15 the unphased folds

# Command(s) to run:
# Run folds and configs separately on different nodes to speed up preprocessing
module load anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
source activate methylseqnet-reproducibility

DATASET_SPLITS=(
    "fold0"
    "fold1"
    "fold2"
    "fold3"
    "fold4"
    "fold5"
    "fold6"
    "fold7"
)

CONFIG_DIR=$(methylseqnet-repro-path configs)

CONFIG_FILES=(
    "$CONFIG_DIR/preprocess/preprocess_borzoi_longread_cell_lines_phased.gin"
    "$CONFIG_DIR/preprocess/preprocess_borzoi_longread_cell_lines_unphased.gin"
)

# Get the fold and config for this array task
DATASET_SPLIT=${DATASET_SPLITS[$((SLURM_ARRAY_TASK_ID % ${#DATASET_SPLITS[@]}))]}
CONFIG_FILE=${CONFIG_FILES[$((SLURM_ARRAY_TASK_ID / ${#DATASET_SPLITS[@]}))]}

methylseqnet-preprocess --config $CONFIG_FILE --subset $DATASET_SPLIT
