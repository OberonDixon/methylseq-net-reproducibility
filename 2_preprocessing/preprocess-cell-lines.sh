#!/bin/bash
#SBATCH --job-name=preprocess-cell-lines
#SBATCH --account=fc_nilah
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A_%a.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A_%a.err
#SBATCH --array=0-7,16-23 # 8 folds x 3 configs, in the order of CONFIG_FILES below: tasks 0-7 phased unmasked, 8-15 unphased, 16-23 phased masked.
# TODO: this skips the unphased tasks (8-15), which were run already, so a plain sbatch does not reproduce all datasets;
# restore the full range, 0-23, here and in the other scripts like this one together, at the end of the revision work.

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

# The _rnadepth configs also correct the RNA for copy number (rna_depth_correction) and write to separate
# *_rnadepth output directories (the phased ones are the unmasked version, with no masking of poorly phased loci, and the masked version,
# which masks the bins where fewer than half of the reads could be phased). The plain configs without the correction are
# preprocess_borzoi_longread_cell_lines_{phased,unphased}.gin.
CONFIG_FILES=(
    "$CONFIG_DIR/preprocess/preprocess_borzoi_longread_cell_lines_phased_rnadepth_unmasked.gin"
    "$CONFIG_DIR/preprocess/preprocess_borzoi_longread_cell_lines_unphased_rnadepth.gin"
    "$CONFIG_DIR/preprocess/preprocess_borzoi_longread_cell_lines_phased_rnadepth_masked.gin"
)

# Get the fold and config for this array task
DATASET_SPLIT=${DATASET_SPLITS[$((SLURM_ARRAY_TASK_ID % ${#DATASET_SPLITS[@]}))]}
CONFIG_FILE=${CONFIG_FILES[$((SLURM_ARRAY_TASK_ID / ${#DATASET_SPLITS[@]}))]}

methylseqnet-preprocess --config $CONFIG_FILE --subset $DATASET_SPLIT
