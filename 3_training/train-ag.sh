#!/bin/bash
#SBATCH --job-name=train_methylseqnet_ag
#SBATCH --account=fc_nilah
#
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/train_methylseqnet_ag_%A_%a.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/train_methylseqnet_ag_%A_%a.err
#
#SBATCH --partition=savio3_gpu
## SBATCH --qos=savio_lowprio
#SBATCH --qos=a40_gpu3_normal
#SBATCH --exclude=n0214.savio3,n0215.savio3,n0211.savio3
#
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:2
#
#SBATCH --time=72:00:00
#SBATCH --requeue
#SBATCH --array=4-5
#
module load anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /global/scratch/projects/vector_streetslab/oberon/envs/methylseqnet-ag

# Fail fast if GPU is not actually usable — avoids silently running 72h on CPU
python3 -c "import torch; assert torch.cuda.is_available(), 'CUDA unavailable'" || { echo "ERROR: No GPU detected on $(hostname), exiting so SLURM can requeue to a working node."; exit 1; }

CONFIG_DIR=$(methylseqnet-repro-path configs)
CONFIG_FILES=(
    "$CONFIG_DIR/train/ag-rep0-128bp_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/ag-rep0-128bp_factorized_atlas+longread_imputed.gin"
    "$CONFIG_DIR/train/ag-rep0-1bp_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/ag-rep0-1bp_factorized_atlas+longread_imputed.gin"
    "$CONFIG_DIR/train/ag-rep1-128bp_factorized_atlas+longread_true.gin"
    "$CONFIG_DIR/train/ag-rep1-128bp_factorized_atlas+longread_imputed.gin"
)
CONFIG_FILE=${CONFIG_FILES[$SLURM_ARRAY_TASK_ID]}
UNIQUE_IDENTIFIER="slurm${SLURM_ARRAY_JOB_ID}task${SLURM_ARRAY_TASK_ID}"
methylseqnet-train --config $CONFIG_FILE --unique-identifier $UNIQUE_IDENTIFIER --batch-size 1