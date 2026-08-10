#!/bin/bash
#SBATCH --job-name=train_methylseqnet_ag
#SBATCH --account=fc_nilah
#
#SBATCH --output=/clusterfs/nilah/ayesha/lightning_logs/ag-methylseqnet/train_methylseqnet_ag_%j.out
#SBATCH --error=/clusterfs/nilah/ayesha/lightning_logs/ag-methylseqnet/train_methylseqnet_ag_%j.err
#
#SBATCH --partition=savio3_gpu
#SBATCH --qos=a40_gpu3_normal
#SBATCH --exclude=n0214.savio3,n0215.savio3
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#
#SBATCH --time=72:00:00
#SBATCH --requeue
#
#SBATCH --mail-user=arbajwa@berkeley.edu
#SBATCH --mail-type=ALL
source activate /clusterfs/nilah/ayesha/envs/methylseq-repro
export WANDB_API_KEY=$(grep -A2 "machine api.wandb.ai" ~/.netrc | awk '/password/{print $2}')

# Fail fast if GPU is not actually usable — avoids silently running 72h on CPU
python3 -c "import torch; assert torch.cuda.is_available(), 'CUDA unavailable'" || { echo "ERROR: No GPU detected on $(hostname), exiting so SLURM can requeue to a working node."; exit 1; }

CONFIG_DIR=$(methylseqnet-repro-path configs)
UNIQUE_IDENTIFIER="slurm${SLURM_JOB_ID}"
methylseqnet-train --config "$CONFIG_DIR/train/ag-rep0_factorized_atlas+longread_true.gin" --unique-identifier $UNIQUE_IDENTIFIER --start-from-checkpoint slurm35220010 --wandb-run-id slurm35220010 --batch-size 1
