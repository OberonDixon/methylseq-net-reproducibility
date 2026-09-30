#!/bin/bash
#SBATCH --job-name=predict_methylseqnet_longread
#SBATCH --account=fc_streets
#SBATCH --partition=savio3_gpu
## SBATCH --qos=savio_lowprio
#SBATCH --qos=a40_gpu3_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:A40:1
#SBATCH --requeue
#SBATCH --time=6:00:00
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/predict_longread_methylseqnet_%A_%a.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/predict_longread_methylseqnet_%A_%a.err
#SBATCH --exclude=n0214.savio3
## n0264.savio3,n0005.savio3,n0217.savio3,n0215.savio3,n0130.savio4,n0132.savio4,n0134.savio3,n0135.savio3,n0136.savio3,n0137.savio3,n0138.savio3,n0143.savio3,n0144.savio3,n0145.savio3,n0158.savio3,n0159.savio3,n0160.savio3,n0161.savio3,n0174.savio3,n0175.savio3,n0176.savio3
#SBATCH --array=1-3 # Specify the range of array jobs (e.g., 0-2 for 3 configs)
# Command(s) to run:
# Define an array of config files

source activate methylseqnet-ag

PREPROCESSED=$(methylseqnet-repro-path preprocessed_datasets)

KWARGS_ARRAY=(
    "--model-identifier slurm32260895task0 --checkpoint-type temp --dataset-files ${PREPROCESSED}/gm12878_haplotyped/fold4.h5"
    "--model-identifier slurm32260895task1 --checkpoint-type temp --dataset-files ${PREPROCESSED}/gm12878_haplotyped/fold4.h5"
    "--model-identifier slurm37760379task4 --checkpoint-type temp --dataset-files ${PREPROCESSED}/gm12878_haplotyped/fold4.h5"
    "--model-identifier slurm37760379task5 --checkpoint-type temp --dataset-files ${PREPROCESSED}/gm12878_haplotyped/fold4.h5"
)

# Get the config file for this array task
KWARGS=${KWARGS_ARRAY[$SLURM_ARRAY_TASK_ID]}

methylseqnet-predict $KWARGS \
    --model-source local \
    --dataset-keys longread
