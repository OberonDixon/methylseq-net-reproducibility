#!/bin/bash
#SBATCH --job-name=preprocess-longread
#SBATCH --account=fc_nilah
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --output=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A.out
#SBATCH --error=/clusterfs/nilah/oberon/lightning/sbatch_logs/preprocess_methylseqnet_%A.err

module load anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
source activate methylseqnet-reproducibility

CONFIG_DIR=$(methylseqnet-repro-path configs)

methylseqnet-preprocess --config "$CONFIG_DIR/preprocess/preprocess_ag_longread.gin"
methylseqnet-preprocess --config "$CONFIG_DIR/preprocess/preprocess_ag_longread_chrX.gin"
methylseqnet-preprocess --config "$CONFIG_DIR/preprocess/preprocess_ag_longread_haplotyped.gin"