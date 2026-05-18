# MethylSeqNet Reproducibility
A collection of scripts and notebooks to re-generate results and figures for MethylSeqNet

## Installation

Clone the repository and navigate into the top-level directory containing `environment.yml` and `pyproject.toml`. 

```
git clone https://github.com/OberonDixon/methylseq-net-reproducibility
cd methylseq-net-reproducibility
```

Create a conda environment by running this following command. This will install the `methylseqnet_repro` python package in editable mode so any code changes will be reflected whenever the python kernel is re-started and the `methylseqnet` core package from the latest git commit from https://github.com/OberonDixon/methylseq-net. Conda dependencies necessary for preprocessing and downstream analysis are also installed.

```
conda env create -f environment.yml
```

You can also update your existing environment:

```
conda env update -f environment.yml
```

And you can update while removing any unnecessary dependencies:

```
conda env update -f environment.yml --prune
```

## Pipeline

In order to reproduce all results end-to-end, the following pieces must be run in order. You can skip the Train steps if you wish to use pre-trained checkpoints provided with the manuscript. Preprocessed and prediction output `.h5` files can be made available upon request but are not currently downloadable by the public.

### Data

1. Run each `download` script to download all the input data for analysis.
2. Run the `prepare` scripts to create locally-processed files necessary for preprocessing and predict steps

### Preprocess

1. Generate cell type (short-read) preprocessed `.h5` files using `preprocess-cell-types.sh`
2. Generate haplotyped (long-read) preprocessed `.h5` files using `preprocess-haplotypes.sh`
3. Generate motif insertion preprocessed `.h5` files using `preprocess-motif-inserted.gin`

### Train

1. Train all models using `train-all.sh`. Reference `train-all-continue.sh` for code to continue training from a pre-existing checkpoint, as may be required in HPC environments with maximum wall time or for cases where training is interrupted. 

### Predict

1. Generate prediction `.h5` output files for the cell types (short-read) test set using `predict-cell-types.sh`
2. Generate prediction `.h5` output files for the haplotypes (long-read) test set using `predict-haplotypes.sh`
3. Generate prediction `.h5` output files for the motif inserted synthetic dataset using `predict-motif-inserted.sh`

### Figures
#### 1_overview
 - `io_figure.ipynb` creates DNA sequence, CpG methylation, and output track visualizations used in Figure 1a
 - `model_performances_from_h5.ipynb` generates Pearson correlation bar plots used in Figure 1 and 3 and in Supplementary Figures 1, 2, and 9
 - `embedding_factorization.ipynb` generates bar plots and scatter plots of conditional and unconditional representation composition used in Figure 1 and Supplementary Figure 2
 - `model_performances_manual.ipynb` generates bar plots for training time and memory, based on wandb logs, used in Supplementary Figure 2
 - `methylation_predictions.ipynb` generates bar plots and scatter plots for methylation prediction performance used in Supplementary Figure 2

#### 2_cell_types
 - `cell_type_stats_optimized.ipynb` generates short read cell-type-wise performance plots of various kinds used in Figure 2 and Supplementary Figure 3
#### 3_locus_plots
 - `plot_loci.ipynb` generates multi-tracked single-locus plots with attributions used in Figure 2, 3, and 4 and in Supplementary Figures 4, 5, and 6
#### 4_allelic_imbalance
 - `visualize_imbalance.ipynb` generates allelic imbalance scatterplots used in Figure 3
#### 5_udn_patient
 - `chrX_v_chr13.ipynb` generates chromosome arm performance plots used in Figure 4
#### 6_motif_insertion
 - `plot_motif_insertion.ipynb` generates numerous motif insertion analysis plots used in Figure 5 and Supplementary Figures 7 and 8
 - `explore_cpg_landscape.ipynb` generates histograms of methylation fraction and CpG density across peaks, used in Supplementary Figure 7