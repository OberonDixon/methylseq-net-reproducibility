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

### Data

### Preprocessing

### Training

### Figures