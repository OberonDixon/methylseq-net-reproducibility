## Setup
Download scripts will download data to the folders specified in `../../methylseqnet_repro/configs/paths.toml`; this file must be created and populated before running anything. The `methylseqnet-reproducibility` conda environment should be created and activated, `conda activate methylseqnet-reproducibility`.

## Run scripts
From within this directory you can run the following:
```
./download_atac.sh
python download_cage.py
./download_methylation.sh
./download_longread.sh
```