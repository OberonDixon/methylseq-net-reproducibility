# Pre-packaged reference data
This data is included in the reproducibility repository to simplify re-analysis workflows.
## Data sources
`borzoi-human-splits` contains hg38 regions separated into 8 data folds from https://github.com/calico/borzoi/tree/main/data

`pwms` contains TF position-weight matrices for hundreds of known TFs, created using https://github.com/calico/scBasset/blob/main/examples/ISM/make_pwm.R, with additional manually-created files for CpG-only motifs (CG1,32,3,5,10) and random motifs (UNIFORMRANDOM5,10,15,20)

`imprinted_genes.csv` contains human imprinted genes curated at https://www.geneimprint.com/site/genes-by-species

`transcription_factor_metadata.csv` contains TF methylation sensitivity information synthesized from https://pubmed.ncbi.nlm.nih.gov/28473536/ with additional curated motifs as described in Methods

`41586_2022_5580_MOESM4_ESM.xlsx` contains methylation atlas measurement metadata from Loyfer, N., Magenheim, J., Peretz, A. et al. A DNA methylation atlas of normal human cell types. Nature 613, 355–364 (2023). https://doi.org/10.1038/s41586-022-05580-6

`NIHMS828671-supplement-2.xlsx` contains chromatin compartment maps across cell types from Schmitt AD, Hu M, Jung I, Xu Z, Qiu Y, Tan CL, Li Y, Lin S, Lin Y, Barr CL, Ren B. A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome. Cell Rep. 2016 Nov 15;17(8):2042-2059. doi: 10.1016/j.celrep.2016.10.061. PMID: 27851967; PMCID: PMC5478386

`tabula_sapiens_2.0_table4_TFs-by-celltype.xlsx` contains TF expression levels across dozens of cell types from Stephen R Quake, The Tabula Sapiens Consortium. Tabula Sapiens reveals transcription factor expression, senescence effects, and sex-specific features in cell types from 28 human organs and tissues. bioRxiv 2024.12.03.626516; doi: https://doi.org/10.1101/2024.12.03.626516
