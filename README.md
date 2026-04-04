# Cancer-Associated Fibroblast Gene Signatures as Novel Prognostic Indicators in Cervical Cancer

This repository contains R scripts and data used to construct and validate a prognostic model based on cancer-associated fibroblast (CAF)-related genes in cervical cancer patients.

## Data
- TCGA cervical cancer dataset
- GEO datasets: GSE63514, GSE7803

## Methods
- Differential gene expression analysis (DEGs)
- XCell algorithm (https://comphealth.ucsf.edu/app/xcell) for fibroblast enrichment score
- ESTIMATE algorithm for stromal cell enrichment score
- WGCNA analysis
- Validation using Ridge Cox regression model and Time-dependent ROC analysis

## Files Description
- `63514 WGCNA.R` – WGCNA analysis for GSE63514
- `7803 WGCNA.R` – WGCNA analysis for GSE7803
- `ESTIMATES63514.R` – ESTIMATE score calculation for GSE63514
- `ESTIMATES7803.R` – ESTIMATE score calculation for GSE7803
- `GSE63514 DEGs analysis.R` – DEG analysis for GSE63514
- `GSE7803 DEGs analysis.R` – DEG analysis for GSE7803
- `Validate (Prognostic model).R` – Model validation

## How to Run
You can run each steps by using the raw data and other supplement files in each step 

## Requirements
- R version ≥ 4.0
- Packages: survival, glmnet, ggplot2, survivalROC, tidyverse

## Contact
For questions, please contact: Khaohom2802@gmail.com
