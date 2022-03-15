Output
==============

Processed output is saved here:

- `processed_catnap_assay.csv` or `processed_catnap_assay.rds`: post-processed CATNAP data. [See processing code](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/analysis/process_catnap.Rmd).
- `sim_results/` - results from the sensitivity analysis. 
  - [See code to run simulations.](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/code/simulation-batch.R).
  - See README in `sim_results/` for data storage and processing notes.
- `empirical-opt/`: results from the [triple bnab combination case study](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/analysis/empirical-case-study.Rmd).
- `bispecific-data.rda`: input data for [bispecific case study](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/analysis/bispecific.Rmd). 
  - Data was subset from [theoretical simulation study](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/analysis/process_sims.Rmd).
