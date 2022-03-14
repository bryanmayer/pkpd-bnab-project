Simulation results output
=====================

This directory exists to store the simulation data as needed. These are relatively big files to store and access locally (100-500MB) and cumbersome to version control.

- [process-sims.Rmd](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/analysis/process-sims.Rmd) pulls the raw simulation data from Dataverse (500 mb file)
- If process-sims.Rmd is executed, three additional datasets are stored in here. A long dataset containing pkpd output optimum and the corresponding optimal ratios. That dataset is processed into two wide datasets for correlation analysis: one for the optimum and one for the ratios.

Raw simulation data citation:

Mayer, Bryan, 2022, "Replication Data for: Optimizing clinical dosing of combination broadly neutralizing antibodies for HIV prevention", https://doi.org/10.7910/DVN/I6PX7D, Harvard Dataverse, V1

R API call: 

```
library(readr)
sim_results = read_csv(url('https://dataverse.harvard.edu/api/access/datafile/6092468'))
```