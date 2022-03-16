Main analysis scripts
=====================

A large portion of the study focused on ratio optimization over realistic PD parameters and 1-cmpt PK. The simulations were conducted using parallelization on a cluster running [`code/simulation-batch.R`](https://github.com/bryanmayer/pkpd-bnab-project/blob/main/code/simulation-batch.R). The results of that simulation are stored in a relatively large datafile stored on the Dataverse. [See here](https://github.com/bryanmayer/pkpd-bnab-project/tree/main/output/sim_results) for more details. That datafile can be accessed, processed, and saved into this repository by running `process-sims.Rmd`. Code to analyze that data is located [here](https://github.com/bryanmayer/pkpd-bnab-project/tree/main/code/python). 


The order of the remaining analysis follows the workflowr build order and is self-contained (reproducible) without running `process-sims.Rmd`.

```
wflow_publish(
  c(
    "analysis/_site.yml",
    "analysis/index.Rmd",
    "analysis/opt-sim-background.Rmd",
    "analysis/sensitivity-analysis.Rmd",
    "analysis/bispecific.Rmd",
    "analysis/empirical-case-study.Rmd",
    "analysis/hill-slope-meta.Rmd",
    "analysis/titer-analysis.Rmd",
  )
)
```

The above analysis is also self-contained without re-running `process-catnap.Rmd`. The CATNAP data used in the analysis--stored in `output/`--was pulled on the date displayed in `process-catnap.Rmd`. 

Additional documentation:

 - `opt-sim-background.Rmd` is high level documentation of the R scripts used throughout the project.
 - `sensitivity-analysis.Rmd` is high level documentation of the python scripts used for the sensitivity analysis.
