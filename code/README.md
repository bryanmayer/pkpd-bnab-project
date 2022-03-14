R scripts and functions for bnab optimization analysis
====================================

Functions within the scripts are documented. 

For analysis function, [more details here]():
 - `neutralization-funs.R`: functions performing PD calculations
 - `onecmptPKPD-sims.R`: functions for main 1-cmpt theoretical optimization simulations
 - `empirical-opt-funs.R`: functions to perform 3-bnab combination optimizations with 2-cmpt PK model

 `simulation-batch.R` was the code used on the cluster to perform the 1-cmpt simulations utilizing parallelization. Requires neutralization-funs.R and onecmptPKPD-sims.R.
 