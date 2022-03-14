args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")

}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  print(args)
}

increment = 307
sim_set = c(start + 1, start + increment)
print(sim_set)


# this optimizes ratio by outcome using grid search search, and more systemtic parameter search
# because this generates all of the sims, everything can be done at the same time and the opt ratio can be selected after

library(tidyverse)
library(doParallel)

registerDoParallel(8)
getDoParWorkers()

source("onecmptPKPD-sims.R")
source("neutralization-funs.R")


## ------------------------------------------------------------------------

# dose (c0): 1; ratio: 1; half-lives: 2; omega: 2, mu: 2, sigma: 2; omega_cor: 1; ic50_cor: 1
# 10 or 12 (if using correlations)

doses = c(150, 300, 600, 1200, 2400)

# parameter grid

hl = c(7, 28, 42, 84)
omega = c(0.33, 0.67, 1)
muIC50 = c(-3, -2, -1)
sigIC50 = c(0.25, 0.5, 1)


full_combinations = crossing(
  hl_1 = hl,
  hl_2 = hl,
  omega_1 = omega,
  omega_2 = omega,
  muIC50_1 = muIC50,
  muIC50_2 = muIC50,
  sigIC50_1 = sigIC50,
  sigIC50_2 = sigIC50
) %>%
  filter(hl_1 <= hl_2)

trough_time = 84
total_sims = nrow(full_combinations) #set to total_grid to match lhs
total_virus_sims = 1000

# since the iterations are determined externally and the total dimensions could change
if(sim_set[2] > total_sims) sim_set[2] = total_sims
if(sim_set[2] < sim_set[1]) stop("start simulation index greater than total simulations")


seed = 125
set.seed(seed)
virus_seeds = sample(1e9, total_sims, replace = F)

system.time({
  sim_resx = plyr::ldply(sim_set[1]:sim_set[2], function(i){
    #print(i)
    map_df(doses, function(dose){

      wide_parms = full_combinations[i,] %>%
        mutate(c0 = dose)

      pk_parms_input =  select(wide_parms, c0, contains("hl_"))

      pd_parms = wide_parms %>%
        select(contains("mu"), contains("omega"), contains("sig"))

      ratio_grid_sims(pk_parms_input = pk_parms_input,
                 pd_parms = pd_parms,
                 virus_seed = virus_seeds[i],
                 total_virus_sims = total_virus_sims,
                 trough_time = trough_time,
                 ratio_grid = seq(0, 1, by = .1)) %>%
        add_column(!!! wide_parms) %>%
        select(names(wide_parms), everything()) %>%
        mutate(sim = i)
    })
    }, .parallel = T)
  })

print("the simulations ran!")

grid_str = paste0("grid", sim_set[1], "_", sim_set[2])

out_time = format(Sys.time(), "-%Y-%m-%d-%H.%M-")
print(out_time)
fout = paste0(grid_str, out_time, seed, ".csv")
print(fout)
write_csv(sim_resx, fout)
