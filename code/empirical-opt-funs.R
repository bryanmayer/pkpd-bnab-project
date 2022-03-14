
# ----- These are PD functions are highly specific to this project --------
#' Calculate neutralization from drc model and titer
#' @param titer titer input
#' @param mod5pl drc model object
#'
#' @return scalar
empirical_pd_fun = function(titer, mod5pl) {
  unname(suppressWarnings({drc::PR(mod5pl, titer)}))
}

#' Calculate BH neutralization from drc model and two titers
#' @param titer1 titer input, bnab 1
#' @param titer2 titer input, bnab 2
#' @param mod5pl drc model object
#'
#' @return scalar
emp_prot_combo2 = function(titer1, titer2, mod5pl) {
  1 - (1 - empirical_pd_fun(titer1, mod5pl = mod5pl)) * (1 - empirical_pd_fun(titer2, mod5pl = mod5pl))
}

#' Calculate BH neutralization from drc model and three titers
#' @param titer1 titer input, bnab 1
#' @param titer2 titer input, bnab 2
#' @param titer3 titer input, bnab 3
#' @param mod5pl drc model object
#'
#' @return scalar
emp_prot_combo3 = function(titer1, titer2, titer3, mod5pl) {
  1 - (1 - empirical_pd_fun(titer1, mod5pl = mod5pl)) *
    (1 - empirical_pd_fun(titer2, mod5pl = mod5pl)) *
    (1 - empirical_pd_fun(titer3, mod5pl = mod5pl))
}


#' Error catching function checking numerical issues
check_dose_calc = function(dat, total_dose = 600){
  if(abs(sum(dat$dose) - total_dose) > 1e-3) browser()
  if(any(dat$dose < 0)) browser()

  dat
}

#' Convert 2-cmpt PK parms to bi-exponential parms
#'
#' @param CL
#' @param Vc
#' @param Q
#' @param Vp
#'
#' @return data.frame(A, B, alpha_parm, beta_parm)
convert_2cmp_pkparms = function (CL, Vc, Q, Vp){
  k_sum = Q / Vc + Q / Vp + CL / Vc
  beta_parm = 0.5 * (k_sum - sqrt(k_sum ^ 2 - 4 * Q * CL / (Vp *
                                                              Vc)))
  alpha_parm = Q * CL / (beta_parm * Vp * Vc)
  A = 1 / Vc * (alpha_parm - Q / Vp) / (alpha_parm - beta_parm)
  B = 1 / Vc * (beta_parm - Q / Vp) / (beta_parm - alpha_parm)
  data.frame(A, B, alpha_parm, beta_parm)
}

# ----------- PK functions ------------
#' 2-cmpt absorption PK model
#'
#' @param mtime PK time
#' @param parms data.frame(A, B, alpha_parm, beta_parm)
#' @param ka absorption rate
#' @param dose input dose
#' @param bioavail biovailability
#' @param llod lower limit
#'
#' @return vector of concentration matching input times
absorp_2comp = function (mtime, parms, ka, dose, bioavail = 1, llod = 0) {

  out_sort_index = rank(mtime)
  mtime = sort(mtime)
  model_parms = parms

  A = model_parms$A
  B = model_parms$B
  alpha_parm = model_parms$alpha_parm
  beta_parm = model_parms$beta_parm

  if (bioavail > 1)
    warning("Bioavailability fraction > 1.")
  pmax(.absorption_single_2comp(mtime, A = A, B = B, alpha_parm = alpha_parm,
                                beta_parm = beta_parm, ka = ka, bioavail = bioavail,
                                dose = dose), llod)[out_sort_index]
}


#' Backend function for absorp_2comp
.absorption_single_2comp = function (mtime, A, B, alpha_parm, beta_parm, ka, bioavail, dose) {
  Aabs = A * (ka/(ka - alpha_parm))
  Babs = B * (ka/(ka - beta_parm))
  dose * bioavail * (Aabs * exp(-alpha_parm * mtime) + Babs *
                       exp(-beta_parm * mtime) - (Aabs + Babs) * exp(-ka * mtime))
}


#' 2-cmpt PK over multiple analyte PK parms (stacks absorp_2comp calls)
#'
#' @param pk_parms tibble(analyte, A, B, alpha_parm, beta_parm, ka, F1) each row is different mab (analyte)
#' @param time_grid vector of time
#'
#' @return tibble of PK concentration curves, long by analyte
create_pk = function(pk_parms, time_grid){
  pk_parms %>%
    rowwise() %>%
    mutate(pk_sim = list(tibble(
      time = time_grid,
      conc = absorp_2comp(
        time_grid,
        dose = dose,
        ka = ka,
        bioavail = F1,
        parms = tibble(
          A = A,
          B = B,
          alpha_parm = alpha_parm,
          beta_parm = beta_parm
        )
      )
    ))) %>%
    unnest(cols = pk_sim) %>%
    dplyr::select(analyte, time, conc)
}

# ----------- PD functions ------------
#
# same for pd_data, wide IC50s for each of these mabs

#' Summarized triple combination PKPD under empirical PD assumptions of the case study
#'
#' @param nested_pk very specific input, nested_pk$data = nested wide pk data with specific columns
#' @param pd_data wide IC50 data. columns names are 3BNC117, 10-1074, VRC07-523-LS
#' @param mod5pl drc model object
#'
#' @return mean outcomes across viruses (IC50s) over time
triple_pd_calc = function(nested_pk, pd_data, mod5pl) {
  nested_pk %>%
    mutate(
      protection_res = map(data, ~ (tibble(
        protection = emp_prot_combo3(
          .x$mab_3bnc / pd_data$`3BNC117`,
          .x$mab_1074 / pd_data$`10-1074`,
          .x$mab_vrc07 / pd_data$`VRC07-523-LS`,
          mod5pl = mod5pl
        )
      ))),
      neut_summary = map(
        protection_res,
        summarize_all,
        list(
          mean_iip = ~ (mean(-log10(1 - .x))),
          coverage_50 = ~ mean(. > 0.5),
          coverage_95  = ~ mean(. > .95)
        )
      )
    )
}

#---------- opt functions ----------

#' this is a logistic partitioner
#'
#' @param par_input vector of numbers, n for n+1 analytes
#'
#'  for n ratios, we need n - 1 partitions: this takes an n - 1 vector of reals and returns n ratios summing to 1
#'  first mab = expit(par_input[1])
#'  2nd mab = expit(par_input[1] + exp(par_input[2])) - expit(par_input[1])
#'  kth mab (recursive) = expit(par_input[1] + sum i=2:k exp(par_input[i])) - k-1 mab ratio (recursive)
#'  nth mab =  1 - n-1 mab ratio
#'
#' the first value can be negative, but the rest must be positive, hence exp(par_input[-1])
#' @return n reals summing to 1
#'
calc_opt_ratios = function(par_input){
  trans_par = exp(par_input[-1])
  input = c(-Inf, cumsum(c(par_input[1], trans_par)), Inf)
  out = diff(binomial()$linkinv(input))
  stopifnot(abs(sum(out) - 1) < 1e-5)
  out
}

#' This function organizes 2-bnab optimizations within the 3-bnab framework
#'
#' @param r ratios returned from calc_opt_ratios(scalar)
#' @param remove_mab which mab not in combination
#'
#' @return
#' @export
#'
#' @examples
organize_opt_ratios = function(r, remove_mab){
  if(remove_mab == 1) return(c(0, r[1], r[2]))
  if(remove_mab == 2) return(c(r[1], 0, r[2]))
  if(remove_mab == 3) return(c(r[1], r[2], 0))

  stop("something wrong in organize_opt_ratios")

}

# 3-bnab (3BNC117-T, j = 10-1074-T, k = VRC07-523LS) optimization function
#
#
# see triple_pd_calc
#' Title
#'
#' @param par numeric vector being optimized, length 1 (2-bnab) or 2 (3-bnabs)
#' @param pk_parms see create_pk() above, always analyte levels to include "3BNC117-T", "10-1074-T" "VRC07-523LS
#' @param total_dose total dose
#' @param pk_window passed as time_grid in create_pk() above
#' @param pd_data wide IC50 data. columns names are 3BNC117, 10-1074, VRC07-523-LS
#' @param outcome_var PKPD output for optimization (can only do one at a time), see below
#' @param pk_endpoint PKPD time target for optimization (can only do one at a time), see below
#' @param mod5pl 5PL model for empirical protection estimates (pegu or AMP in empitical-case-study.Rmd)
#' @param remove_mab which mab excluded in 2-bnab combination when length(ratios) == 2
#'
#'  i = 3BNC117-T, j = 10-1074-T, k = VRC07-523LS
#'  if single ratio, remove_mab removes mab at given level: 1 = i, 2 = j, 3 = k (default)
#'  opt_outcomes = crossing(outcome_var = c("coverage_50", "coverage_95"), pk_endpoint = c("trough", "auc"))
#'
#' @return
#'
combo_opt_fun = function(par, pk_parms, total_dose, pk_window, pd_data,
                         outcome_var, pk_endpoint, mod5pl,
                         remove_mab = 3){

  # this expects 2 or 3 mab combinations
  stopifnot(length(par) == 1 | length(par) == 2)
  ratios = calc_opt_ratios(par)
  if(length(ratios) == 2){
    stopifnot(remove_mab %in% 1:3)
    ratios = organize_opt_ratios(ratios, remove_mab)
  }
  if(any(ratios < 0)) browser()

  i = ratios[1]
  j = ratios[2]
  k = ratios[3]
  stopifnot(abs(sum(c(i, j, k)) - 1) < 1e-5)

  pkpd_res = pk_parms %>%
    mutate(dose = if_else(
      analyte == "3BNC117-T",
      total_dose * i,
      if_else(analyte == "10-1074-T", total_dose * j,
              total_dose * k)
    )) %>%
    check_dose_calc(total_dose) %>%
    create_pk(pk_window) %>%
    pivot_wider(values_from = conc, names_from = analyte) %>%
    rename(mab_3bnc = `3BNC117-T`, mab_1074 = `10-1074-T`, mab_vrc07 = `VRC07-523LS`) %>%
    group_by(time) %>%
    nest() %>%
    triple_pd_calc(pd_data, mod5pl = mod5pl) %>%
    select(-data, -protection_res) %>%
    unnest(neut_summary) %>%
    ungroup() %>%
    gather(outcome, value, mean_iip, coverage_50, coverage_95) %>%
    group_by(outcome) %>%
    arrange(time) %>%
    summarize(
      trough_time = max(time),
      trough = value[which.max(time)],
      auc = pracma::trapz(time, value)/trough_time
    )

  #browser()
  -dplyr::filter(pkpd_res, outcome == outcome_var)[[pk_endpoint]]

}
