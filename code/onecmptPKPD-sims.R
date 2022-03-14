# This must be loaded with neutralization-funs.R

#' Main function: summarize PKPD of two bnabs with given PK and either given PD or simulated
#'
#' @param pk_parms pk_input = tibble(c0_1, c0_2, hl_1, hl_2, optional: V_1, V_2)
#' @param pd_parms pd_input = tibble(omega_1, muIC50_1, sigIC50_1, omega_2, muIC50_2, sigIC50_2)
#' @param final_time final time of simulation
#' @param virus_seed simulation seed for reproducibility
#' @param total_virus_sims total sims
#' @param IC50s this will skip virus sampling if raw data is provided
#' @param time_res time grid
#' @param return_PK flag to return PK once it is generated
#' @param return_pkpd flag to return PKPD over time data
#'
#' @return tibble of summarized PKPD using both auc and trough wide by interaction and pkpd outcomes
#'
summarize_2mab1cmpt_PKPD = function(pk_parms,
                                    pd_parms = tibble(),
                                    final_time,
                                    virus_seed = as.numeric(Sys.time()),
                                    total_virus_sims = 1000,
                                    IC50s = tibble(),
                                    time_res = 0.1,
                                    return_PK = F,
                                    return_pkpd = F){

  assertthat::assert_that(
    all(map_lgl(list(pk_parms, pd_parms, IC50s), is.data.frame))
    , msg = "all input sets must be data.frame or tibbles")

  # PD setup
  if(nrow(IC50s) == 0){
    virus_sims = .sample_IC50s(total_virus_sims = total_virus_sims, pd_parms = pd_parms, virus_seed = virus_seed)
  } else{
    assertthat::assert_that(all(c("IC50_1", "IC50_2") %in% names(IC50s))
                            ,msg = "IC50s should be wide with IC50_1 and IC50_2")

    virus_sims = IC50s
  }

  # simulate PK data
  pk_data = .pk_1cmpt(pk_parms, final_time, time_res = time_res)
  if(return_PK) return(pk_data)

  pk_summary = pk_data %>% summarize_at(c("conc_1", "conc_2"),
                                        list(
                                          auc = ~ pracma::trapz(time, .x) / max(time),
                                          trough = ~(.x[trough])
                                        ))

  pkpd_dat = pk_data %>%
    group_by(time, trough) %>%
    nest() %>%
    mutate(
      ID50_dat = map(data, ~ (
        tibble(
          ID50_1 = .x$conc_1 / virus_sims$IC50_1,
          ID50_2 = .x$conc_2 / virus_sims$IC50_2
        )
      )),
      neut_res = map(ID50_dat, .generate_interactions),
      neut_summary = map(neut_res, summarize_all, list(mean = mean, median = median))
    )

  if (return_pkpd)
    return(pkpd_dat %>%
             select(-data) %>%
             mutate(virus_seed = virus_seed) %>%
             add_column(!!! pk_summary)
           )

  pkpd_dat %>%
    select(-data, -ID50_dat, -neut_res) %>%
    unnest(neut_summary) %>%
    ungroup() %>%
    summarize_at(vars(-time, -trough), list(
      auc = ~ pracma::trapz(time, .x) / max(time),
      trough = ~(.x[trough])
    )) %>%
    mutate(
      virus_seed = virus_seed,
      total_virus_sims = total_virus_sims
    ) %>%
    add_column(!!! pk_summary)

}

#' Backend function to simulate viral population
#'
#' @param total_virus_sims # total viral population
#' @param pd_parms tibble(omega_1, muIC50_1, sigIC50_1, omega_2, muIC50_2, sigIC50_2)
#' @param virus_seed sampling seed for reproducibility
#'
#' @return tibble of IC50s wide by bnab
#'
.sample_IC50s = function(total_virus_sims, pd_parms, virus_seed = as.numeric(Sys.time())){
  set.seed(virus_seed)
  with(pd_parms, tibble(
    IC50_1 = (1/rbinom(total_virus_sims, size = 1, prob = omega_1)) *
      10^(rnorm(total_virus_sims, muIC50_1, sigIC50_1)),
    IC50_2 = (1/rbinom(total_virus_sims, size = 1, prob = omega_2)) *
      10^(rnorm(total_virus_sims, muIC50_2, sigIC50_2))
  ))

}

#' Backend function to get 1-cmpt PK
#' hl = log(2)/k
#'
#' @param pk_parms pk_input = tibble(c0_1, c0_2, hl_1, hl_2, optional: V_1, V_2)
#' @param final_time end of PK curve
#' @param V Volume of distribution, default 3
#' @param time_res granularity of time grid
#' @param start_time #time grid start
#'
#' @return tibble of concentrations over time wide by bnab
#'
.pk_1cmpt = function(pk_parms, final_time, V = 3, time_res = 0.1, start_time = 0){
  assertthat::assert_that(
    all(flatten_chr(map(c("c0", "hl"), ~paste(., 1:2, sep ="_"))) %in%
          names(pk_parms))
    , msg = "pk_parms missing required variables")

  # PK simulations
  if(any(!c("V_1", "V_2") %in% names(pk_parms))){
    pk_parms$V_1 = V
    pk_parms$V_2 = V
  }

  pk_data = pk_parms  %>%
    gather(parm, value) %>%
    separate(parm, into = c("parm", "mab"), sep = "_", fill = "right") %>%
    filter(!is.na(mab)) %>%
    pivot_wider(names_from = parm, values_from =  value) %>%
    group_by(mab) %>%
    nest() %>%
    mutate(
      pk = map(data, ~tibble(
        time = seq(start_time, final_time, by = time_res),
        conc = .x$c0/.x$V * exp(- time* log(2)/.x$hl)
      ))
    ) %>%
    select(-data) %>%
    unnest(pk)  %>%
    pivot_wider(names_from = mab, values_from = conc, names_prefix = "conc_") %>%
    mutate(trough = time == max(time))

  if(nrow(subset(pk_data, trough)) != 1) stop("issue with trough time assignment")
  if(any(is.na(pk_data))) stop("pk transpose issue")
  pk_data

}

#' Backend helper function to implementing interaction models and PKPD outcomes
#' This function calculates all of the combination PKPD results
#'
#' @param dat tibble(ID50_1, ID50_2), created by summarize function above
#'
#' @return wide titer in -> every outcome by interaction model
#'
.generate_interactions = function(dat) {
  dat %>%
    mutate(
      min_ID50 = pmin(ID50_1, ID50_2),
      minNeut = pmin(titer2neut(ID50_1), titer2neut(ID50_2)),
      max_ID50 = pmax(ID50_1, ID50_2),
      maxNeut = pmax(titer2neut(ID50_1), titer2neut(ID50_2)),
      add_ID50 = ID50_1 + ID50_2,
      additivity = 1 - 1 / (1 + (ID50_1 + ID50_2)),
      BH_ID50 = calc_BHID50(ID50_1, ID50_2),
      BH = 1 - (1 - titer2neut(ID50_1)) * (1 - titer2neut(ID50_2))
    ) %>%
    select(-ID50_1, -ID50_2) %>%
    mutate_at(vars(contains("ID50")), list(
      log10 = function(x)
        pmax(log10(x), -5)
    )) %>%
    rename_at( vars( contains( "_log10") ), list( ~ str_replace(., "_log10", "log10")) ) %>%
    mutate_at(vars(BH, additivity, maxNeut, minNeut), list(neut = identity, iip = neut2iip)) %>%
    select(-BH, -additivity, -maxNeut, -minNeut) %>%
    mutate_at(c("min_ID50", "max_ID50", "add_ID50", "BH_ID50"), list(
      coverage1 = ~(.x > 1),
      coverage10 = ~(.x > 10),
      coverage100 = ~(.x > 100),
      coverage1000 = ~(.x > 1e3),
      coverage1e4 = ~(.x > 1e4)
    )) %>%
    mutate_at(vars(contains("iip")), list(
      coverage1 = ~(.x > 1),
      coverage2 = ~(.x > 2),
      coverage3 = ~(.x > 3),
      coverage4 = ~(.x > 4)
    )) %>%
    mutate_at(vars(contains("coverage")), as.numeric)
}

#' Wrapper to perform ratio optimization by looping over summarize_2mab1cmpt_PKPD
#'
#' @param pk_parms_input pk_input = tibble(c0, hl_1, hl_2, optional: V_1, V_2)
#' @param pd_parms pd_input = tibble(omega_1, muIC50_1, sigIC50_1, omega_2, muIC50_2, sigIC50_2)
#' @param virus_seed simulation seed for reproducibility
#' @param total_virus_sims total sims
#' @param trough_time final time for simulations
#' @param ratio_grid grid for ratio optimization
#'
#' @return stacked tibbles returned from summarize_2mab1cmpt_PKPD, long over ratio
#'
ratio_grid_sims = function(pk_parms_input, pd_parms, virus_seed, total_virus_sims,
                           trough_time, ratio_grid = seq(0, 1, by = 0.01)){


  map_df(ratio_grid, function(ratio){
    pk_parms = pk_parms_input %>%
      mutate(
        ratio = ratio,
        c0_1 = ratio * c0, c0_2 = (1 - ratio) * c0
      )

    summarize_2mab1cmpt_PKPD(
      pk_parms = pk_parms,
      pd_parms = pd_parms,
      final_time = trough_time,
      total_virus_sims = total_virus_sims,
      virus_seed = virus_seed,
      time_res = 1
    ) %>%
      mutate(ratio = ratio)

  })

}
