#' ---
#' title: Simulation study for asSAMs
#' abstract: For additional smaller simulation concerning estimation and uncertainty quantification, and a Bayesian approach to fitting SAMs is included
#' author: FKCH
#' ---
rm(list = ls())
library(tidyverse)
library(mvtnorm)
library(doParallel)
library(tictoc)
library(cluster)
library(ecomix)
library(foreach)
library(cmdstanr)
library(label.switching)
library(assam)


##----------------------
#' # Simulation parameters
##----------------------
load(file = "simsetting1_truemodel.RData")

num_X <- 9
num_spp <- 50
num_archetype <- 3

new_formula <- ~ GBR_BATHY + GBR_TS_BSTRESS + GA_CRBNT + GA_GRAVEL + GA_MUD + CRS_O2_AV + CRS_S_AV + CRS_T_AV + SW_CHLA_AV

true_spp_effects <- true_spp_effects[1:num_spp,,drop=FALSE]
true_betas <- true_betas[1:num_archetype,]
true_dispparam <- true_dispparam[1:num_spp]
true_powerparam <- true_powerparam[1:num_spp]
true_mixprop <- c(0.25,0.5,0.25)

##----------------------
#' # Set up simulation function
##----------------------
simfn <- function(seed,
                  num_units = 250) {
    
    num_boot <- 100
    
    ##----------------------
    # Generate some multivariate abundance (count) data from a SAM
    ##----------------------
    set.seed(052025)
    sel_sites <- sample(1:nrow(environ_dat), num_units, replace = ifelse(num_units > nrow(environ_dat), TRUE, FALSE))
    covariate_dat <- environ_dat[sel_sites,, drop = FALSE]
    set.seed(NULL)
    
    simdat <- create_samlife(family = binomial(),
                             formula = new_formula, 
                             data = covariate_dat,
                             betas = true_betas,
                             spp_effects = true_spp_effects,
                             spp_dispparam = true_dispparam,
                             spp_powerparam = true_powerparam,
                             mixing_proportion = true_mixprop,
                             seed = seed)
    
    tictoc::tic.clearlog()
    
    
    ##----------------------
    #' # Method 1: asSAMs
    ##----------------------
    #' ## Construct the initial stacked species model fits that will be used throughout the model selection process below
    tictoc::tic("assams_prefit")
    samfit_prefit <- assam(y = simdat$y,
                           formula = new_formula,
                           data = covariate_dat,
                           family = binomial(),
                           num_archetypes = 2, #' This is arbitrary and does not matter
                           num_cores = detectCores() - 4,
                           do_assam_fit = FALSE)
    tictoc::toc(log = TRUE)
    

    tictoc::tic("assams_final_pointestimate_only")
    samfit_final_pointest <- try(assam(y = simdat$y,
                                       formula = new_formula,
                                       data = covariate_dat,
                                       family = binomial(),
                                       num_archetypes = num_archetype,
                                       uncertainty_quantification = FALSE,
                                       supply_quadapprox = samfit_prefit,
                                       num_cores = detectCores() - 4),
                                 silent = TRUE)
    samfit_final_pointest$sdmTMB_fits <- samfit_final_pointest$linear_predictor <- NULL
    tictoc::toc(log = TRUE)

    
    tictoc::tic("assams_final_fastbootstrap")
    samfit_final_boot_fast <- try(assam(y = simdat$y,
                                    formula = new_formula,
                                    data = covariate_dat,
                                    family = binomial(),
                                    num_archetypes = num_archetype,
                                    uncertainty_quantification = TRUE,
                                    bootstrap_control = list(method = "fast", num_boot = num_boot),
                                    supply_quadapprox = samfit_prefit,
                                    num_cores = detectCores() - 4),
                                    silent = TRUE)
    samfit_final_boot_fast$sdmTMB_fits <- samfit_final_boot_fast$linear_predictor <- NULL
    samfit_final_boot_fast$bootstrap_posterior_probability <- NULL
    samfit_final_boot_fast$bootstrap_parameters <- NULL
    tictoc::toc(log = TRUE)
        
    assams_tictoclog_txt <- tictoc::tic.log(format = TRUE)
    assams_tictoclog_lst <- tictoc::tic.log(format = FALSE)
    tictoc::tic.clearlog()
    rm(samfit_prefit)
    gc()

    
    ##----------------------
    #' # Method 2: species_mix
    ##----------------------
    sam_form <- paste0('cbind(',paste("spp", 1:num_spp, sep = "", collapse = ','), ") ~ ", as.character(new_formula)[2]) %>% 
        as.formula
    
    tictoc::tic("speciesmix_pointestimateonly")
    speciesmix_pointestimate <- try(species_mix(archetype_formula = sam_form,
                                                species_formula = ~ 1,
                                                data = data.frame(simdat$y, covariate_dat), 
                                                family = "bernoulli", 
                                                nArchetypes = num_archetype),
                                    silent = TRUE)
    tictoc::toc(log = TRUE)
    
    tictoc::tic("speciesmix_addonbootstrap")
    speciesmix_boot <- ecomix::species_mix.bootstrap(speciesmix_pointestimate,
                                                     nboot = num_boot,
                                                     mc.cores = detectCores() - 4,
                                                     quiet = FALSE)
    speciesmix_boot
    speciesmix_pointestimate$vcov <- vcov(speciesmix_pointestimate, speciesmix_boot)
    speciesmix_pointestimate$summary <- summary(speciesmix_pointestimate) %>%
        as.data.frame() %>%
        rownames_to_column(var = "effect") %>%
        filter(str_detect(effect, "beta")) %>%
        mutate(significance = ifelse(p < 0.05, 1, 0))
    tictoc::toc(log = TRUE)
    rm(speciesmix_boot)
    speciesmix_pointestimate$titbits <- NULL
    speciesmix_pointestimate$mus <- NULL
    speciesmix_pointestimate$terms <- NULL

    speciesmix_tictoclog_txt <- tictoc::tic.log(format = TRUE)
    speciesmix_tictoclog_lst <- tictoc::tic.log(format = FALSE)
    
    
    ##----------------------
    #' # Method 3: Fitting a SAM using Bayesian MCMC via stan
    #' Then account for label switching using the Stephens (2000) algorithm as implemented in the label.switching package
    ##----------------------
    mod <- cmdstan_model("sam.stan")
    
    tictoc::tic("stanfit")
    stan_data <- list(
        num_sites = num_units,
        num_species = num_spp,
        num_X = num_X,
        num_archetypes = num_archetype,
        y = simdat$y,
        X = model.matrix(new_formula, data = covariate_dat)[,-1],
        alpha_dirichlet = rep(1, num_archetype) # Uniform prior on mixing
        )
    init_fun <- function() {
        list(beta = true_betas, pi = true_mixprop)
        }
    
    
    fit_stan <- mod$sample(data = stan_data,
                           parallel_chains = 4,
                           init = init_fun,
                           iter_warmup = 1000,
                           iter_sampling = 2000)
    z_draws <- fit_stan$draws(variables = "z_prob") %>%
        posterior::as_draws_matrix(z_draws)
    P_array <- array(NA, dim = c(nrow(z_draws), num_spp, num_archetype))
    for (k in 1:num_archetype) {
        for (j in 1:num_spp) {
            var_name <- sprintf("z_prob[%d,%d]", j, k)
            P_array[, j, k] <- z_draws[, var_name]
            }
        }
    ls_run <- stephens(P_array)
    perms <- ls_run$permutations
    rm(ls_run, z_draws, P_array)
    
    beta_raw <- fit_stan$draws("beta") %>%
        posterior::as_draws_matrix()
    pi_raw <- fit_stan$draws("mixprop") %>%
        posterior::as_draws_matrix()
    beta_fixed <- beta_raw
    pi_fixed <- pi_raw
    for (l0 in 1:nrow(beta_raw)) {
        for (k0 in 1:num_X) {
            # Extract the original values for all archetypes for variable k
            orig_vals <- sapply(1:num_archetype, function(g) { beta_raw[l0, sprintf("beta[%d,%d]", g, k0)] })
            
            # Reassign them based on the permutation
            for (g0 in 1:num_archetype) {
                beta_fixed[l0, sprintf("beta[%d,%d]", g0, k0)] <- orig_vals[perms[l0,g0]]
                }
            }
        
        rm(orig_vals)
        pi_fixed[l0,] <- pi_raw[l0,][perms[l0,]]
        }
    rm(beta_raw, pi_raw, l0, k0, g0)
    fit_stan_summary <- list(beta = posterior::summarize_draws(beta_fixed,
                                                               mean,
                                                               median,
                                                               sd,
                                                               ~quantile(.x, probs = c(0.025, 0.975))),
                             spp_effects = posterior::summarize_draws(fit_stan$draws(variables = "alpha") %>% posterior::as_draws_matrix(z_draws),
                                                                      mean,
                                                                      median,
                                                                      sd,
                                                                      ~quantile(.x, probs = c(0.025, 0.975))),
                             mixture_proportions = posterior::summarize_draws(pi_fixed,
                                                                              mean,
                                                                              median,
                                                                              sd,
                                                                              ~quantile(.x, probs = c(0.025, 0.975))))
    tictoc::toc(log = TRUE)
    rm(fit_stan, stan_data)
    stanfit_tictoclog_txt <- tictoc::tic.log(format = TRUE)
    stanfit_tictoclog_lst <- tictoc::tic.log(format = FALSE)
    tictoc::tic.clearlog()
    gc()
    
    
    
    ##----------------------
    #' # Final return
    ##----------------------
    true_archetype_label <- simdat$archetype_label
    
    save(samfit_final_pointest,
         samfit_final_boot_fast,
         assams_tictoclog_txt,
         assams_tictoclog_lst,
         speciesmix_pointestimate,
         speciesmix_tictoclog_txt,
         speciesmix_tictoclog_lst,
         fit_stan_summary,
         stanfit_tictoclog_txt,
         stanfit_tictoclog_lst,
         true_archetype_label,
         file = paste0("simdat_numunits", num_units, "_dataset", seed, ".RData"))
    
    tictoc::tic.clearlog()
    }



##----------------------
#' # Run the simulation
##----------------------
registerDoParallel(cores = detectCores() - 5)

for(k in 1:200)
    simfn(seed = k, num_units = 250)

for(k in 1:200)
    simfn(seed = k, num_units = 500)



##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 (2024-04-24)
# os       Linux Mint 20.3
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# rstudio  2024.12.0+467 Kousa Dogwood (desktop)
# pandoc   NA
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package         * version    date (UTC) lib source
# abind             1.4-8      2024-09-12 [1] CRAN (R 4.4.0)
# assam           * 0.1        2025-10-06 [1] local
# assertthat        0.2.1      2019-03-21 [1] CRAN (R 4.4.0)
# cli               3.6.4      2025-02-13 [1] CRAN (R 4.4.0)
# cluster         * 2.1.6      2023-12-01 [4] CRAN (R 4.3.2)
# coda              0.19-4.1   2024-01-31 [1] CRAN (R 4.4.0)
# codetools         0.2-19     2023-02-01 [4] CRAN (R 4.2.2)
# collapse          2.0.16     2024-08-21 [1] CRAN (R 4.4.0)
# colorspace        2.1-1      2024-07-26 [1] CRAN (R 4.4.0)
# combinat          0.0-8      2012-10-29 [1] CRAN (R 4.4.0)
# doParallel      * 1.0.17     2022-02-07 [1] CRAN (R 4.4.0)
# dplyr           * 1.1.4      2023-11-17 [1] CRAN (R 4.4.0)
# ecomix          * 1.0.0      2025-04-11 [1] Github (skiptoniam/ecomix@15b5e51)
# emmeans           1.10.4     2024-08-21 [1] CRAN (R 4.4.0)
# estimability      1.5.1      2024-05-12 [1] CRAN (R 4.4.0)
# forcats         * 1.0.0      2023-01-29 [1] CRAN (R 4.4.0)
# foreach         * 1.5.2      2022-02-02 [1] CRAN (R 4.4.0)
# generics          0.1.3      2022-07-05 [1] CRAN (R 4.4.0)
# ggplot2         * 3.5.2      2025-04-09 [1] CRAN (R 4.4.0)
# glue              1.8.0      2024-09-30 [1] CRAN (R 4.4.0)
# gtable            0.3.6      2024-10-25 [1] CRAN (R 4.4.0)
# hms               1.1.3      2023-03-21 [1] CRAN (R 4.4.0)
# iterators       * 1.0.14     2022-02-05 [1] CRAN (R 4.4.0)
# label.switching   1.8        2019-07-01 [1] CRAN (R 4.4.0)
# lattice           0.22-5     2023-10-24 [4] CRAN (R 4.3.1)
# lifecycle         1.0.4      2023-11-07 [1] CRAN (R 4.4.0)
# lpSolve           5.6.20     2023-12-10 [1] CRAN (R 4.4.0)
# lubridate       * 1.9.3      2023-09-27 [1] CRAN (R 4.4.0)
# magrittr          2.0.3      2022-03-30 [1] CRAN (R 4.4.0)
# Matrix            1.6-5      2024-01-11 [4] CRAN (R 4.3.3)
# mgcv              1.9-1      2023-12-21 [4] CRAN (R 4.3.2)
# munsell           0.5.1      2024-04-01 [1] CRAN (R 4.4.0)
# mvtnorm         * 1.3-3      2025-01-10 [1] CRAN (R 4.4.0)
# nlme              3.1-163    2023-08-09 [4] CRAN (R 4.3.1)
# pillar            1.10.2     2025-04-05 [1] CRAN (R 4.4.0)
# pkgconfig         2.0.3      2019-09-22 [1] CRAN (R 4.4.0)
# purrr           * 1.0.2      2023-08-10 [1] CRAN (R 4.4.0)
# quadprog          1.5-8      2019-11-20 [1] CRAN (R 4.4.0)
# R6                2.6.1      2025-02-15 [1] CRAN (R 4.4.0)
# Rcpp              1.0.14     2025-01-12 [1] CRAN (R 4.4.0)
# RcppHungarian     0.3        2023-09-05 [1] CRAN (R 4.4.0)
# readr           * 2.1.5      2024-01-10 [1] CRAN (R 4.4.0)
# RhpcBLASctl       0.23-42    2023-02-11 [1] CRAN (R 4.4.0)
# rlang             1.1.5      2025-01-17 [1] CRAN (R 4.4.0)
# rstudioapi        0.16.0     2024-03-24 [1] CRAN (R 4.4.0)
# scales            1.3.0      2023-11-28 [1] CRAN (R 4.4.0)
# sdmTMB          * 0.6.0.9004 2024-09-16 [1] Github (pbs-assess/sdmTMB@cb83a62)
# sessioninfo       1.2.2      2021-12-06 [1] CRAN (R 4.4.0)
# stringi           1.8.4      2024-05-06 [1] CRAN (R 4.4.0)
# stringr         * 1.5.1      2023-11-14 [1] CRAN (R 4.4.0)
# tdigest           0.4.1      2022-10-04 [1] CRAN (R 4.4.0)
# tibble          * 3.2.1      2023-03-20 [1] CRAN (R 4.4.0)
# tictoc          * 1.2.1      2024-03-18 [1] CRAN (R 4.4.0)
# tidyr           * 1.3.1      2024-01-24 [1] CRAN (R 4.4.0)
# tidyselect        1.2.1      2024-03-11 [1] CRAN (R 4.4.0)
# tidyverse       * 2.0.0      2023-02-22 [1] CRAN (R 4.4.0)
# timechange        0.3.0      2024-01-18 [1] CRAN (R 4.4.0)
# TMB               1.9.11     2024-04-03 [1] CRAN (R 4.4.0)
# tweedie           2.3.5      2022-08-17 [1] CRAN (R 4.4.0)
# tzdb              0.4.0      2023-05-12 [1] CRAN (R 4.4.0)
# vctrs             0.6.5      2023-12-01 [1] CRAN (R 4.4.0)
# withr             3.0.2      2024-10-28 [1] CRAN (R 4.4.0)
# xtable            1.8-4      2019-04-21 [1] CRAN (R 4.4.0)
# 
# [1] /home/fhui28/R/x86_64-pc-linux-gnu-library/4.4
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
