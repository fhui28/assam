#' ---
#' title: Simulation study for asSAMs
#' abstract: For setting 1 concerning estimation and uncertainty quantification -- constructing some profile likelihoods for select coefficients in the stacked SDM and comparing them to the quadratic/Gaussian approximation used in asSAMs.
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
library(sdmTMB)
library(assam)
library(ProfileLikelihood)


##----------------------
#' # Simulation parameters
##----------------------
load(file = "simsetting1_truemodel.RData")

num_X <- 9
num_spp <- length(true_spp_effects)
num_archetype <- nrow(true_betas)

new_formula <- ~ GBR_BATHY + GBR_TS_BSTRESS + GA_CRBNT + GA_GRAVEL + GA_MUD + CRS_O2_AV + CRS_S_AV + CRS_T_AV + SW_CHLA_AV

#' For the additional simulation involving datasets with reduced prevalence, use the following modification which sets all species-specific effects are set to -2.5. 
#' On average, this makes the mean prevalence around 9-10% with a range from 7 to 12% for binary and negative binomial responses.
# true_spp_effects <- matrix(runif(num_spp, -2.5, -2.5), ncol = 1) 


##----------------------
#' # Set up simulation function
##----------------------
simfn <- function(seed,
                  family = binomial(),
                  num_units = 250) {
    
    ##----------------------
    # Generate some multivariate abundance data from a SAM
    ##----------------------
    set.seed(052025)
    sel_sites <- sample(1:nrow(environ_dat), num_units, replace = ifelse(num_units > nrow(environ_dat), TRUE, FALSE))
    covariate_dat <- environ_dat[sel_sites,, drop = FALSE]
    set.seed(NULL)
    
    simdat <- create_samlife(family = family,
                             formula = new_formula, 
                             data = covariate_dat,
                             betas = true_betas,
                             spp_effects = true_spp_effects,
                             spp_dispparam = true_dispparam,
                             spp_powerparam = true_powerparam,
                             mixing_proportion = true_mixprop,
                             seed = seed)
    while(min(table(simdat$archetype_label)) < 5) {
        simdat <- create_samlife(family = family,
                                 formula = new_formula,
                                 data = covariate_dat,
                                 betas = true_betas,
                                 spp_effects = true_spp_effects,
                                 spp_dispparam = true_dispparam,
                                 spp_powerparam = true_powerparam,
                                 mixing_proportion = true_mixprop)
        }
    
    
    
    ##----------------------
    #' # Construct profile likelihood of select coefficients in stacked SDM and compare it to the Gaussian approximation used in asSAMs
    ##----------------------
    allspp_proflik_dat <- vector("list", num_spp)
    
    for(k0 in 1:num_spp) {
        new_full_formula <- update.formula(new_formula, resp ~ .)
        # resp ~ GBR_BATHY + GBR_TS_BSTRESS + GA_CRBNT + GA_GRAVEL + GA_MUD + CRS_O2_AV + CRS_S_AV + CRS_T_AV + SW_CHLA_AV
        if(family$family == "binomial") 
            stackedsdm <- glm(new_full_formula, 
                              data = data.frame(resp = simdat$y[,k0], covariate_dat), 
                              family = family)
        if(family$family == "nbinom2")
            stackedsdm <- glm(new_full_formula, 
                              data = data.frame(resp = simdat$y[,k0], covariate_dat), 
                              family = negative.binomial(theta = true_dispparam[k0], link = "log"))
        # if(family$family == "tweedie")
        #     stackedsdm <- glm(new_full_formula, 
        #                       data = data.frame(resp = simdat$y[,k0], covariate_dat), 
        #                       family = statmod::tweedie(var.power = true_powerparam[k0], link.power = 0))
        
        current_terms <- labels(terms(new_formula))
        
        all_profile_lik_dat <- NULL
        for(k1 in 1:num_X) {
            paste0("Fitting profile likelihood for species ", k0, " and covariate ", current_terms[k1]) %>% message
            
            new_subtract_formula <- reformulate(termlabels = current_terms[-k1], response = "resp")
            if(family$family == "binomial")
                getprofLik <- profilelike.glm(new_subtract_formula, 
                                              data = data.frame(resp = simdat$y[,k0], covariate_dat), 
                                              profile.theta = current_terms[k1], 
                                              family = family, 
                                              length = 500)
            if(family$family == "nbinom2")
                getprofLik <- profilelike.glm(new_subtract_formula, 
                                              data = data.frame(resp = simdat$y[,k0], covariate_dat), 
                                              profile.theta = current_terms[k1], 
                                              family = negative.binomial(theta = true_dispparam[k0], link = "log"), 
                                              length = 500)
            # if(family$family == "tweedie")
            #     getprofLik <- tweedieprofilelike.glm(new_subtract_formula, 
            #                                          data = data.frame(resp = simdat$y[,k0], covariate_dat), 
            #                                          profile.theta = current_terms[k1], 
            #                                          family = statmod::tweedie(var.power = true_powerparam[k0], link.power = 0), 
            #                                          length = 500)
            
            define_quad_fn <- function(x) {
                # Get the MLE and standard error for the coefficient of interest
                mle <- coef(stackedsdm)[current_terms[k1]]
                infoinv <- vcov(stackedsdm)[current_terms[k1], current_terms[k1]]
                
                out <- dnorm(x, mean = mle, sd = sqrt(infoinv), log = FALSE) 
                }
        
            #' Now we have the profile likelihood values for a range of theta values (the coefficient of interest) and we can compare it to the Gaussian approximation. 
            #' Note both likelihoods are normalized to have a maximum of 1 so we can compare them directly.
            proflogL_dat <- data.frame(theta = getprofLik$theta,
                                       assams_quadapprox = sapply(getprofLik$theta, define_quad_fn),
                                       profile_lik = getprofLik$profile.lik.norm) %>% 
                mutate(assams_quadapprox = assams_quadapprox/max(assams_quadapprox)) %>% 
                pivot_longer(-theta, names_to = "method", values_to = "likelihood") %>% 
                mutate(method = recode(method, 
                                       assams_quadapprox = "asSAMs Gaussian approximation", 
                                       profile_lik = "Profile likelihood"))
            all_profile_lik_dat <- bind_rows(all_profile_lik_dat,
                                             data.frame(covariate = current_terms[k1], proflogL_dat))
            
            rm(getprofLik, define_quad_fn, proflogL_dat)
            }
    
        allspp_proflik_dat[[k0]] <- all_profile_lik_dat %>% 
            mutate(species = paste0("Species ", k0))
        }
        
    save(simdat,
         allspp_proflik_dat,
         file = paste0("simdat_", family$family, "_numunits", num_units, "_profilelikelihood_seed", seed, ".RData"))

    }


##----------------------
#' # Run the simulation
##----------------------
registerDoParallel(cores = detectCores() - 5)

n_seq <- c(250, 500, 1000, 2000)
foreach(k = 1:length(n_seq)) %dopar% simfn(seed = 1, 
                                           num_units = n_seq[k], 
                                           family = binomial())

n_seq <- c(250, 500, 1000, 2000)
foreach(k = 1:length(n_seq)) %dopar% simfn(seed = 1, 
                                           num_units = n_seq[k], 
                                           family = nbinom2())

# ggplot(allspp_proflik_dat[[num_spp]] %>%
#            mutate(covariate = fct_inorder(covariate),
#                   covariate = fct_recode(covariate, 
#                                          "Bathymetry" = "GBR_BATHY", 
#                                          "Bottom stress" = "GBR_TS_BSTRESS", 
#                                          "Carbonate" = "GA_CRBNT", 
#                                          "Gravel" = "GA_GRAVEL", 
#                                          "Mud" = "GA_MUD", 
#                                          "Oxygen" = "CRS_O2_AV", 
#                                          "Salinity" = "CRS_S_AV", 
#                                          "Temperature" = "CRS_T_AV", 
#                                          "Chlorophyll-a" = "SW_CHLA_AV")),
#        aes(x = theta, y = likelihood, color = method, group = method)) +
#     geom_line() +
#     facet_wrap(~covariate, scales = "free_x") +
#     scale_color_discrete_qualitative() +
#     labs(x = "Coefficient value", y = "Normalized likelihood", color = "Method",
#          title = "Count responses (prevalence = 0.094)") +
#     theme_bw() +
#     theme(legend.position = "bottom")
# 
# 
# ggsave("setting1_countreponses_n2000_profilelikelihood.pdf", width = 10, height = 10)



##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.3 (2024-02-29)
# os       Linux Mint 22
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2026-04-10
# rstudio  2026.01.0+392 Apple Blossom (desktop)
# pandoc   3.1.3 @ /usr/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package           * version    date (UTC) lib source
# abind               1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# assam             * 0.1        2026-03-05 [1] local
# assertthat          0.2.1      2019-03-21 [1] CRAN (R 4.3.3)
# cli                 3.6.5      2025-04-23 [1] CRAN (R 4.3.3)
# cluster           * 2.1.6      2023-12-01 [4] CRAN (R 4.3.3)
# coda                0.19-4     2020-09-30 [1] CRAN (R 4.3.0)
# codetools           0.2-19     2023-02-01 [4] CRAN (R 4.2.2)
# collapse            2.0.10     2024-02-17 [3] CRAN (R 4.3.2)
# combinat            0.0-8      2012-10-29 [1] CRAN (R 4.3.3)
# dichromat           2.0-0.1    2022-05-02 [3] CRAN (R 4.2.0)
# doParallel        * 1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
# dplyr             * 1.1.2      2023-04-20 [1] CRAN (R 4.3.0)
# ecomix            * 1.0.0      2026-03-06 [1] Github (skiptoniam/ecomix@15b5e51)
# emmeans             1.9.0      2023-12-18 [1] CRAN (R 4.3.1)
# estimability        1.4.1      2022-08-05 [1] CRAN (R 4.3.1)
# farver              2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
# foreach           * 1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
# generics            0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# ggplot2           * 4.0.1      2025-11-14 [1] CRAN (R 4.3.3)
# glue                1.8.0      2024-09-30 [1] CRAN (R 4.3.3)
# gtable              0.3.6      2024-10-25 [1] CRAN (R 4.3.3)
# hms                 1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
# iterators         * 1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
# label.switching     1.8        2019-07-01 [1] CRAN (R 4.3.3)
# lattice             0.22-5     2023-10-24 [4] CRAN (R 4.3.3)
# lifecycle           1.0.5      2026-01-08 [1] CRAN (R 4.3.3)
# lpSolve             5.6.20     2023-12-10 [3] CRAN (R 4.3.2)
# lubridate         * 1.9.2      2023-02-10 [1] CRAN (R 4.3.0)
# magrittr            2.0.4      2025-09-12 [1] CRAN (R 4.3.3)
# MASS                7.3-60.0.1 2024-01-13 [4] CRAN (R 4.3.2)
# Matrix              1.6-0      2023-07-08 [1] CRAN (R 4.3.1)
# mgcv                1.9-1      2023-12-21 [4] CRAN (R 4.3.2)
# multcomp            1.4-25     2023-06-20 [3] CRAN (R 4.3.1)
# mvtnorm           * 1.2-2      2023-06-08 [1] CRAN (R 4.3.0)
# nlme                3.1-164    2023-11-27 [4] CRAN (R 4.3.3)
# pillar              1.11.1     2025-09-17 [1] CRAN (R 4.3.3)
# pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# ProfileLikelihood * 1.3        2023-08-25 [1] CRAN (R 4.3.3)
# purrr             * 1.2.1      2026-01-09 [1] CRAN (R 4.3.3)
# quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.3.0)
# R6                  2.6.1      2025-02-15 [1] CRAN (R 4.3.3)
# RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                1.1.1      2026-01-10 [1] CRAN (R 4.3.3)
# RcppHungarian       0.3        2023-09-05 [1] CRAN (R 4.3.3)
# readr             * 2.1.4      2023-02-10 [1] CRAN (R 4.3.0)
# RhpcBLASctl         0.23-42    2023-02-11 [3] CRAN (R 4.3.2)
# rlang               1.1.7      2026-01-09 [1] CRAN (R 4.3.3)
# rsconnect           1.2.0      2023-12-15 [3] CRAN (R 4.3.2)
# rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.3.3)
# S7                  0.2.1      2025-11-14 [1] CRAN (R 4.3.3)
# sandwich            3.0-2      2022-06-15 [1] CRAN (R 4.3.1)
# scales              1.4.0      2025-04-24 [1] CRAN (R 4.3.3)
# sdmTMB            * 0.8.1      2026-01-08 [1] CRAN (R 4.3.3)
# sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.3.1)
# stringi             1.8.7      2025-03-27 [1] CRAN (R 4.3.3)
# stringr           * 1.6.0      2025-11-04 [1] CRAN (R 4.3.3)
# survival            3.5-8      2024-02-14 [4] CRAN (R 4.3.2)
# tdigest             0.4.2      2026-03-05 [1] Github (hrbrmstr/tdigest@e52b746)
# TH.data             1.1-2      2023-04-17 [3] CRAN (R 4.3.1)
# tibble            * 3.3.1      2026-01-11 [1] CRAN (R 4.3.3)
# tictoc            * 1.2.1      2024-03-18 [1] CRAN (R 4.3.3)
# tidyr             * 1.3.0      2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect          1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse         * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
# timechange          0.2.0      2023-01-11 [1] CRAN (R 4.3.0)
# TMB                 1.9.4      2023-04-18 [1] CRAN (R 4.3.0)
# tweedie             2.3.5      2022-08-17 [1] CRAN (R 4.3.0)
# tzdb                0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
# vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.3.3)
# withr               3.0.2      2024-10-28 [1] CRAN (R 4.3.3)
# xtable              1.8-4      2019-04-21 [1] CRAN (R 4.3.1)
# zoo                 1.8-12     2023-04-13 [1] CRAN (R 4.3.1)
# 
# [1] /home/fhui28/R/x86_64-pc-linux-gnu-library/4.3
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────