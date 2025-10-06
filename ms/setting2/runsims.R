#' ---
#' title: Simulation study for asSAMs
#' abstract: For setting 2 concerning model selection 
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


##----------------------
#' # Simulation parameters
##----------------------
num_X <- 12
num_spp <- 100
num_archetype <- 5

set.seed(052025)
true_betas <- runif(num_archetype * num_X, -1, 1) %>% matrix(nrow = num_archetype)
true_betas[which(abs(true_betas) < 0.2)] <- 0 # Making archetypal coefficients sparse
true_betas[,9:12] <- 0
true_betas

true_spp_effects <- matrix(runif(num_spp, -3, -1), ncol = 1)
true_dispparam <- 1/runif(num_spp, 1, 5)
true_powerparam <- runif(num_spp, 1.2, 1.8)
true_mixprop <- c(0.2, 0.2, 0.35, 0.15, 0.1)


##----------------------
#' # Set up simulation function
##----------------------
simfn <- function(seed,
                  family = binomial(),
                  num_units = 500) {
    
    if(!family$family %in% c("poisson", "binomial", "nbinom2", "tweedie")) {
        stop("Family is currently not permitted...sorry!")
        }
    
    ##----------------------
    # Generate some multivariate abundance (count) data from a sparse SAM
    ##----------------------
    set.seed(052025)
    H <- outer(1:num_X, 1:num_X, "-")
    H <- 0.7^abs(H)
    covariate_dat <- rmvnorm(num_units, sigma = H) %>%
        as.data.frame %>%
        rename_with(., .fn = function(x) paste0("covariate", x))
    rm(H)
    set.seed(NULL)
    
    simdat <- create_samlife(family = family,
                             formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                             data = covariate_dat,
                             betas = true_betas,
                             spp_effects = true_spp_effects,
                             spp_dispparam = true_dispparam,
                             spp_powerparam = true_powerparam,
                             mixing_proportion = true_mixprop,
                             seed = seed)

    tictoc::tic.clearlog()
    
    ##----------------------
    #' # Method 1: pasSAMs
    #' Note this includes both pasSAMs that perform order selection only and order selection + archetypal regression coefficient selection
    ##----------------------
    #' ## Construct the initial stacked species model fits that will be used throughout the model selection process below
    tictoc::tic("passams_prefit")
    samfit_prefit <- assam(y = simdat$y,
                           formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                           data = covariate_dat,
                           family = family,
                           num_archetypes = 2, #' This is arbitrary and does not matter
                           num_cores = detectCores() - 4,
                           do_assam_fit = FALSE)
    tictoc::toc(log = TRUE)
    
    
    #' ## asSAMs/pasSAMs for order selection only 
    tictoc::tic("passams_orderselect")
    assam_fitK <- function(k) {
        out <- assam(y = simdat$y,
                     formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                     data = covariate_dat,
                     family = family,
                     num_archetypes = k,
                     uncertainty_quantification = FALSE,
                     supply_quadapprox = samfit_prefit,
                     num_cores = 1)
        
        out$sdmTMB_fits <- NULL
        return(out)
        }
    
    registerDoParallel(cores = detectCores() - 4)
    all_assamK <- foreach(k = 1:10) %dopar% assam_fitK(k = k)
    samfit_orderselect <- list(BIC = sapply(all_assamK, function(x) -2*x$logL + log(num_spp)*(length(x$mixing_proportion)-1 + length(x$betas)))) 
    tictoc::toc(log = TRUE)    
    
    select_num_archetypes <- which.min(samfit_orderselect$BIC)
    
    #' ## Now perform selection on the archetypal regression coefficients
    #' Minimum tuning parameter is such that there are at least five non-zero coefficients
    tictoc::tic("passams_betaselect")
    samfit_betaselect <- try(passam(y = simdat$y,
                                    formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                                    data = covariate_dat,
                                    family = family,
                                    num_archetypes = select_num_archetypes,
                                    selection_on = "betas",
                                    supply_quadapprox = samfit_prefit,
                                    num_cores = detectCores() - 4,
                                    beta_selection_control = list(min_df = 5)),
                             silent = TRUE)
    tictoc::toc(log = TRUE)

    
    #' ## Now fit the final asSAMs given a chosen value of the tuning parameter
    choose_lambda <- which.min(samfit_betaselect$BIC2)
    
    tictoc::tic("passams_final_BIC_BIC2")
    samfit_final_BIC_BIC2 <- assam(y = simdat$y,
                                   formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                                   data = covariate_dat,
                                   family = family,
                                   beta_selection = TRUE,
                                   num_archetypes = select_num_archetypes,
                                   uncertainty_quantification = FALSE,
                                   supply_quadapprox = samfit_prefit,
                                   beta_selection_control = list(lambda = samfit_betaselect$lambda[choose_lambda],
                                                                 warm_start = list(betas = samfit_betaselect$betas_path[,,choose_lambda-1],
                                                                                   spp_effects = samfit_betaselect$spp_effects_path[,,choose_lambda-1],
                                                                                   spp_nuisance = samfit_betaselect$spp_nuisance_path[,,choose_lambda-1],
                                                                                   mixing_proportions = samfit_betaselect$mixing_proportions_path[,choose_lambda-1],
                                                                                   posterior_probability = samfit_betaselect$posterior_probability_path[,,choose_lambda-1])),
                                   num_cores = detectCores() - 4)
    samfit_final_BIC_BIC2$sdmTMB_fits <- samfit_final_BIC_BIC2$linear_predictor <- NULL
    tictoc::toc(log = TRUE)
    
    rm(samfit_prefit)
    samfit_betaselect$posterior_probability_path <- NULL
    
    
    ##----------------------
    #' # Method 2: Stacked penalized glms using glmnet, using 10-fold cross-validation with deviance as the loss function and selecting lambda.1se. Then apply K-means or K-medoids assuming either the number of clusters is unknown, or unknown as determined using the gap statistic
    ##----------------------
    registerDoParallel(cores = detectCores() - 4)
    
    tictoc::tic("stacked_glmnet")
    MM <- model.matrix(paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                       data = covariate_dat)[,-1]
    stacked_glmnet_coefs <- matrix(0, nrow = num_spp, ncol = ncol(MM)+1)
    for(k1 in 1:num_spp) {
        message("glmnet onto species ", k1)
        
        make_family <- family
        if(family$family == "tweedie")
            make_family <- statmod::tweedie(var.power = true_powerparam[k1], link.power = 0)
        if(family$family == "nbinom2")
            make_family <- MASS::negative.binomial(theta = 1/true_dispparam[k1])
        
        find_cv <- glmnet::cv.glmnet(y = simdat$y[,k1], 
                                     x = MM, 
                                     parallel = TRUE, 
                                     type.measure = "default", 
                                     family = make_family)
        stacked_glmnet_coefs[k1,] <- as.vector(coef(find_cv, s = "lambda.1se", exact = TRUE))
        }
    tictoc::toc(log = TRUE)
    rm(find_cv)
    
    
    #' ### Perform K-medoids + gap statistic
    tictoc::tic("Gap statistic with K-medoids on stacked glmnet coefficients")
    select_num_archetypes <- clusGap(stacked_glmnet_coefs[,-1],
                                     FUN = pam, 
                                     K.max = 10,
                                     B = 500)
    get_clusters <- maxSE(select_num_archetypes$Tab[,"gap"],
                          select_num_archetypes$Tab[,"SE.sim"],
                          method = "firstSEmax")
    stacked_glmnet_pam_gap <- pam(stacked_glmnet_coefs[,-1], k = get_clusters) 
    tictoc::toc(log = TRUE)
    
    #' ### Perform K-medoids assuming number of clusters is known
    tictoc::tic("K-medoids on stacked glmnet coefficients, assuming number of archetypes is known") #' nstart = 25 not needed given the help file info
    stacked_glmnet_pam_knownclusters <- pam(stacked_glmnet_coefs[,-1], k = num_archetype)
    tictoc::toc(log = TRUE)
    
    rm(select_num_archetypes, get_clusters, MM)

    
    ##----------------------
    #' # Method 3: species_mix.multifit
    #' Since starting value used to same time for a simulation!
    #' Currently does not do Tweedie due to a bug (https://github.com/skiptoniam/ecomix/issues/27) as of May 2025
    ##----------------------
    if(family$family != "tweedie") {
        sam_form <- paste0('cbind(',paste("spp", 1:num_spp, sep = "", collapse = ','), ") ~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula
        
        pick_response_type <- NULL
        if(family$family == "binomial")
            pick_response_type <- "bernoulli"
        if(family$family == "nbinom2")
            pick_response_type <- "negative.binomial"
        if(family$family == "tweedie")
            pick_response_type <- "tweedie"
        
        tictoc::tic("speciesmix_multifit")
        speciesmix_fitK <- function(k) {
            speciesmix_chooseK <- try(species_mix(archetype_formula = sam_form,
                                                  species_formula = ~ 1,
                                                  data = data.frame(simdat$y, covariate_dat), 
                                                  family = pick_response_type, 
                                                  nArchetypes = k),
                                      silent = TRUE)
            
            return(speciesmix_chooseK)
            }
        
        registerDoParallel(cores = detectCores() - 4)
        all_speciesmix <- foreach(k = 1:10) %dopar% speciesmix_fitK(k = k)
        
        select_num_archetypes <- sapply(all_speciesmix, function(x) try(BIC(x), silent = TRUE))
        if(inherits(select_num_archetypes, "try-error")) {
            return()
            }     
        select_num_archetypes <- which.min(select_num_archetypes)
        speciesmix_chosen <- all_speciesmix[[select_num_archetypes]]
        speciesmix_chosen$titbits <- NULL
        speciesmix_chosen$terms<- NULL
        speciesmix_chosen$mus <- NULL
        tictoc::toc(log = TRUE)
        
        rm(all_speciesmix)
        }
    
    if(family$family == "tweedie") {
        speciesmix_chosen <- NA
        }
    
    
    ##----------------------
    #' # Final return
    ##----------------------
    gc()
    
    alltimes_tictoclog_txt <- tictoc::tic.log(format = TRUE)
    alltimes_tictoclog_lst <- tictoc::tic.log(format = FALSE)
    true_archetype_label <- simdat$archetype_label
    
    save(samfit_orderselect,
         samfit_betaselect,
         samfit_final_BIC_BIC2,
         stacked_glmnet_coefs,
         stacked_glmnet_pam_gap,
         stacked_glmnet_pam_knownclusters,
         speciesmix_chosen,
         alltimes_tictoclog_txt,
         alltimes_tictoclog_lst,
         true_archetype_label,
         file = paste0("simdat_", family$family, "_numunits", num_units, "_dataset", seed, ".RData"))
    }


##----------------------
#' # Run the simulation
##----------------------
for(k in 1:200)
    simfn(seed = k, num_units = 250, family = binomial())

for(k in 1:200)
    simfn(seed = k, num_units = 500, family = binomial())

for(k in 1:200)
    simfn(seed = k, num_units = 1000, family = binomial())

for(k in 1:200)
    simfn(seed = k, num_units = 2000, family = binomial())


for(k in 1:200)
    simfn(seed = k, num_units = 250, family = nbinom2())

for(k in 1:200)
    simfn(seed = k, num_units = 500, family = nbinom2())

for(k in 1:200)
    simfn(seed = k, num_units = 1000, family = nbinom2())

for(k in 1:200)
    simfn(seed = k, num_units = 2000, family = nbinom2())


for(k in 1:200)
    simfn(seed = k, num_units = 250, family = tweedie())

for(k in 1:200)
    simfn(seed = k, num_units = 500, family = tweedie())

for(k in 1:200)
    simfn(seed = k, num_units = 1000, family = tweedie())

for(k in 1:200)
    simfn(seed = k, num_units = 2000, family = tweedie())



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
