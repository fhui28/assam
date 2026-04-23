#' ---
#' title: Simulation study for asSAMs -- assessment of performance
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
library(colorspace)
library(ROCR)
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
#' # Source results
##----------------------
response_types <- c("binomial")
n_seq <- c(250, 500)
num_dataset <- 200
method_names <- c("assam_pointestimateonly", "assam_bootfast", "speciesmix_pointestimateonly", "speciesmix_boot", "stan")

all_computation_times <- all_RMSE <- all_bias <- all_coverage <- array(NA, 
                                                                       dim = c(length(n_seq), num_dataset, length(method_names)), 
                                                                       dimnames = list(n_seq, 1:num_dataset, method_names))



for(k0 in 1:length(n_seq)) {
    for(k1 in 1:num_dataset) {
        cw_filename <- paste0("simdat_numunits", n_seq[k0], "_dataset", k1, ".RData")
        if(!file.exists(cw_filename))
            next;
        load(cw_filename)
        
        
        #' ## asSAM results
        all_computation_times[k0, k1, 1:2] <- c(assams_tictoclog_lst[[1]]$toc - assams_tictoclog_lst[[1]]$tic + assams_tictoclog_lst[[2]]$toc - assams_tictoclog_lst[[2]]$tic,
                                                assams_tictoclog_lst[[1]]$toc - assams_tictoclog_lst[[1]]$tic + assams_tictoclog_lst[[3]]$toc - assams_tictoclog_lst[[3]]$tic)
        
        get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, apply(samfit_final_pointest$posterior_probability, 1, which.max)) %>% apply(., 1, which.max)))
        all_bias[k0, k1, 1:2] <- c(mean(samfit_final_pointest$betas[get_pairs$pairs[,2],] - true_betas),
                                   mean(samfit_final_boot_fast$betas[get_pairs$pairs[,2],] - true_betas))
        all_RMSE[k0, k1, 1:2] <- c(sqrt(mean((samfit_final_pointest$betas[get_pairs$pairs[,2],] - true_betas)^2)),
                                   sqrt(mean((samfit_final_boot_fast$betas[get_pairs$pairs[,2],] - true_betas)^2)))
        all_coverage[k0, k1, 1:2] <- c(NA,
                                       mean(samfit_final_boot_fast$confidence_intervals$betas$upper[get_pairs$pairs[,2],] > true_betas & 
                                                samfit_final_boot_fast$confidence_intervals$betas$lower[get_pairs$pairs[,2],] < true_betas))
        rm(list = ls(pattern = "samfit"))
        
        
        #' ## species_mix results
        all_computation_times[k0, k1, 3:4] <- c(speciesmix_tictoclog_lst[[1]]$toc - speciesmix_tictoclog_lst[[1]]$tic,
                                                speciesmix_tictoclog_lst[[1]]$toc - speciesmix_tictoclog_lst[[1]]$tic + speciesmix_tictoclog_lst[[2]]$toc - speciesmix_tictoclog_lst[[2]]$tic)
        get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, apply(speciesmix_pointestimate$tau, 1, which.max)) %>% apply(., 1, which.max)))
        all_bias[k0, k1, 3:4] <- c(mean(speciesmix_pointestimate$coefs$beta[get_pairs$pairs[,2],] - true_betas),
                                   mean(speciesmix_pointestimate$coefs$beta[get_pairs$pairs[,2],] - true_betas))
        all_RMSE[k0, k1, 3:4] <- c(sqrt(mean((speciesmix_pointestimate$coefs$beta[get_pairs$pairs[,2],] - true_betas)^2)),
                                   sqrt(mean((speciesmix_pointestimate$coefs$beta[get_pairs$pairs[,2],] - true_betas)^2)))
        form_cis <- list(lower = (speciesmix_pointestimate$summary$Estimate - qnorm(0.975)*speciesmix_pointestimate$summary$SE) %>% matrix(nrow = num_archetype, ncol = num_X),
                         upper = (speciesmix_pointestimate$summary$Estimate + qnorm(0.975)*speciesmix_pointestimate$summary$SE) %>% matrix(nrow = num_archetype, ncol = num_X))
        all_coverage[k0, k1, 3:4] <- c(NA,
                                       mean(form_cis$upper[get_pairs$pairs[,2],] > true_betas & form_cis$lower[get_pairs$pairs[,2],] < true_betas))
        rm(speciesmix_pointestimate, form_cis)
        
        
        #' ## stan results
        all_computation_times[k0, k1, 5] <- stanfit_tictoclog_lst[[1]]$toc - stanfit_tictoclog_lst[[1]]$tic
        get_pairs <- RcppHungarian::HungarianSolver(-tcrossprod(true_betas, fit_stan_summary$beta$mean %>% matrix(nrow = num_archetype, ncol = num_X))) #' Need to use Hungarian algorithm to match archetypes since posterior probabilities are not saved/easily available for the stanfit object
        fit_stan_summary$estimate$beta <- matrix(fit_stan_summary$beta$mean, nrow = num_archetype, ncol = num_X)
        all_bias[k0, k1, 5] <- mean(fit_stan_summary$estimate$beta[get_pairs$pairs[,2],] - true_betas)
        all_RMSE[k0, k1, 5] <- sqrt(mean((fit_stan_summary$estimate$beta[get_pairs$pairs[,2],] - true_betas)^2))
        form_cis <- list(lower = fit_stan_summary$beta$`2.5%` %>% matrix(nrow = num_archetype, ncol = num_X),
                         upper = fit_stan_summary$beta$`97.5%` %>% matrix(nrow = num_archetype, ncol = num_X))
        all_coverage[k0, k1, 5] <- mean(form_cis$upper[get_pairs$pairs[,2],] > true_betas & form_cis$lower[get_pairs$pairs[,2],] < true_betas)
        rm(fit_stan_summary, form_cis)
        }
    }



apply(all_bias*10, c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_RMSE, c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_coverage, c(1,3), mean, na.rm = TRUE) %>% round(3)
apply(all_computation_times, c(1,3), mean, na.rm = TRUE) %>% round(3)

