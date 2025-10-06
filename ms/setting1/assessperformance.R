#' ---
#' title: Simulation study for asSAMs -- assessment of performance
#' abstract: For setting 1 concerning estimation and uncertainty quantification
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
library(sdmTMB)
library(assam)


##----------------------
#' # Simulation parameters
##----------------------
load(file = "simsetting1_truemodel.RData")

num_X <- 9
num_spp <- length(true_spp_effects)
num_archetype <- nrow(true_betas)

new_formula <- ~ GBR_BATHY + GBR_TS_BSTRESS + GA_CRBNT + GA_GRAVEL + GA_MUD + CRS_O2_AV + CRS_S_AV + CRS_T_AV + SW_CHLA_AV



##----------------------
#' # Source results
##----------------------
response_types <- c("tweedie")
n_seq <- c(250, 500, 1000, 2000)
num_dataset <- 200
method_names <- c("assam_pointestimateonly", "assam_bootfast", "speciesmix_pointestimateonly")

all_computation_times <- all_RMSE <- all_bias <- all_coverage <- array(NA, 
                                                                       dim = c(length(n_seq), num_dataset, length(method_names)), 
                                                                       dimnames = list(n_seq, 1:num_dataset, method_names))



for(k0 in 1:length(n_seq)) {
    for(k1 in 1:num_dataset) {
        
        #' ## (p)asSAM results
        cw_filename <- paste0("simdat_", response_types, "_numunits", n_seq[k0], "_assam_dataset", k1, ".RData")
        if(!file.exists(cw_filename))
            next;
        load(cw_filename)
        
        all_computation_times[k0, k1, 1:2] <- c(alltimes_tictoclog_lst[[1]]$toc - alltimes_tictoclog_lst[[1]]$tic + alltimes_tictoclog_lst[[2]]$toc - alltimes_tictoclog_lst[[2]]$tic,
                                                alltimes_tictoclog_lst[[1]]$toc - alltimes_tictoclog_lst[[1]]$tic + alltimes_tictoclog_lst[[3]]$toc - alltimes_tictoclog_lst[[3]]$tic)
        
        get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, apply(samfit_final_pointest$posterior_probability, 1, which.max)) %>% apply(., 1, which.max)))
        all_bias[k0, k1, 1:2] <- c(mean(samfit_final_pointest$betas[get_pairs$pairs[,2],] - true_betas),
                                   mean(samfit_final_boot_fast$betas[get_pairs$pairs[,2],] - true_betas))
        all_RMSE[k0, k1, 1:2] <- c(sqrt(mean((samfit_final_pointest$betas[get_pairs$pairs[,2],] - true_betas)^2)),
                                   sqrt(mean((samfit_final_boot_fast$betas[get_pairs$pairs[,2],] - true_betas)^2)))
        all_coverage[k0, k1, 1:2] <- c(NA,
                                       mean(samfit_final_boot_fast$confidence_intervals$betas$upper[get_pairs$pairs[,2],] > true_betas & 
                                                samfit_final_boot_fast$confidence_intervals$betas$lower[get_pairs$pairs[,2],] < true_betas))

        rm(list = ls(pattern = "tictoclog"))
        rm(list = ls(pattern = "samfit"))
        
        
        #' ## species_mix results
        if(response_types != "tweedie") {
            all_computation_times[k0, k1, 3] <- c(alltimes_tictoclog_lst[[4]]$toc - alltimes_tictoclog_lst[[4]]$tic)
            
            get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, apply(speciesmix_pointestimate$tau, 1, which.max)) %>% apply(., 1, which.max)))
            all_bias[k0, k1, 3] <- c(mean(speciesmix_pointestimate$coefs$beta[get_pairs$pairs[,2],] - true_betas))
            all_RMSE[k0, k1, 3] <- c(sqrt(mean((speciesmix_pointestimate$coefs$beta[get_pairs$pairs[,2],] - true_betas)^2)))
            }
        
        rm(speciesmix_pointestimate)
        }
    }





apply(all_bias, c(1,3), mean, na.rm = TRUE) %>% round(4)
apply(all_RMSE, c(1,3), mean, na.rm = TRUE) %>% round(4)
apply(all_coverage, c(1,3), mean, na.rm = TRUE) %>% round(4)
apply(all_computation_times, c(1,3), mean, na.rm = TRUE) %>% round(4)


