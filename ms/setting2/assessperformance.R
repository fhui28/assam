#' ---
#' title: Simulation study for asSAMs -- assessing performance of model selection
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
library(ROCR)
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
true_powerparam <- runif(num_spp, 1.4, 1.8)
true_mixprop <- c(0.2, 0.2, 0.35, 0.15, 0.1)


##----------------------
#' # Source results
##----------------------
response_types <- c("binomial")
n_seq <- c(250, 500, 1000, 2000)
num_dataset <- 200
method_names <- c("assam_orderselect", "samfit_orderbetaselect", 
                  "stackedglmnet_PAM_chooseK", "stackedglmnet_PAM_knownclusters", 
                  "speciesmix")

all_computation_times <- all_order_selection <- array(NA, 
                                                      dim = c(length(n_seq), num_dataset, length(method_names)), 
                                                      dimnames = list(n_seq, 1:num_dataset, method_names))
all_sensitivity <- all_specificity <- all_F1score <- all_phi <- array(NA, 
                                                                      dim = c(length(n_seq), num_dataset, length(method_names[c(2,4)])), 
                                                                      dimnames = list(n_seq, 1:num_dataset, method_names[c(2,4)]))


for(k0 in 1:length(n_seq)) {
    for(k1 in 1:num_dataset) {
        
        cw_filename <- paste0("simdat_", response_types, "_numunits", n_seq[k0], "_assam_dataset", k1, ".RData")
        if(!file.exists(cw_filename))
            next;
        load(cw_filename)
        
        #' ## (p)asSAM results
        all_computation_times[k0, k1, 1:2] <- c(alltimes_tictoclog_lst[[1]]$toc-alltimes_tictoclog_lst[[1]]$tic + 
                                                    alltimes_tictoclog_lst[[2]]$toc-alltimes_tictoclog_lst[[2]]$tic,
                                                alltimes_tictoclog_lst[[1]]$toc-alltimes_tictoclog_lst[[1]]$tic + 
                                                    alltimes_tictoclog_lst[[2]]$toc-alltimes_tictoclog_lst[[2]]$tic +
                                                    alltimes_tictoclog_lst[[3]]$toc-alltimes_tictoclog_lst[[3]]$tic + 
                                                    alltimes_tictoclog_lst[[4]]$toc-alltimes_tictoclog_lst[[4]]$tic)
        
        all_order_selection[k0, k1, 1:2] <- c(samfit_betaselect$num_archetypes,
                                              samfit_betaselect$num_archetypes)
        
        
        if(samfit_final_BIC_BIC2$num_archetypes == 5) {
            get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, apply(samfit_final_BIC_BIC2$posterior_probability, 1, which.max)) %>% apply(., 1, which.max)))
            
            makepred_obj <- ROCR::prediction(predictions = as.vector(1*(samfit_final_BIC_BIC2$betas[get_pairs$pairs[,2],] != 0)), labels = as.vector(1*(true_betas != 0)))
            all_sensitivity[k0,k1,1] <- ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2]
            all_specificity[k0,k1,1] <- ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2]
            all_F1score[k0,k1,1] <- ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2]
            all_phi[k0,k1,1] <- ROCR::performance(makepred_obj, measure = "phi")@y.values[[1]][2]
            }
        
        rm(list = ls(pattern = "samfit"))
        
        
        #' ## Stacked GLMnet and speciesmix results
        if(response_types != "tweedie") {
            all_computation_times[k0, k1, 3:5] <- c(alltimes_tictoclog_lst[[5]]$toc - alltimes_tictoclog_lst[[5]]$tic +
                                                        alltimes_tictoclog_lst[[6]]$toc - alltimes_tictoclog_lst[[6]]$tic,
                                                    alltimes_tictoclog_lst[[5]]$toc - alltimes_tictoclog_lst[[5]]$tic +
                                                        alltimes_tictoclog_lst[[7]]$toc - alltimes_tictoclog_lst[[7]]$tic,
                                                    alltimes_tictoclog_lst[[8]]$toc - alltimes_tictoclog_lst[[8]]$tic)
            
            all_order_selection[k0, k1, 3:5] <- c(nrow(stacked_glmnet_pam_gap$medoids),
                                                  nrow(stacked_glmnet_pam_knownclusters$medoids),
                                                  nrow(speciesmix_chosen$coefs$beta))
            
            
            get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, stacked_glmnet_pam_knownclusters$clustering) %>% apply(., 1, which.max)))
            makepred_obj <- ROCR::prediction(predictions = as.vector(1*(stacked_glmnet_pam_knownclusters$medoids[get_pairs$pairs[,2],] != 0)), labels = as.vector(1*(true_betas != 0)))
            all_sensitivity[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2]
            all_specificity[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2]
            all_F1score[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2]
            all_phi[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "phi")@y.values[[1]][2]
            }
        
        if(response_types == "tweedie") {
            all_computation_times[k0, k1, 3:4] <- c(alltimes_tictoclog_lst[[5]]$toc - alltimes_tictoclog_lst[[5]]$tic +
                                                        alltimes_tictoclog_lst[[6]]$toc - alltimes_tictoclog_lst[[6]]$tic,
                                                    alltimes_tictoclog_lst[[5]]$toc - alltimes_tictoclog_lst[[5]]$tic +
                                                        alltimes_tictoclog_lst[[7]]$toc - alltimes_tictoclog_lst[[7]]$tic)
            
            all_order_selection[k0, k1, 3:4] <- c(nrow(stacked_glmnet_pam_gap$medoids),
                                                  nrow(stacked_glmnet_pam_knownclusters$medoids))
            
            
            get_pairs <- list(pairs = cbind(1:num_archetype, table(true_archetype_label, stacked_glmnet_pam_knownclusters$clustering) %>% apply(., 1, which.max)))
            makepred_obj <- ROCR::prediction(predictions = as.vector(1*(stacked_glmnet_pam_knownclusters$medoids[get_pairs$pairs[,2],] != 0)), labels = as.vector(1*(true_betas != 0)))
            all_sensitivity[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "sens")@y.values[[1]][2]
            all_specificity[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "spec")@y.values[[1]][2]
            all_F1score[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "f")@y.values[[1]][2]
            all_phi[k0,k1,2] <- ROCR::performance(makepred_obj, measure = "phi")@y.values[[1]][2]
            }
        
        rm(list = ls(pattern = "stacked_glmnet"))
        rm(speciesmix_chosen)
        rm(true_archetype_label)
        }
    }





all_order_selection %>% 
    as.data.frame.table %>% 
    group_by(Var1, Var3) %>%
    reframe(prop_correct = sum(Freq == 5, na.rm = TRUE),
            num_converged_datasets = sum(!is.na(Freq))) %>% 
    print(n = Inf)

# apply(all_sensitivity, c(1,3), mean, na.rm = TRUE) %>% round(4)
# apply(all_specificity, c(1,3), mean, na.rm = TRUE) %>% round(4)
apply(all_F1score, c(1,3), mean, na.rm = TRUE) %>% round(4)
# apply(all_phi, c(1,3), mean, na.rm = TRUE) %>% round(4)

apply(all_computation_times, c(1,3), mean, na.rm = TRUE) %>% round(4)
