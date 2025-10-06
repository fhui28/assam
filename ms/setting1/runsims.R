#' ---
#' title: Simulation study for asSAMs
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
#' # Set up simulation function
##----------------------
simfn <- function(seed,
                  family = binomial(),
                  num_units = 500) {
    
    if(!family$family %in% c("poisson", "binomial", "nbinom2", "tweedie")) {
        stop("Family is currently not permitted...sorry!")
        }
    
    
    num_boot <- 100
    
    ##----------------------
    # Generate some multivariate abundance (count) data from a SAM
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
    
    tictoc::tic.clearlog()
    
    ##----------------------
    #' # Method 1: asSAMs
    ##----------------------
    #' ## Construct the initial stacked species model fits that will be used throughout the model selection process below
    tictoc::tic("assams_prefit")
    samfit_prefit <- assam(y = simdat$y,
                           formula = new_formula,
                           data = covariate_dat,
                           family = family,
                           num_archetypes = 2, #' This is arbitrary and does not matter
                           num_cores = detectCores() - 4,
                           do_assam_fit = FALSE)
    tictoc::toc(log = TRUE)
    

    tictoc::tic("assams_final_pointestimate_only")
    samfit_final_pointest <- try(assam(y = simdat$y,
                                       formula = new_formula,
                                       data = covariate_dat,
                                       family = family,
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
                                    family = family,
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
        
    rm(samfit_prefit)
    gc()

    
    ##----------------------
    #' # Method 2: species_mix
    #' Note (bootstrap) standard errors are not attempted as computationally they are too expensive. 
    #' Currently does not do Tweedie due to a bug (https://github.com/skiptoniam/ecomix/issues/27) as of April 2025
    ##----------------------
    if(family$family != "tweedie") {
        sam_form <- paste0('cbind(',paste("spp", 1:num_spp, sep = "", collapse = ','), ") ~ ", as.character(new_formula)[2]) %>% 
            as.formula
        
        pick_response_type <- NULL
        if(family$family == "binomial")
            pick_response_type <- "bernoulli"
        if(family$family == "nbinom2")
            pick_response_type <- "negative.binomial"
        # if(family$family == "tweedie")
        #     pick_response_type <- "tweedie"

        tictoc::tic("speciesmix_pointestimateonly")
        speciesmix_pointestimate <- try(species_mix(archetype_formula = sam_form,
                                                    species_formula = ~ 1,
                                                    data = data.frame(simdat$y, covariate_dat), 
                                                    family = pick_response_type, 
                                                    nArchetypes = num_archetype),
                                        silent = TRUE)
        tictoc::toc(log = TRUE)
        
        speciesmix_pointestimate$titbits <- NULL
        speciesmix_pointestimate$mus <- NULL
        speciesmix_pointestimate$terms <- NULL
        }
    
    if(family$family == "tweedie") {
        speciesmix_pointestimate <- NA
        }
    
    
    ##----------------------
    #' # Final return
    ##----------------------
    gc()
    
    alltimes_tictoclog_txt <- tictoc::tic.log(format = TRUE)
    alltimes_tictoclog_lst <- tictoc::tic.log(format = FALSE)
    true_archetype_label <- simdat$archetype_label
    
    save(samfit_final_pointest,
         samfit_final_boot_fast,
         speciesmix_pointestimate,
         alltimes_tictoclog_txt,
         alltimes_tictoclog_lst,
         true_archetype_label,
         file = paste0("simdat_", family$family, "_numunits", num_units, "_dataset", seed, ".RData"))
    
    tictoc::tic.clearlog()
    }



##----------------------
#' # Run the simulation
##----------------------
registerDoParallel(cores = detectCores() - 5)

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
    

