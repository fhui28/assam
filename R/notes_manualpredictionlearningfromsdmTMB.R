rm(list = ls())
library(sdmTMB)
library(tidyverse)

spatial_train <- 1:500
pcod_2011_train <- pcod_2011[spatial_train,]
pcod_2011_test <- pcod_2011[-spatial_train,]

mesh <- make_mesh(pcod_2011_train, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
m <- sdmTMB(formula = density ~ 1 + depth_scaled + depth_scaled2,
            data = pcod_2011_train, 
            spatial = FALSE,
            mesh = mesh, 
            family = tweedie())

summary(m)

predictions_gold <- predict(m, 
                       newdata = pcod_2011_test)
str(predictions_gold)


#' # Testing manual method for prediction
use_pars <- .get_pars(m)
use_pars[["b_j"]] <- rnorm(length(use_pars[["b_j"]]))
use_pars[["ln_phi"]] <- 0
use_pars[["thetaf"]] <- 0
use_pars[["ln_tau_O"]] <- 1
use_pars[["ln_kappa"]] <- matrix(1, 2, 1)

use_map <- list()
use_map$b_j <- as.factor(rep(NA, length(use_pars[["b_j"]])))
use_map$ln_phi <- as.factor(NA)
use_map$thetaf <- as.factor(NA)
use_map$ln_tau_O <- as.factor(NA)
use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 


#' Method 1 of prediction/recovery of spatial random effects -- NEW
new_tmb_obj <- TMB::MakeADFun(data = predict(m, newdata = pcod_2011_test, return_tmb_object = TRUE)$pred_tmb_data,
                              profile = m$control$profile,
                              parameters = use_pars, 
                              map = use_map,
                              random = m$tmb_random,
                              DLL = "sdmTMB",
                              silent = TRUE)

new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
r <- new_tmb_obj$report(new_tmb_obj$env$last.par) # last.par taken since it is the newest set of parameters

data.frame(tmb_predict = r$proj_eta[,1],
           manual_predict = model.matrix(m$formula[[1]], data = pcod_2011_test) %*% use_pars$b_j)
           #predictions_gold$est)



