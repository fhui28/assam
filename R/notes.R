#' Used in assam.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat)[-c(1:2)], collapse = "+")) %>% as.formula
    family <- binomial()
    y = simdat$y
    data = covariate_dat
    mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 50)
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    num_archetypes = 5
    num_cores = 8
    uncertainty_quantification <- TRUE
    control = list(max_iter = 500, tol = 1e-5, temper_prob = 0.8, trace = TRUE)
    bootstrap_control = list(num_boot = 100, ci_alpha = 0.05, seed = NULL, ci_type = "percentile")

    
    samfit <- assam(y = simdat$y,
                    formula = paste("~ ", paste0(colnames(covariate_dat)[-(1:2)], collapse = "+")) %>% as.formula,
                    data = covariate_dat,
                    mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 50),
                    family = binomial(),
                    do_parallel = FALSE,
                    uncertainty_quantification = TRUE,
                    num_archetypes = num_archetype,
                    num_cores = 8,
                    control = list(trace = 1),
                    bootstrap_control = list(num_boot = 20))
    }

#'
#'
#' Used in .quadapprox_fn.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat)[-(1:2)], collapse = "+")) %>% as.formula
    family <- binomial()
    resp = simdat$y
    data = covariate_dat
    mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 50)
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    num_nuisance_perspp = NULL
    start_archetype_labels = NULL
    start = NULL
    
    l = 30
    tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
    fit0_orig <- sdmTMB::sdmTMB(tmp_formula,
                            data = data.frame(response = simdat$y[,l], covariate_dat, trial_size = trial_size), 
                            spatial = TRUE,
                            mesh = mesh,
                            family = family)
    summary(fit0_orig)

    
    use_pars <- .get_pars(fit0_orig)
    use_pars[["b_j"]] <- rep(0, length(use_pars[["b_j"]]))
    #use_pars[["ln_phi"]] <- do_em$new_nuisance[l, grep("ln_phi", colnames(do_em$new_nuisance))]
    use_pars[["ln_tau_O"]] <- 1
    use_pars[["ln_kappa"]] <- matrix(1, 2, 1)

    use_map <- fit0_orig$tmb_map
    use_map$b_j <- as.factor(rep(NA, length(use_pars[["b_j"]])))
    use_map$ln_tau_O <- as.factor(NA)
    use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
    
    #' Method 1 of prediction/recovery of spatial random effects -- NEW
    new_tmb_obj <- TMB::MakeADFun(data = predict(fit0_orig, return_tmb_object = TRUE)$pred_tmb_data, #fit0_orig$tmb_data,
                                  profile = fit0_orig$control$profile,
                                  parameters = use_pars, 
                                  map = use_map,
                                  random = fit0_orig$tmb_random,
                                  DLL = "sdmTMB",
                                  silent = TRUE)
    
    new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
    new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
    r <- new_tmb_obj$report(new_tmb_obj$env$last.par) # last.par taken since it is the newest set of parameters
    
    
    #' #' Method 2 of prediction/recovery of spatial random effects -- OLD
    #' use_map <- list()
    #' use_map$b_j = as.factor(rep(NA, length(use_pars[["b_j"]])))
    #' use_map$ln_tau_O <- as.factor(NA)
    #' use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
    #' 
    #' new_tmb_obj2 <- TMB::MakeADFun(data = fit0_orig$tmb_data,
    #'                               profile = fit0_orig$control$profile,
    #'                               parameters = use_pars,
    #'                               map = use_map,
    #'                               random = fit0_orig$tmb_random,
    #'                               DLL = "sdmTMB",
    #'                               silent = TRUE)
    #' 
    #' new_fit0 <- nlminb(start = new_tmb_obj2$par,
    #'                    objective = new_tmb_obj2$fn,
    #'                    gradient = new_tmb_obj2$gr,
    #'                    control = fit0_orig$nlminb_control)
    #' 
    #' r2 <- new_tmb_obj2$report(new_tmb_obj2$env$last.par.best) #' Use TMB's report function to return the estimated field of assam

    #' Compare methods 1 and 2 
    # r %>% str
    # r2 %>% str
    # log(sqrt(8)/exp(use_pars$ln_kappa))
    # 
    # data.frame(true = simdat$spatial_fields[,l], r$omega_s_A, r2$omega_s_A)    
    }
    
