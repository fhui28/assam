#' Used in assam.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat)[-c(1:2)], collapse = "+")) %>% as.formula
    family <- nbinom2()
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
    }

#'
#'
#' Used in .quadapprox_fn.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat)[-(1:2)], collapse = "+")) %>% as.formula
    family <- nbinom2()
    resp = simdat$y
    data = covariate_dat
    mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 50)
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    num_nuisance_perspp = NULL
    start_archetype_labels = NULL
    start = NULL
    
    tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
    fit0_orig <- sdmTMB::sdmTMB(tmp_formula,
                            data = data.frame(response = simdat$y[,l], covariate_dat, trial_size = trial_size), 
                            spatial = TRUE,
                            mesh = mesh,
                            family = nbinom2())


    use_pars <- fit0_orig$tmb_obj$env$parList()
    use_pars[["b_j"]] <- rnorm(use_pars[["b_j"]]) #c(do_em$new_spp_intercept[l], do_em$new_betas[which.max(do_em$post_prob[l,]),])
    #use_pars[["ln_phi"]] <- do_em$new_nuisance[l, grep("ln_phi", colnames(do_em$new_nuisance))]
    #use_pars[["ln_tau_O"]] <- do_em$new_nuisance[l, grep("ln_tau_O", colnames(do_em$new_nuisance))]
    #use_pars[["ln_kappa"]] <- matrix(do_em$new_nuisance[l, grep("ln_kappa", colnames(do_em$new_nuisance))], 2, 1)
    new_tmbobj <- TMB::MakeADFun(
        data = fit0_orig$tmb_data,
        profile = fit0_orig$control$profile,
        parameters = use_pars, 
        map = list(b_j = as.factor(rep(NA, 1+num_X)), ln_phi = as.factor(NA), ln_tau_O = as.factor(NA), ln_kappa = as.factor(matrix(NA, 2, 1))),
        random = fit0_orig$tmb_random,
        DLL = "sdmTMB",
        silent = TRUE
        )
    new_fit0 <- nlminb(start = new_tmbobj$par,
                       objective = new_tmbobj$fn,
                       gradient = new_tmbobj$gr,
                       control = fit0_orig$nlminb_control)
    r <- new_tmbobj$report(new_tmbobj$env$last.par.best)
    cbind(predict(fit0_orig)$omega_s, r$omega_s_A)
    
    
    new_fit0 <- sdmTMB(tmp_formula,
                       data = data.frame(response = simdat$y[,l], covariate_dat, trial_size = trial_size), 
                       spatial =  TRUE,
                        mesh = mesh,
                        control = sdmTMBcontrol(
                            start = list(
                                 b_j = as.vector(c(do_em$new_spp_intercept[l], do_em$new_betas[which.max(do_em$post_prob[l,]),])),
                                 ln_phi = 1, #do_em$new_nuisance[l, grep("ln_phi", colnames(do_em$new_nuisance))],
                                 ln_tau_O = 0.5, #do_em$new_nuisance[l, grep("ln_tau_O", colnames(do_em$new_nuisance))],
                                 ln_kappa = matrix(1, 2, 1)), #do_em$new_nuisance[l, grep("ln_kappa", colnames(do_em$new_nuisance))]
                            map = list(ln_phi = as.factor(NA), ln_tau_O = as.factor(NA), ln_kappa = as.factor(matrix(NA, 2, 1)), b_j = as.factor(rep(NA, 1+num_X))),
                            lower = list(b_j = as.vector(c(do_em$new_spp_intercept[l], do_em$new_betas[which.max(do_em$post_prob[l,]),]) - 1e-1)),
                            upper = list(b_j = as.vector(c(do_em$new_spp_intercept[l], do_em$new_betas[which.max(do_em$post_prob[l,]),]) + 1e-1)),
                            ),
                       family = family)
}
    
#'
#'
#' Used in predict.assam.R
function() {
    object <- testfit
    newdata = covariate_dat
    newoffset = NULL
    type = "archetype"
    se_fit = TRUE
    num_cores = 8
    coverage = 0.95
}


