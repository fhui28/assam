#' General testing
function(){
    library(sdmTMB)
    library(glmmTMB)
    library(mgcv)
    
    sdmfit <- sdmTMB(formula = paste("response ~ ", paste0("s(", colnames(covariate_dat)[1:3], ")", collapse = "+"), "+", paste0(colnames(covariate_dat)[4:7], collapse = "+")) %>% as.formula,
                     data = data.frame(response = simdat$y[,1], covariate_dat),
                     family = nbinom2(),
                     spatial = FALSE)
    
    summary(sdmfit)
    sdmfit_covariance <- sdmfit$sd_report$jointPrecision %>% solve

    
    gamfit <- gam(formula = paste("response ~ ", paste0("s(", colnames(covariate_dat)[1:3], ")", collapse = "+"), "+", paste0(colnames(covariate_dat)[4:7], collapse = "+")) %>% as.formula,
                 data = data.frame(response = simdat$y[,1], covariate_dat),
                 family = nb(),
                 method = "ML")
    
    summary(gamfit)
    summary(sdmfit)
    vcov(gamfit) %>% diag %>% sqrt
    sdmfit_covariance %>% diag %>% sqrt #' Need to rearrange such that all the smoothers are put together?
    
    }


#' Used in assam.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula
    family <- tweedie()
    y = simdat$y
    data = covariate_dat
    mesh = NULL #sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 50)
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
    family <- family
    resp = simdat$y
    data = covariate_dat
    mesh = mesh
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    num_nuisance_perspp = NULL
    start_archetype_labels = NULL
    start = NULL
    
    
    l = 1
    tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
    fit0_orig <- sdmTMB(tmp_formula,
                   data = data.frame(response = resp[,l], data, trial_size = trial_size), 
                   spatial =  TRUE,
                   mesh = mesh,
                   family = family)
    

    use_pars <- fit0_orig$tmb_obj$env$parList()
    use_pars[["b_j"]] <- c(do_em$new_spp_intercept[l], do_em$new_betas[which.max(do_em$post_prob[l,]),])
    use_pars[["ln_phi"]] <- do_em$new_nuisance[l, grep("ln_phi", colnames(do_em$new_nuisance))]
    use_pars[["ln_tau_O"]] <- do_em$new_nuisance[l, grep("ln_tau_O", colnames(do_em$new_nuisance))]
    use_pars[["ln_kappa"]] <- matrix(do_em$new_nuisance[l, grep("ln_kappa", colnames(do_em$new_nuisance))], 2, 1)
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
                       gradient = new_tmbobj$gr)

    new_fit0 <- sdmTMB(tmp_formula,
                        data = data.frame(response = resp[,l], data, trial_size = trial_size), 
                        spatial =  TRUE,
                        mesh = mesh,
                        control = sdmTMB::sdmTMBcontrol(
                            start = list(
                                b_j = c(do_em$new_spp_intercept[l], do_em$new_betas[which.max(do_em$post_prob[l,]),]),
                                ln_phi = do_em$new_nuisance[l, grep("ln_phi", colnames(do_em$new_nuisance))],
                                ln_tau_O = do_em$new_nuisance[l, grep("ln_tau_O", colnames(do_em$new_nuisance))],
                                ln_kappa = matrix(do_em$new_nuisance[l, grep("ln_kappa", colnames(do_em$new_nuisance))], 2, 1)),
                            map = list(b_j = as.factor(rep(NA, 1+num_X)), ln_phi = as.factor(NA), ln_tau_O = as.factor(NA), ln_kappa = as.factor(matrix(NA, 2, 1)))
                            ),
                       do_fit = TRUE,
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


