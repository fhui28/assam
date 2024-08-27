#' **TODO**
#' - Speed up glmmTMB [why does initial values not make much of a difference? Probably because of overhead on glmmTMB, but hard to change this givne things like mgcv::gam will be slower in the actual fitting process...]
#' 
#' 
#' 
#' Used in assam.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula
    family <- tweedie()
    y = simdat$y
    data = covariate_dat
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    num_archetypes = 5
    num_cores = NULL
    uncertainty_quantification <- TRUE
    control = list(max_iter = 500, tol = 1e-5, temper_prob = 0.8, trace = TRUE)
    bootstrap_control = list(num_boot = 100, ci_alpha = 0.05, seed = NULL, ci_type = "percentile")
    }

#'
#'
#' Used in .quadapprox_fm.R
function() {
    formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula
    family <- nbinom2()
    resp = simdat$y
    data = covariate_dat
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    num_nuisance_perspp = NULL
    start_archetype_labels = NULL
    start = NULL
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


