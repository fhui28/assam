#' @title Simulate data from approximate and scalable SAM object
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#'
#' Simulate multivariate abundance data based on a fitted \code{assam} object.
#' 
#' @param object An object of class \code{assam}.
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). This has to be supplied since \code{object} may not (necessarily) contain a data argument to leverage.
#' @param nsim A positive integer specifying the number of simulated datasets. Defaults to 1.
#' @param do_parallel Should parallel computing be used to fit the asSAM. Defaults to \code{TRUE}, and should be kept this way as much as possible as parallel computing is one of the key ingredients in making asSAMs scalable. 
#' @param num_cores If \code{do_parallel = TRUE}, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to \code{parallel::detectCores() - 2}.
#' @param seed An integer to set seed number. Defaults to a random seed number.
#' @param ... Not used.
#' 
#' @details 
#' Simulates multivariate abundance data from a species archetype model (SAM). For the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_k,}
#' 
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_k} denotes the corresponding regression coefficients for archetype \eqn{k}. Based on the mean model given above, responses \eqn{y_{ij}} are then simulated from the assumed distribution, using the additional dispersion and power parameters as appropriate.
#' 
#' @return A matrix and \code{nsim} columns, where each column contains the output from one run of [create_samlife()].
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' @examples
#' \dontrun{
#' #' Please see the help file for assam for example.
#' }
#' 
#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats as.formula plogis
#' @md


simulate.assam <- function(object, 
                           data,
                           nsim = 1, 
                           do_parallel = TRUE, 
                           num_cores = NULL, 
                           seed = NULL, 
                           ...) {
    
    if(!inherits(object, "assam")) 
        stop("object must be of class assam.")
    if(do_parallel) {
        if(is.null(num_cores))
            registerDoParallel(cores = detectCores() - 2)
        if(!is.null(num_cores))
            registerDoParallel(cores = num_cores)
        }
    
    use_model <- list(family = object$family, 
                     formula = object$formula,
                     data = data,
                     which_spp_effects = object$which_spp_effects,
                     offset = object$offset,
                     betas = object$betas, 
                     spp_effects = object$spp_effects, 
                     mixture_proportion = object$mixture_proportion, 
                     mesh = object$mesh,
                     spp_spatial_sd = object$spp_nuisance$spatial_SD,
                     spp_spatial_range = object$spp_nuisance$spatial_range,
                     trial_size = object$trial_size, 
                     spp_dispparam = NULL, 
                     spp_powerparam = NULL) 
    if(object$family$family[1] %in% c("Beta","nbinom2","Gamma","gaussian","tweedie"))
        use_model$spp_dispparam <- object$spp_nuisance$dispersion
    if(object$family$family[1] %in% c("tweedie"))
        use_model$spp_powerparam <- object$spp_nuisance$power
    if(object$add_spatial) {
        use_model$spp_spatial_sd[!is.finite(use_model$spp_spatial_sd)] <- 1e-6 # This is arbitrary!!!
        use_model$spp_spatial_sd[use_model$spp_spatial_sd < 1e-6] <- 1e-6
        use_model$spp_spatial_sd[use_model$spp_spatial_sd > 1e6] <- 1e6
        use_model$spp_spatial_range[!is.finite(use_model$spp_spatial_range)] <- 1e6 # This is arbitrary!!!
        use_model$spp_spatial_range[use_model$spp_spatial_range < 1e-6] <- 1e-6
        use_model$spp_spatial_range[use_model$spp_spatial_range > 1e6] <- 1e6
        }
    
    create_seeds <- NULL
    if(!is.null(seed))
        create_seeds <- seed + 1:nsim

    if(do_parallel)
        out <- foreach::foreach(l = 1:nsim,
                                .combine = "cbind") %dopar% create_samlife(family = use_model$family, 
                                                                           formula = use_model$formula,
                                                                           data = use_model$data,
                                                                           which_spp_effects = use_model$which_spp_effects,
                                                                           betas = use_model$betas, 
                                                                           offset = use_model$offset,
                                                                           spp_effects = use_model$spp_effects, 
                                                                           spp_dispparam = use_model$spp_dispparam, 
                                                                           spp_powerparam = use_model$spp_powerparam, 
                                                                           mixture_proportion = use_model$mixture_proportion, 
                                                                           mesh = use_model$mesh,
                                                                           spp_spatial_sd = use_model$spp_spatial_sd,
                                                                           spp_spatial_range = use_model$spp_spatial_range,
                                                                           trial_size = use_model$trial_size,
                                                                           seed = create_seeds[l])
    if(!do_parallel)
        out <- foreach::foreach(l = 1:nsim,
                                .combine = "cbind") %do% create_samlife(family = use_model$family, 
                                                                        formula = use_model$formula,
                                                                        data = use_model$data,
                                                                        which_spp_effects = use_model$which_spp_effects,
                                                                        betas = use_model$betas, 
                                                                        offset = use_model$offset,
                                                                        spp_effects = use_model$spp_effects, 
                                                                        spp_dispparam = use_model$spp_dispparam, 
                                                                        spp_powerparam = use_model$spp_powerparam, 
                                                                        mixture_proportion = use_model$mixture_proportion, 
                                                                        mesh = use_model$mesh,
                                                                        spp_spatial_sd = use_model$spp_spatial_sd,
                                                                        spp_spatial_range = use_model$spp_spatial_range,
                                                                        trial_size = use_model$trial_size,
                                                                        seed = create_seeds[l])
    
    #doParallel::stopImplicitCluster()
    return(out)
    }   


#' @rdname simulate.assam
#' @export 
simulate <- function(object, ...) {
        UseMethod("simulate")
    }


#' @noRd
#' @noMd
.fastsimulate_assam <- function(object,
                                nsim,
                                do_parallel,
                                num_cores,
                                seed) {
    
    if(!inherits(object, "assam")) 
        stop("object must be of class assam.")
    if(do_parallel) {
        if(is.null(num_cores))
            registerDoParallel(cores = detectCores() - 2)
        if(!is.null(num_cores))
            registerDoParallel(cores = num_cores)
        }

    use_model <- list(family = object$family, 
                      formula = object$formula,
                      which_spp_effects = object$which_spp_effects,
                      offset = object$offset,
                      betas = object$betas, 
                      spp_effects = object$spp_effects, 
                      mixture_proportion = object$mixture_proportion, 
                      mesh = object$mesh,
                      spp_spatial_sd = object$spp_nuisance$spatial_SD,
                      spp_spatial_range = object$spp_nuisance$spatial_range,
                      trial_size = object$trial_size, 
                      spp_dispparam = NULL, 
                      spp_powerparam = NULL) 
    if(object$family$family[1] %in% c("Beta","nbinom2","Gamma","gaussian","tweedie"))
        use_model$spp_dispparam <- object$spp_nuisance$dispersion
    if(object$family$family[1] %in% c("tweedie"))
        use_model$spp_powerparam <- object$spp_nuisance$power
    if(object$add_spatial) {
        use_model$spp_spatial_sd[!is.finite(use_model$spp_spatial_sd)] <- 1e-6 # This is arbitrary!!!
        use_model$spp_spatial_sd[use_model$spp_spatial_sd < 1e-6] <- 1e-6
        use_model$spp_spatial_sd[use_model$spp_spatial_sd > 1e6] <- 1e6
        use_model$spp_spatial_range[!is.finite(use_model$spp_spatial_range)] <- 1e6 # This is arbitrary!!!
        use_model$spp_spatial_range[use_model$spp_spatial_range < 1e-6] <- 1e-6
        use_model$spp_spatial_range[use_model$spp_spatial_range > 1e6] <- 1e6
        }
    
    create_seeds <- NULL
    if(!is.null(seed))
        create_seeds <- seed + 1:nsim
    
    
    fastsim_fn <- function(family,
                           mixture_proportion,
                           which_spp_effects,
                           betas,
                           spp_effects,
                           spp_dispparam,
                           spp_powerparam,
                           add_spatial,
                           spp_spatial_sd,
                           spp_spatial_range,
                           seed) {
        
        num_archetypes <- nrow(betas)
        num_spp <- nrow(spp_effects)
        bootstrap_parameters_dataset <- matrix(NA, nrow = num_spp, ncol = ncol(get_qa$parameters)) 
        
        #' Sample bootstrap archetype labels
        set.seed(seed)
        archetype_label <- sample(1:num_archetypes, size = num_spp, prob = mixture_proportion, replace = TRUE)
        set.seed(NULL)
        
        #' Sample bootstrap parameters from normal approximation with mean vector and covariance matrix evaluated at the asSAM estimates
        for(l1 in 1:num_spp) {
            cw_seed <- NULL
            if(!is.null(seed))
                cw_seed <- seed + l1
        #     
        #     
        #     if(!add_spatial) {
        #         spp_bootstrap_mean <- object$sdmTMB_fits[[l1]]$parameters
        #         cw_b_j <- c(spp_effects[l1,], betas[archetype_label[l1],])
        #         cw_b_j[which_spp_effects] <- spp_effects[l1,]
        #         cw_b_j[-which_spp_effects] <- betas[archetype_label[l1],]
        #         spp_bootstrap_mean[grep("b_j", names(spp_bootstrap_mean))] <- cw_b_j
        #         rm(cw_b_j)
        #         if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
        #             use_pars[grep("ln_phi", names(spp_bootstrap_mean))] <- log(spp_dispparam[l1])
        #         if(family$family[1] == "tweedie") 
        #             use_pars[grep("thetaf", names(spp_bootstrap_mean))] <- qlogis(spp_powerparam[l1] - 1)
        #         
        #         spp_bootstrap_covariance <- solve(object$sdmTMB_fits[[l1]]$tmb_obj$he(spp_bootstrap_mean)) #' Not guaranteed to be positive definite!!!
        #         }
                
            
            #' #' Set up TMB object evaluated at asSAM estimates
            #' use_pars <- .get_pars2(object = object$sdmTMB_fits[[l1]])
            #' cw_b_j <- c(spp_effects[l1,], betas[archetype_label[l1],])
            #' cw_b_j[which_spp_effects] <- spp_effects[l1,]
            #' cw_b_j[-which_spp_effects] <- betas[archetype_label[l1],]
            #' use_pars[["b_j"]] <- as.vector(unlist(cw_b_j))
            #' rm(cw_b_j)
            #' 
            #' if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
            #'     use_pars[["ln_phi"]] <- log(spp_dispparam[l1])
            #' if(family$family[1] == "tweedie") 
            #'     use_pars[["thetaf"]] <- qlogis(spp_powerparam[l1] - 1)
            #' if(add_spatial) {
            #'     use_pars[["ln_kappa"]] <- matrix(log(1/spp_spatial_range[l1]), nrow = 2, ncol = 1) # This has two rows as set up as sdmTMB
            #'     use_pars[["ln_tau_O"]] <- 1/(sqrt(4*pi) * object$spp_spatial_sd[l1] * exp(use_pars[["ln_kappa"]][1,1])) 
            #'     }
            #' 
            #' #' Truncate spatial parameters that are very large in magnitude 
            #' if(use_pars[["ln_tau_O"]] < -30) use_pars[["ln_tau_O"]] <- -30
            #' if(use_pars[["ln_tau_O"]] > 30) use_pars[["ln_tau_O"]] <- 30
            #' if(any(use_pars[["ln_kappa"]] < -30)) use_pars[["ln_kappa"]] <- matrix(-30, nrow = 2, ncol = 1)
            #' if(any(use_pars[["ln_kappa"]] > 30)) use_pars[["ln_kappa"]] <- matrix(30, nrow = 2, ncol = 1)
            #' 
            #' #' Set up new map to constrain all parameters at asSAM estimates
            #' use_map <- object$sdmTMB_fits[[l1]]$tmb_map
            #' use_map$b_j <- as.factor(rep(NA, length(use_pars[["b_j"]])))
            #' if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie"))
            #'     use_map$ln_phi <- as.factor(NA)
            #' if(family$family[1] == "tweedie") 
            #'     use_map$thetaf <- as.factor(NA)
            #' if(add_spatial) {
            #'     use_map$ln_tau_O <- as.factor(NA)
            #'     use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
            #'     }
            #' 
            #'     
            #' #' Construct new TMB object
            #' new_tmb_obj <- TMB::MakeADFun(data = object$sdmTMB_fits[[l1]]$tmb_data, #make_pred_tmb_data,
            #'                               profile = object$sdmTMB_fits[[l1]]$control$profile,
            #'                               parameters = use_pars,
            #'                               map = use_map,
            #'                               random = object$sdmTMB_fits[[l1]]$tmb_random,
            #'                               DLL = "sdmTMB",
            #'                               silent = TRUE)
            #' 
            #' new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
            #' new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
            #' 
            #' spp_bootstrap_mean <- new_tmb_sdreport$par.fixed
            #' spp_bootstrap_covariance <- new_tmb_sdreport$cov.fixed
            
            set.seed(cw_seed)
            bootstrap_parameters_dataset[l1,] <- mvtnorm::rmvnorm(n = 1, mean = spp_bootstrap_mean, sigma = spp_bootstrap_covariance)
            set.seed(NULL)
            }
            
        return(bootstrap_parameters_dataset)
        }
        

    out <- foreach::foreach(l = 1:nsim, 
                            .combine = "cbind") %do% fastsim_fn(family = use_model$family, 
                                                                mixture_proportion = use_model$mixture_proportion, 
                                                                which_spp_effects = use_model$which_spp_effects,
                                                                betas = use_model$betas, 
                                                                spp_effects = use_model$spp_effects, 
                                                                spp_dispparam = use_model$spp_dispparam, 
                                                                spp_powerparam = use_model$spp_powerparam, 
                                                                add_spatial = object$add_spatial,
                                                                spp_spatial_sd = use_model$spp_spatial_sd,
                                                                spp_spatial_range = use_model$spp_spatial_range,
                                                                seed = create_seeds[l])
    
    
    }

