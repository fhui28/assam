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
                     mixing_proportion = object$mixing_proportion,
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
                                                                           mixing_proportion = use_model$mixing_proportion,
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
                                                                        mixing_proportion = use_model$mixing_proportion,
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
.fastsimulate_assam <- function(qa_object,
                                em_object,
                                num_X, 
                                which_spp_effects,
                                nsim,
                                do_parallel,
                                num_cores,
                                seed) {
    
    if(do_parallel) {
        if(is.null(num_cores))
            registerDoParallel(cores = detectCores() - 2)
        if(!is.null(num_cores))
            registerDoParallel(cores = num_cores)
        }

    create_seeds <- NULL
    if(!is.null(seed))
        create_seeds <- seed + 1:nsim
    
    
    fastsim_fn <- function(qa_object,
                           em_object,
                           which_spp_effects,
                           num_X,
                           seed) {
        
        num_archetypes <- nrow(em_object$new_betas)
        num_spp <- nrow(qa_object$parameters)
        bootstrap_parameters_dataset <- matrix(NA, nrow = num_spp, ncol = ncol(qa_object$parameters)) 
        
        #' Sample bootstrap archetype labels
        set.seed(seed)
        archetype_label <- sample(1:num_archetypes, size = num_spp, prob = em_object$new_mixprop, replace = TRUE)

        #' Sample bootstrap parameters from normal approximation with mean vector and covariance matrix evaluated at the asSAM estimates
        spp_bootstrap_mean <- cbind(em_object$new_spp_effects, em_object$new_betas[archetype_label,], em_object$new_nuisance)
        spp_bootstrap_mean[, which_spp_effects] <- em_object$new_spp_effects
        spp_bootstrap_mean[, (1:num_X)[-which_spp_effects]] <- em_object$new_betas[archetype_label,]
        
        for(l1 in 1:num_spp) {
            bootstrap_parameters_dataset[l1,] <- mvtnorm::rmvnorm(n = 1, mean = spp_bootstrap_mean[l1,], sigma = solve(qa_object$hessian[[l1]]))
            }
        set.seed(NULL)
        
        return(list(bootstrap_parameters = bootstrap_parameters_dataset, 
                    archetype_label = archetype_label))
        }

    out <- foreach::foreach(l = 1:nsim) %dopar% fastsim_fn(qa_object = qa_object,
                                                           em_object = em_object,
                                                           which_spp_effects = which_spp_effects,
                                                           num_X = num_X,
                                                           seed = create_seeds[l])
    
    return(out)
    }

