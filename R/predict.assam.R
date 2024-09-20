#' @title Construct predictions from an approximate and scalable SAM object
#' 
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' Takes a fitted \code{assam} object and produces predictions given (potentially) a new set of observational units with their corresponding covariates. Predictions can be accompanied by uncertainty intervals. 
#' 
#' @param object An object of class \code{assam}.
#' @param newdata A data frame containing the values of the covariates at which predictions are to be calculated, that is, a model matrix from this and \code{object$formula}.
#' @param newoffset A set of offset terms. If supplied, then it must either be a vector (if \code{type = "archetype"}) or a matrix (if \code{type = "species_max"} or \code{"species_mean"}). If the former, the length should be equal to \code{nrow(newdata)}. If the latter, the number of rows should be equal \code{nrow(newdata)} and the number of columns should be equal to \code{length(object$spp_intercepts) i.e., the number of species.}
#' @param type The type of prediction required:  
#' If \code{type = "archetype"}, then archetypal predictions are constructed on the scale of the linear predictors i.e., \eqn{x_i^\top\beta_k};
#' If \code{type = "species_max"}, then species-specific predictions on the scale of the responses i.e., \eqn{\mu_{ij}}, are constructed based on the most likely archetype the species belongs to, as judged by the posterior probabilities.
#' If \code{type = "species_mean"}, then species-specific predictions on the scale of the responses i.e., \eqn{\mu_{ij}}, are constructed based on a weighted mean of the predicted responses from each archetype, where the weights are given by the posterior probabilities.
#' @param se_fit When this is set to \code{TRUE} (not default), then uncertainty intervals are returned for the point predictions.
#' @param num_cores To speed up calculation of the uncertainty intervals, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-2}.
#' @param coverage The coverage probability of the uncertainty intervals for prediction. Defaults to 0.95, which corresponds to 95% uncertainty intervals.
#' @param ... Not used.
#' 
#' @details 
#' Uncertainty intervals produced by \code{predict.assam} are based on the bootstrapped parameter estimates, if available from the \code{object} itself. That is, predictions are constructed using the bootstrapped estimates and then quantiles as appropriate as used to construct the uncertainty intervals. **Note particularly with asSAMs including species-specific spatial fields, there is no general guarantee that the point prediction necessarily falls inside the corresponding uncertainty intervals constructing using this approach.**
#' 
#' Note archetypal predictions are constructed on the scale of the linear predictions, while species-specific predictions are constructed on the scale of the responses. The latter is analogous to what is available from [fitted.assam()]. 
#' 
#' @return If \code{se_fit = FALSE}, then a matrix of point predictions. If \code{se_fit = TRUE}, then a list with the following components is returned:
#' \item{point_prediction:}{A matrix of predicted values.}
#' \item{lower:}{A matrix of the lower limits for the uncertainty intervals.}
#' \item{upper:}{A matrix of the upper limits for the uncertainty intervals.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @examples
#' \dontrun{
#' ##----------------------
#' # Generate some multivariate abundance data from a SAM
#' ##----------------------
#' library(tidyverse)
#' library(mvtnorm)
#' library(GGally)
#' 
#' set.seed(092024)
#' 
#' num_X <- 10
#' num_units <- 1000
#' num_spp <- 80
#' num_archetype <- 5
#' H <- outer(1:num_X, 1:num_X, "-")
#' H <- 0.5^abs(H)
#' covariate_dat <- rmvnorm(num_units, sigma = H) %>% 
#'     as.data.frame %>% 
#'     rename_with(., .fn = function(x) paste0("covariate", x))
#' rm(H)
#' 
#' true_betas <- runif(num_archetype * num_X, -1, 1) %>% matrix(nrow = num_archetype)
#' true_intercepts <- runif(num_spp, -3, 0)  
#' true_dispparam <- 1/runif(num_spp, 0, 5) 
#' true_powerparam <- runif(num_spp, 1.4, 1.8)
#' true_mixprop <- c(0.2, 0.25, 0.3, 0.1, 0.15)
#'  
#' simdat <- create_samlife(family = nbinom2(), 
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula, 
#' data = covariate_dat, 
#' betas = true_betas, 
#' spp_intercept = true_intercepts, 
#' spp_dispparam = true_dispparam, 
#' spp_powerparam = true_powerparam, 
#' mixture_proportion = true_mixprop,
#' seed = 092024)
#'  
#'  
#' ##----------------------
#' # Fit asSAM and assess results 
#' #' **Most users should start here**
#' ##----------------------
#' tic <- proc.time()
#' samfit <- assam(y = simdat$y,
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' family = nbinom2(),
#' num_archetypes = num_archetype,
#' num_cores = 8)
#' toc <- proc.time
#'  
#' #' Archetype-level predictions
#' predict(samfit, newdata = covariate_dat, type = "archetype", se_fit = TRUE) 
#' 
#' #' Species-level predictions
#' predict(samfit,  newdata = covariate_dat, type = "species_max", se_fit = TRUE) 
#'  
#' }
#' 
#' 
#' @export
#' 
#' @importFrom foreach foreach %dopar%
#' @import Matrix
#' @importFrom abind abind
#' @importFrom collapse fquantile
#' @importFrom doParallel registerDoParallel
#' @importFrom sdmTMB sdmTMB
#' @importFrom parallel detectCores
#' @importFrom stats as.formula model.matrix predict
#' @importFrom TMB MakeADFun sdreport
#' @md

    
predict.assam <- function(object, 
                          newdata, 
                          newoffset = NULL,
                          type = c("species_max", "species_mean", "archetype"), 
                          se_fit = FALSE, 
                          num_cores = NULL,
                          coverage = 0.95, 
                          ...) {
    
    ##-----------------------
    #' # Checks and balances
    ##-----------------------
    num_spp <- length(object$spp_intercepts)
    num_units <- nrow(newdata)
    
    if(!inherits(object, "assam")) 
        stop("`object' is not of class \"assam\"")
    
    if(is.null(num_cores))
        registerDoParallel(cores = detectCores()-2)
    if(!is.null(num_cores))
        registerDoParallel(cores = num_cores)
    
    if(!is.null(newoffset)) {
        if(type == "archetype" & !as.vector(newoffset))
            stop("For predictions of archetypal respones on the linear predictor scale, newoffset should be a vector of a same length as nrow(newdata).")
        if(type == "archetype" & length(newoffset) != nrow(newdata))
            stop("For predictions of archetypal respones on the linear predictor scale, newoffset should be a vector of a same length as nrow(newdata).")
        if(type %in% c("species_max", "species_mean") & nrow(newoffset) != nrow(newdata))
            stop("For predictions of species responses, newoffset should have the same number of rows as newdata.")
        if(type %in% c("species_max", "species_mean") & ncol(newoffset) != num_spp)
            stop("For predictions of species responses, newoffset should have the same number of columns as length(object$spp_intercepts) i.e., the number of species.")
        }

    if(se_fit == TRUE & object$uncertainty_quantification == FALSE)
        stop("Standard errors for prediction can not be calculated since the object$uncertainty_quantification = FALSE.") 
    
    type <- match.arg(type, choices = c("species_max", "species_mean", "archetype"))

        
    ##-----------------------
    #' # Construct X and associated point predictions
    ##-----------------------
    if(type == "archetype") {
        tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse = " ") ) )
        nullfit <- sdmTMB(tmp_formula, 
                          spatial = FALSE,
                          data = data.frame(newdata, response = rnorm(nrow(newdata)))) #' This may not work in the future if smootihng terms are included, say, due to the standardization that needs to be applied    
        X <- model.matrix(nullfit$formula[[1]], data = nullfit$data)[,-1] # Remove the intercept term
        
        get_eta <- tcrossprod(X, object$betas)
        if(!is.null(newoffset))
            get_eta <- get_eta + newoffset
        
        pt_pred <- get_eta 
        rownames(pt_pred) <- rownames(newdata)
        colnames(pt_pred) <- names(object$mixture_proportion)
        rm(tmp_formula, nullfit)
        }
    
    if(type %in% c("species_max", "species_mean")) {
        pred_per_archetype <- function(k0) {
            out <- .predict_eta_sdmTMB(object = object, newdata = newdata, k0 = k0)
            if(!is.null(newoffset))
                out <- out + newoffset
            
            return(out)
            }
        all_spp_mu <- foreach(l = 1:object$num_archetypes) %dopar% pred_per_archetype(k0 = l) 
        rm(pred_per_archetype)
        all_spp_mu <- abind::abind(all_spp_mu, along = 3)
        dimnames(all_spp_mu) <- list(units = rownames(newdata), spp = names(object$spp_intercepts), archetype = names(object$mixture_proportion))
        
        all_spp_mu <- object$family$linkinv(all_spp_mu)
        
        if(type == "species_max") {
            pt_pred <- sapply(1:num_spp, function(j) all_spp_mu[,j,which.max(object$posterior_probability[j,])])
            }
        
        if(type == "species_mean") {
            pt_pred <- sapply(1:num_spp, function(j) rowSums(matrix(object$posterior_probability[j,], nrow = num_units, ncol = object$num_archetypes, byrow = TRUE) * all_spp_mu[,j,]))
            }
        
        rownames(pt_pred) <- rownames(newdata)
        colnames(pt_pred) <- names(object$spp_intercepts)
        }
    
    
    if(!se_fit)
        return(pt_pred)

        
    ##-----------------------
    #' # Uncertainty quantification
    ##-----------------------
    if(se_fit) {
        message("Using bootstrap samples to construct uncertainty quantification for predictions...this will take a while so go a brew a cup of tea (or two)!")
        
        construction_predictions_per_bootstrap <- function(k0, newdata, newoffset) {
            if(type == "archetype") {
                cw_bootstrap_parameters <- object$bootstrap_parameters[k0,]
                cw_spp_intercepts <- cw_bootstrap_parameters[grep("spp_intercept", names(cw_bootstrap_parameters))]
                cw_mixture_proportions <- cw_bootstrap_parameters[grep("mixture_proportion", names(cw_bootstrap_parameters))]
                cw_betas <- cw_bootstrap_parameters[grep("beta[1-9]", names(cw_bootstrap_parameters))]             
                cw_betas <- matrix(cw_betas, nrow = object$num_archetypes, byrow = TRUE)
                
                get_cw_eta <- tcrossprod(X, cw_betas)
                if(!is.null(newoffset))
                    get_eta <- get_eta + newoffset
                
                pt_pred <- get_cw_eta 
                rownames(pt_pred) <- rownames(newdata)
                colnames(pt_pred) <- names(object$mixture_proportion)
                }
            
            if(type %in% c("species_max", "species_mean")) {
                all_spp_mu <- .predict_eta_sdmTMB_bootstrap(object = object, 
                                                            newdata = newdata, 
                                                            newoffset = newoffset, 
                                                            k0 = k0)
                all_spp_mu <- object$family$linkinv(all_spp_mu)
                
                if(type == "species_max") {
                    pt_pred <- sapply(1:num_spp, function(j) all_spp_mu[,j,which.max(object$bootstrap_posterior_probability[j,,k0])])
                    }
                
                if(type == "species_mean") {
                    pt_pred <- sapply(1:num_spp, function(j) rowSums(matrix(object$bootstrap_posterior_probability[j,,k0], nrow = num_units, ncol = object$num_archetypes, byrow = TRUE) * all_spp_mu[,j,]))
                    }
                
                rownames(pt_pred) <- rownames(newdata)
                colnames(pt_pred) <- names(object$spp_intercepts)
                }
            
            return(pt_pred)
            }
        
        all_predictions <- foreach(l = 1:nrow(object$bootstrap_parameters)) %dopar% construction_predictions_per_bootstrap(k0 = l, newdata = newdata, newoffset = newoffset)
        all_predictions <- abind::abind(all_predictions, along = 3)
        gc()
        
        ci_alpha <- (1 - coverage)/2
        quantile_predictions <- apply(all_predictions, c(1,2), collapse::fquantile, probs = c(ci_alpha, 1 - ci_alpha))
        lower_predictions <- quantile_predictions[1,,]
        upper_predictions <- quantile_predictions[2,,]
        rm(quantile_predictions)
        
        if(type == "archetype") {
            rownames(lower_predictions) <- rownames(upper_predictions) <- rownames(newdata)
            colnames(lower_predictions) <- colnames(upper_predictions) <- names(object$mixture_proportion)
            }
        if(type %in% c("species_max", "species_mean")) {
            rownames(lower_predictions) <- rownames(upper_predictions) <- rownames(newdata)
            colnames(lower_predictions) <- colnames(upper_predictions) <- names(object$spp_intercepts)
            }
        
        return(list(point_prediction = pt_pred, 
                    lower = lower_predictions,
                    upper = upper_predictions))
        }
    }



# Hidden function to construct all species-specific predictions for the k0-th archetype, potentially to new data
# The sdmTMB object is extracted and then manipulated into a new TMB object fixed at the asSAM parameter estimates. This is then used as the basis for prediction
#' @noMd
#' @noRd
.predict_eta_sdmTMB <- function(object, newdata, k0) {
    num_units <- nrow(newdata)
    num_spp <- length(object$spp_intercepts)
    
    out <- NULL    
    make_pred_tmb_data <- predict(object$sdmTMB_fits[[1]], newdata = newdata, return_tmb_object = TRUE)$pred_tmb_data ## Assume this is the same across all species...do not see any reason why it should not be?!
    
    for(l1 in 1:length(object$spp_intercepts)) {
        #' Set up new parameters as per the asSAM
        use_pars <- .get_pars2(object = object$sdmTMB_fits[[l1]])
        use_pars[["b_j"]] <- c(object$spp_intercepts[l1], object$betas[k0,])
        if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
            use_pars[["ln_phi"]] <- log(object$spp_nuisance$dispersion[l1])
        if(object$family$family[1] == "tweedie") 
            use_pars[["thetaf"]] <- qlogis(object$spp_nuisance$power[l1] - 1)
        if(!is.null(object$mesh)) {
            use_pars[["ln_kappa"]] <- matrix(log(1/object$spp_nuisance$spatial_range[l1]), nrow = 2, ncol = 1) # This has two rows as set up as sdmTMB
            use_pars[["ln_tau_O"]] <- 1/(sqrt(4*pi) * object$spp_nuisance$spatial_SD[l1] * exp(use_pars[["ln_kappa"]][1,1])) 
            }
        
        #' Truncate spatial parameters that are very large in magnitude 
        if(use_pars[["ln_tau_O"]] < -30) use_pars[["ln_tau_O"]] <- -30
        if(use_pars[["ln_tau_O"]] > 30) use_pars[["ln_tau_O"]] <- 30
        if(any(use_pars[["ln_kappa"]] < -30)) use_pars[["ln_kappa"]] <- matrix(-30, nrow = 2, ncol = 1)
        if(any(use_pars[["ln_kappa"]] > 30)) use_pars[["ln_kappa"]] <- matrix(30, nrow = 2, ncol = 1)
        
        #' Set up new MAP to constrain all parameters at asSAM estimates
        use_map <- object$sdmTMB_fits[[l1]]$tmb_map
        use_map$b_j <- as.factor(rep(NA, length(use_pars[["b_j"]])))
        if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie"))
            use_map$ln_phi <- as.factor(NA)
        if(object$family$family[1] == "tweedie") 
            use_map$thetaf <- as.factor(NA)
        if(!is.null(object$mesh)) {
            use_map$ln_tau_O <- as.factor(NA)
            use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
            }
        
        #' Construct new TMB object
        new_tmb_obj <- TMB::MakeADFun(data = make_pred_tmb_data,
                                      profile = object$sdmTMB_fits[[l1]]$control$profile,
                                      parameters = use_pars,
                                      map = use_map,
                                      random = object$sdmTMB_fits[[l1]]$tmb_random,
                                      DLL = "sdmTMB",
                                      silent = TRUE)
        
        new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
        new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
        r <- new_tmb_obj$report(new_tmb_obj$env$last.par) # last.par taken since it is the newest set of parameters
        
        out <- cbind(out, r$proj_eta[,1])
        }

    return(out)
    }
    


# Hidden function to construct all species-specific predictions, across all archetypes, for the k0-th bootstrap dataset
# The sdmTMB object is extracted and then manipulated into a new TMB object fixed at the asSAM parameter estimates. This is then used as the basis for prediction
#' @noMd
#' @noRd
.predict_eta_sdmTMB_bootstrap <- function(object, newdata, newoffset, k0) {
    num_units <- nrow(newdata)
    num_spp <- length(object$spp_intercepts)
    
    cw_bootstrap_parameters <- object$bootstrap_parameters[k0,]
    cw_spp_intercepts <- cw_bootstrap_parameters[grep("spp_intercept", names(cw_bootstrap_parameters))]
    cw_mixture_proportions <- cw_bootstrap_parameters[grep("mixture_proportion", names(cw_bootstrap_parameters))]
    cw_betas <- cw_bootstrap_parameters[grep("beta[1-9]", names(cw_bootstrap_parameters))]             
    cw_betas <- matrix(cw_betas, nrow = object$num_archetypes, byrow = TRUE)
    if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
        cw_ln_phi <- cw_bootstrap_parameters[grep("ln_phi", names(cw_bootstrap_parameters))]
    if(object$family$family[1] == "tweedie") 
        cw_thetaf <- cw_bootstrap_parameters[grep("thetaf", names(cw_bootstrap_parameters))]
    if(!is.null(object$mesh)) {
        cw_ln_kappa <- cw_bootstrap_parameters[grep("ln_kappa", names(cw_bootstrap_parameters))]
        cw_ln_tau_O <- cw_bootstrap_parameters[grep("ln_tau_O", names(cw_bootstrap_parameters))] 
        }
    
    
    out <- array(NA, dim = c(num_units, num_spp, object$num_archetypes))
    make_tmb_pred_data <-predict(object$sdmTMB_fits[[1]], newdata = newdata, return_tmb_object = TRUE)$pred_tmb_data
    
    for(l0 in 1:object$num_archetypes) {
        for(l1 in 1:length(object$spp_intercepts)) {
            #' Set up new parameters as per the asSAM
            use_pars <- .get_pars2(object = object$sdmTMB_fits[[l1]])
            use_pars[["b_j"]] <- c(cw_spp_intercepts[l1], cw_betas[l0,])
            if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
                use_pars[["ln_phi"]] <- cw_ln_phi[l1]
            if(object$family$family[1] == "tweedie") 
                use_pars[["thetaf"]] <- cw_thetaf[l1]
            if(!is.null(object$mesh)) {
                use_pars[["ln_kappa"]] <- matrix(cw_ln_kappa[l1], nrow = 2, ncol = 1) # This has two rows as set up as sdmTMB
                use_pars[["ln_tau_O"]] <- cw_ln_tau_O[l1]
                }
            
            #' Truncate spatial parameters that are very large in magnitude 
            if(use_pars[["ln_tau_O"]] < -30) use_pars[["ln_tau_O"]] <- -30
            if(use_pars[["ln_tau_O"]] > 30) use_pars[["ln_tau_O"]] <- 30
            if(any(use_pars[["ln_kappa"]] < -30)) use_pars[["ln_kappa"]] <- matrix(-30, nrow = 2, ncol = 1)
            if(any(use_pars[["ln_kappa"]] > 30)) use_pars[["ln_kappa"]] <- matrix(30, nrow = 2, ncol = 1)
            
            #' Set up new MAP to constrain all parameters at asSAM estimates
            use_map <- object$sdmTMB_fits[[l1]]$tmb_map
            use_map$b_j <- as.factor(rep(NA, length(use_pars[["b_j"]])))
            if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie"))
                use_map$ln_phi <- as.factor(NA)
            if(object$family$family[1] == "tweedie") 
                use_map$thetaf <- as.factor(NA)
            if(!is.null(object$mesh)) {
                use_map$ln_tau_O <- as.factor(NA)
                use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
                }
        
            
            #' Construct new TMB object
            new_tmb_obj <- TMB::MakeADFun(data = make_tmb_pred_data,
                                          profile = object$sdmTMB_fits[[l1]]$control$profile,
                                          parameters = use_pars,
                                          map = use_map,
                                          random = object$sdmTMB_fits[[l1]]$tmb_random,
                                          DLL = "sdmTMB",
                                          silent = TRUE)
        
            new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
            new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
            r <- new_tmb_obj$report(new_tmb_obj$env$last.par) # last.par taken since it is the newest set of parameters
            
            out[,l1,l0] <- r$proj_eta[,1]
            } 
        
        if(!is.null(newoffset))
            out[,,l0] <- out[,,l0] + newoffset
        }
    
    return(out)
    }



# Hidden function taken directly from sdmTMB utils.R
#' @noMd
#' @noRd
.get_pars2 <- function(object) {
    # based on glmmTMB:
    ee <- object$tmb_obj$env
    x <- ee$last.par.best
    if (length(ee$random) > 0) x <- x[-ee$random]
    p <- ee$parList(x = x)
    p
}

