#' @title Construct predictions from an approximate and scalable SAM object
#' 
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' Takes a fitted \code{assam} object and produces predictions given potentially a new set of observational units with their corresponding covariates. Predictions can be accompanied by uncertainty intervals. 
#' 
#' @param object An object of class \code{assam}.
#' @param newdata A data frame containing the values of the covariates at which predictions are to be calculated, that is, a model matrix from this and \code{object$formula}.
#' @param newoffset A set of offset terms. If supplied, then it must either be a vector (if \code{type = "archetype"}) or a matrix (if \code{type = "species_max"} or \code{"species_mean"}). If the former, the length should be equal to \code{nrow(newdata)}. If the latter, the number of rows should be equal \code{nrow(newdata)} and the number of columns should be equal to \code{nrow(object$spp_effects) i.e., the number of species.}
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
#' The standard point predictor given is constructed based on the asSAM estimates, along with predictions of the species-specific spatial fields if required. If 
#' bootstrapped parameter estimates are also available and \code{se_fit = TRUE}, then two additional point predictors are available, based on the median and mean of the bootstrap predictions. **In our experience, we tend to find that bootstrap median predictor can often be the most stable and accurate, especially if species-specific spatial fields are included in the asSAM**.
#' 
#' Uncertainty intervals produced by \code{predict.assam} are based on the bootstrapped parameter estimates, if available from the \code{object} itself. That is, predictions are constructed using the bootstrapped estimates and then quantiles as appropriate as used to construct the uncertainty intervals. Note particularly with asSAMs including species-specific spatial fields, there is **no general guarantee the standard point predictor necessarily falls inside the corresponding uncertainty intervals constructing using this approach.** This is less of a problem with the bootstrap median/mean predictors though. 
#' 
#' Note archetypal predictions are constructed on the scale of the linear predictions, while species-specific predictions are constructed on the scale of the responses. The latter is analogous to what is available from [fitted.assam()]. 
#' 
#' @return If \code{se_fit = FALSE}, then a matrix of point predictions. If \code{se_fit = TRUE}, then a list with the following components is returned:
#' \item{point_prediction:}{A matrix of predicted values.}
#' \item{bootstrap_median_prediction:}{A matrix of predicted values based on the median of the bootstrap predictions.}
#' \item{bootstrap_mean_prediction:}{A matrix of predicted values based on the mean of the bootstrap predictions.}
#' \item{lower:}{A matrix of the lower limits for the uncertainty intervals.}
#' \item{upper:}{A matrix of the upper limits for the uncertainty intervals.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @examples
#' \dontrun{
#' Please see the help file for assam for example.
#' }
#' 
#' 
#' @export
#' 
#' @importFrom foreach foreach %dopar% %:%
#' @import Matrix
#' @importFrom abind abind
#' @importFrom collapse fquantile
#' @importFrom doParallel registerDoParallel
#' @importFrom sdmTMB sdmTMB
#' @importFrom parallel detectCores
#' @importFrom stats as.formula model.matrix predict median
#' @importFrom TMB MakeADFun sdreport
#' @importFrom utils setTxtProgressBar txtProgressBar
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
    num_spp <- nrow(object$spp_effects)
    num_units <- nrow(newdata)
    
    if(!inherits(object, "assam")) 
        stop("`object' is not of class \"assam\"")
    
    if(is.null(num_cores))
        registerDoParallel(cores = detectCores()-2)
    if(!is.null(num_cores))
        registerDoParallel(cores = num_cores)
    
    if(!is.null(newoffset)) {
        if(type == "archetype") {
            if(!is.vector(newoffset))
                stop("For predictions of archetypal respones on the linear predictor scale, newoffset should be a vector of a same length as nrow(newdata).")
            if(length(newoffset) != nrow(newdata))
                stop("For predictions of archetypal respones on the linear predictor scale, newoffset should be a vector of a same length as nrow(newdata).")
            }
        if(type %in% c("species_max", "species_mean")) {
            if(nrow(newoffset) != nrow(newdata))
                stop("For predictions of species responses, newoffset should have the same number of rows as newdata.")
            if(ncol(newoffset) != num_spp)
                stop("For predictions of species responses, newoffset should have the same number of columns as length(object$spp_intercepts) i.e., the number of species.")
            }
        }

    if(se_fit == TRUE & object$uncertainty_quantification == FALSE)
        stop("Standard errors for prediction can not be calculated since the object$uncertainty_quantification = FALSE.") 
    
    type <- match.arg(type, choices = c("species_max", "species_mean", "archetype"))

        
    ##-----------------------
    #' # Construct X and associated point predictions -- Parallelized across archetypes
    ##-----------------------
    if(type == "archetype") {
        tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse = " ") ) )
        nullfit <- sdmTMB(tmp_formula, 
                          spatial = FALSE,
                          data = data.frame(newdata, response = rnorm(nrow(newdata)))) #' This may not work in the future if smootihng terms are included, say, due to the standardization that needs to be applied    
        X <- model.matrix(nullfit$formula[[1]], data = nullfit$data)
        
        get_eta <- tcrossprod(X[, -object$which_spp_effects, drop = FALSE], object$betas)
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
        dimnames(all_spp_mu) <- list(units = rownames(newdata), spp = rownames(object$spp_effects), archetype = names(object$mixture_proportion))
        
        all_spp_mu <- object$family$linkinv(all_spp_mu)
        
        if(type == "species_max") {
            pt_pred <- sapply(1:num_spp, function(j) all_spp_mu[,j,which.max(object$posterior_probability[j,])])
            }
        
        if(type == "species_mean") {
            pt_pred <- sapply(1:num_spp, function(j) rowSums(matrix(object$posterior_probability[j,], nrow = num_units, ncol = object$num_archetypes, byrow = TRUE) * all_spp_mu[,j,]))
            }
        
        rownames(pt_pred) <- rownames(newdata)
        colnames(pt_pred) <- rownames(object$spp_effects)
        }
    
    
    if(!se_fit)
        return(pt_pred)

        
    ##-----------------------
    #' # Uncertainty quantification -- Parallelized across bootstrap datasets
    ##-----------------------
    if(se_fit) {
        message("Using bootstrap samples to construct uncertainty quantification for predictions...this will take a while so go a brew a cup of tea (or two)!")
        make_pred_tmb_data <- NULL
        if(type %in% c("species_max", "species_mean"))
            make_pred_tmb_data <- lapply(1:num_spp, function(l) predict(object$sdmTMB_fits[[l]], newdata = newdata, return_tmb_object = TRUE)$pred_tmb_data)
        
        construction_predictions_per_bootstrap <- function(k0, newdata, newoffset, pred_tmb_data) {
            if(type == "archetype") {
                cw_bootstrap_parameters <- object$bootstrap_parameters[k0,]
                cw_spp_effects <- matrix(cw_bootstrap_parameters[grep("spp_effects", names(cw_bootstrap_parameters))], nrow = num_spp, byrow = TRUE)
                cw_mixture_proportions <- cw_bootstrap_parameters[grep("mixture_proportion", names(cw_bootstrap_parameters))]
                cw_betas <- cw_bootstrap_parameters[grep("beta[1-9]", names(cw_bootstrap_parameters))]             
                cw_betas <- matrix(cw_betas, nrow = object$num_archetypes, byrow = TRUE)
                
                get_cw_eta <- tcrossprod(X[, -object$which_spp_effects, drop = FALSE], cw_betas)
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
                                                            k0 = k0,
                                                            pred_tmb_data = pred_tmb_data)
                all_spp_mu <- object$family$linkinv(all_spp_mu)
                if(object$num_archetypes == 1)
                    all_spp_mu <- array(all_spp_mu, dim = c(nrow(all_spp_mu), ncol(all_spp_mu), 1))

                
                if(type == "species_max") {
                    pt_pred <- sapply(1:num_spp, function(j) all_spp_mu[,j,which.max(object$bootstrap_posterior_probability[j,,k0])])
                    }
                if(type == "species_mean") {
                    pt_pred <- sapply(1:num_spp, function(j) rowSums(matrix(object$bootstrap_posterior_probability[j,,k0], nrow = num_units, ncol = object$num_archetypes, byrow = TRUE) * all_spp_mu[,j,]))
                    }
                
                rownames(pt_pred) <- rownames(newdata)
                colnames(pt_pred) <- rownames(object$spp_effects)
                }
            
            setTxtProgressBar(pb, k0)
            return(pt_pred)
            }
        
        pb <- txtProgressBar(min = 0, max = nrow(object$bootstrap_parameters), style = 3)
        all_predictions <- foreach(l = 1:nrow(object$bootstrap_parameters)) %do% construction_predictions_per_bootstrap(k0 = l, newdata = newdata, newoffset = newoffset, pred_tmb_data = make_pred_tmb_data)
        all_predictions <- abind::abind(all_predictions, along = 3)
        gc()
        close(pb) 
        rm(pb, make_pred_tmb_data)
        
        ci_alpha <- (1 - coverage)/2
        quantile_predictions <- apply(all_predictions, c(1,2), collapse::fquantile, probs = c(ci_alpha, 1 - ci_alpha))
        median_predictions <- apply(all_predictions, c(1,2), median, na.rm = TRUE)
        mean_predictions <- apply(all_predictions, c(1,2), mean, na.rm = TRUE)
        lower_predictions <- quantile_predictions[1,,]
        upper_predictions <- quantile_predictions[2,,]
        rm(quantile_predictions)
        
        if(type == "archetype") {
            if(object$num_archetypes == 1) {
                lower_predictions <- matrix(lower_predictions, ncol = 1)
                upper_predictions <- matrix(upper_predictions, ncol = 1)
                } 
            rownames(lower_predictions) <- rownames(upper_predictions) <- rownames(median_predictions) <- rownames(mean_predictions) <- rownames(newdata)
            colnames(lower_predictions) <- colnames(upper_predictions) <- colnames(median_predictions) <- colnames(mean_predictions) <- names(object$mixture_proportion)
            }
        if(type %in% c("species_max", "species_mean")) {
            rownames(lower_predictions) <- rownames(upper_predictions) <- rownames(median_predictions) <- rownames(mean_predictions) <- rownames(newdata)
            colnames(lower_predictions) <- colnames(upper_predictions) <- colnames(median_predictions) <- colnames(mean_predictions) <- rownames(object$spp_effects)
            }
        
        return(list(point_prediction = pt_pred, 
                    bootstrap_median_prediction = median_predictions,
                    bootstrap_mean_prediction = mean_predictions,
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
    num_spp <- nrow(object$spp_effects)
    
    out <- NULL    

    for(l1 in 1:num_spp) {
        #' Set up new parameters as per the asSAM
        use_pars <- .get_pars2(object = object$sdmTMB_fits[[l1]])
        cw_b_j <- c(object$spp_effects[l1,], object$betas[k0,])
        cw_b_j[object$which_spp_effects] <- object$spp_effects[l1,]
        cw_b_j[-object$which_spp_effects] <- object$betas[k0,]
        use_pars[["b_j"]] <- cw_b_j
        rm(cw_b_j)
        if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
            use_pars[["ln_phi"]] <- log(object$spp_nuisance$dispersion[l1])
        if(object$family$family[1] == "tweedie") 
            use_pars[["thetaf"]] <- qlogis(object$spp_nuisance$power[l1] - 1)
        if(object$add_spatial) {
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
        if(object$add_spatial) {
            use_map$ln_tau_O <- as.factor(NA)
            use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
            }
        
        #' Construct new TMB object
        new_tmb_obj <- TMB::MakeADFun(data = predict(object$sdmTMB_fits[[l1]], newdata = newdata, return_tmb_object = TRUE)$pred_tmb_data, #make_pred_tmb_data,
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
.predict_eta_sdmTMB_bootstrap <- function(object, newdata, newoffset, k0, pred_tmb_data) {
    num_units <- nrow(newdata)
    num_spp <- nrow(object$spp_effects)
    
    cw_bootstrap_parameters <- object$bootstrap_parameters[k0,]
    cw_spp_effects <- matrix(cw_bootstrap_parameters[grep("spp_effects", names(cw_bootstrap_parameters))], nrow = num_spp, byrow = TRUE)
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
    
    do_fn <- function(l0, l1) {
        #' Set up new parameters as per the asSAM
        use_pars <- .get_pars2(object = object$sdmTMB_fits[[l1]])
        cw_b_j <- c(cw_spp_effects[l1,], cw_betas[l0,])
        cw_b_j[object$which_spp_effects] <- cw_spp_effects[l1,]
        cw_b_j[-object$which_spp_effects] <- cw_betas[l0,]
        use_pars[["b_j"]] <- cw_b_j
        rm(cw_b_j)
        if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
            use_pars[["ln_phi"]] <- cw_ln_phi[l1]
        if(object$family$family[1] == "tweedie") 
            use_pars[["thetaf"]] <- cw_thetaf[l1]
        if(object$add_spatial) {
            use_pars[["ln_kappa"]] <- matrix(cw_ln_kappa[l1], nrow = 2, ncol = 1) # This has two rows as set up as sdmTMB
            use_pars[["ln_tau_O"]] <- cw_ln_tau_O[l1]
            }
        
        #' Truncate spatial parameters that are very large in magnitude 
        if(use_pars[["ln_tau_O"]] < -30) use_pars[["ln_tau_O"]] <- -30
        if(use_pars[["ln_tau_O"]] > 30) use_pars[["ln_tau_O"]] <- 30
        if(any(use_pars[["ln_kappa"]] < -30)) use_pars[["ln_kappa"]] <- matrix(-30, nrow = 2, ncol = 1)
        if(any(use_pars[["ln_kappa"]] > 30)) use_pars[["ln_kappa"]] <- matrix(30, nrow = 2, ncol = 1)
        
        #' Set up new map to constrain all parameters at asSAM estimates
        use_map <- object$sdmTMB_fits[[l1]]$tmb_map
        use_map$b_j <- as.factor(rep(NA, length(use_pars[["b_j"]])))
        if(object$family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie"))
            use_map$ln_phi <- as.factor(NA)
        if(object$family$family[1] == "tweedie") 
            use_map$thetaf <- as.factor(NA)
        if(object$add_spatial) {
            use_map$ln_tau_O <- as.factor(NA)
            use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 
            }
    
        
        #' Construct new TMB object
        new_tmb_obj <- TMB::MakeADFun(data = pred_tmb_data[[l1]],
                                      profile = object$sdmTMB_fits[[l1]]$control$profile,
                                      parameters = use_pars,
                                      map = use_map,
                                      random = object$sdmTMB_fits[[l1]]$tmb_random,
                                      DLL = "sdmTMB",
                                      silent = TRUE)
    
        new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
        new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
        r <- new_tmb_obj$report(new_tmb_obj$env$last.par) # last.par taken since it is the newest set of parameters
        
        #out[,l1,l0] <- r$proj_eta[,1]
        return(r$proj_eta[,1]) 
        }
    
    #' Nested foreach loops used here -- Not sure this saves that much time compared to doing foreach at the bootstrap dataset level but whatever...
    out <- foreach(l2 = 1:num_spp) %:% 
        foreach(l0 = 1:object$num_archetypes, .combine = "cbind") %dopar% {
            do_fn(l0 = l0, l1 = l2)
        } 
    
    out <- abind::abind(out, along = 1.5)
    

    if(!is.null(newoffset)) {
        for(l0 in 1:object$num_archetypes)
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

