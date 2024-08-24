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
#' If \code{type = "species_max"}, then species-specific predictions on the scale of the responses i.e., \eqn{\mu_{ij}}, are constructed based on the most likely archetype the species belongs to (as based on the posterior probabilities);
#' If \code{type = "species_mean"}, then species-specific predictions on the scale of the responses i.e., \eqn{\mu_{ij}}, are constructed based on a weighted mean of the predictions from each archetype, where the weights are given by the posterior probabilities.
#' 
#' Note for both types of species-specific predictions, **currently uncertainty of the posterior probabilities are not taken into account...this is on the to-do list!**
#' @param se_fit When this is set to \code{TRUE} (not default), then uncertainty intervals are returned for the point predictions.
#' @param num_cores To speed up calculation of the uncertainty intervals, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-2}.
#' @param coverage The coverage probability of the uncertainty intervals for prediction. Defaults to 0.95, which corresponds to 95% uncertainty intervals.
#' 
#' @details 
#' Uncertainty intervals produced by \code{predict.assam} are based on the bootstrapped parameter estimates, if available from the \code{object} itself. That is, predictions are constructed using the bootstrapped estimates and then quantiles as appropriate as used to construct the uncertainty intervals.
#' 
#' Note archetypal predictions are constructed on the scale of the linear predictions, while species-specific predictions are constructed on the scale of the respones. The latter is analogous to what is available from [fitted.assam()]. 
#' 
#' @return If \code{se_fit = FALSE}, then a matrix of point predictions. If \code{se_fit = TRUE}, then a list with the following components is returned:
#' \item{point_prediction:}{A matrix of predicted values.}
#' \item{lower:}{A matrix of the lower limits for the uncertainty intervals.}
#' \item{upper:}{A matrix of the upper limits for the uncertainty intervals.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @export
#' 
#' @importFrom foreach foreach %dopar%
#' @import Matrix
#' @importFrom abind abind
#' @importFrom doParallel registerDoParallel
#' @importFrom glmmTMB glmmTMB model.matrix
#' @importFrom parallel detectCores
#' @importFrom collapse fquantile
#' @md


function() {
     object <- testfit
     newdata = covariate_dat
     newoffset = NULL
     type = "archetype"
     se_fit = TRUE
     num_cores = 8
     coverage = 0.95
     }
    
predict.assam <- function(object, 
                          newdata, 
                          newoffset = NULL,
                          type = c("species_max", "species_mean", "archetype"), 
                          se_fit = FALSE, 
                          num_cores = NULL,
                          coverage = 0.95) {
    
    ##-----------------------
    #' # Checks and balances
    ##-----------------------
    num_spp <- length(object$spp_intercepts)
    
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
    tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse = " ") ) )
    nullfit <- glmmTMB(tmp_formula, se = FALSE, data = data.frame(newdata, response = rnorm(nrow(newdata)))) #' This may not work in the future if smootihng terms are included, say, due to the standardization that needs to be applied    
    X <- model.matrix(nullfit)[,-1] # Remove the intercept term
    num_units <- nrow(X)
    
    get_eta <- tcrossprod(X, object$betas)
    if(!is.null(newoffset))
        get_eta <- get_eta + newoffset

    
    if(type == "archetype") {
        pt_pred <- get_eta 
        rownames(pt_pred) <- rownames(newdata)
        colnames(pt_pred) <- names(object$mixture_proportion)
        }
    
    if(type %in% c("species_max", "species_mean")) {
        all_spp_mu <- array(NA, dim = c(num_units, num_spp, object$num_archetypes),
                                 dimnames = list(units = rownames(newdata), spp = names(object$spp_intercepts), archetype = names(object$mixture_proportion)))
        for(k0 in 1:object$num_archetypes) {
            all_spp_mu[,,k0] <- matrix(object$spp_intercepts, nrow = num_units, ncol = num_spp, byrow = TRUE) + matrix(get_eta[,k0], nrow = num_units, ncol = num_spp, byrow = FALSE)
            if(!is.null(newoffset))
                all_spp_mu[,,k0] <- all_spp_mu[,,k0] + newoffset
            }    
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
    #' Currently the posterior probabilities of belonging to species belonging to archetypes, but don't know; see (https://github.com/skiptoniam/ecomix/issues/36)
    ##-----------------------
    if(se_fit) {
        construction_predictions_per_bootsrap <- function(k0) {
            cw_bootstrap_parameters <- object$bootsrap_parameters[k0,]
            cw_spp_intercepts <- cw_bootstrap_parameters[grep("spp_intercept", names(cw_bootstrap_parameters))]
            cw_mixture_proportions <- cw_bootstrap_parameters[grep("mixture_proportion", names(cw_bootstrap_parameters))]
            cw_betas <- cw_bootstrap_parameters[grep("beta[1-9]", names(cw_bootstrap_parameters))]             
            cw_betas <- matrix(cw_betas, nrow = object$num_archetypes, byrow = TRUE)

            get_cw_eta <- tcrossprod(X, cw_betas)
            if(!is.null(newoffset))
                get_eta <- get_eta + newoffset
            
            if(type == "archetype") {
                pt_pred <- get_cw_eta 
                rownames(pt_pred) <- rownames(newdata)
                colnames(pt_pred) <- names(object$mixture_proportion)
                }
            
            if(type %in% c("species_max", "species_mean")) {
                all_spp_mu <- array(NA, dim = c(num_units, num_spp, object$num_archetypes),
                                    dimnames = list(units = rownames(newdata), spp = names(object$spp_intercepts), archetype = names(object$mixture_proportion)))
                for(k0 in 1:object$num_archetypes) {
                    all_spp_mu[,,k0] <- matrix(cw_spp_intercepts, nrow = num_units, ncol = num_spp, byrow = TRUE) + matrix(get_eta[,k0], nrow = num_units, ncol = num_spp, byrow = FALSE)
                    if(!is.null(newoffset))
                        all_spp_mu[,,k0] <- all_spp_mu[,,k0] + newoffset
                    }    
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
        
        all_predictions <- foreach(k0 = 1:nrow(object$bootsrap_parameters)) %dopar% construction_predictions_per_bootsrap(k0 = k0)
        all_predictions <- abind::abind(all_predictions, along = 3)
        
        ci_alpha <- (1 - coverage)/2
        quantile_predictions <- apply(all_predictions, c(1,2), collapse::fquantile, probs = c(ci_alpha, 1 - ci_alpha))
        lower_predictions <- quantile_predictions[1,,]
        upper_predictions <- quantile_predictions[2,,]
        rm(quantile_predictions)
        
        if(type == "archetype") {
            rownames(lower_predictions) <- rownames(upper_predictions) <- rownames(newdata)
            colnames(lower_predictions) <- colnames(upper_predictions) <- names(object$mixture_proportion)
            }
        if(type == "species_max") {
            rownames(lower_predictions) <- rownames(upper_predictions) <- rownames(newdata)
            colnames(lower_predictions) <- colnames(upper_predictions) <- names(object$spp_intercepts)
            }
        
        return(list(point_prediction = pt_pred, 
                    lower = lower_predictions,
                    upper = upper_predictions))
        }
    }


 

