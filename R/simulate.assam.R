#' @title Simulate data from approximate and scalable SAM object
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#'
#' Simulate multivariate abundance data based on a fitted \code{assam} object.
#' 
#' @param object An object of class \code{assam}.
#' @param data The data frame used as part of the \code{assam} object. This needs to be supplied since currently \code{assam} objects do not save the data to save memory. 
#' @param nsim A positive integer specifying the number of simulated datasets. Defaults to 1.
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
#' @return A matrix with two rows and \code{nsim} columns, where the first row contains the simulated multivariate abundance response matrix, and the second row contains the vector of archetype labels for each species.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' @examples
#' \dontrun{
#' #' Please see the help file for assam for example.
#' }
#' 
#' @export
#' @importFrom stats as.formula plogis
#' @md


simulate.assam <- function(object, data, nsim = 1, seed = NULL, ...) {
    if(class(object) != "assam") 
        stop("object must be of class assam.")

    use_model <- list(family = object$family, 
                     formula = object$formula,
                     data = data,
                     offset = object$offset,
                     betas = object$betas, 
                     spp_intercepts = object$spp_intercepts, 
                     mixture_proportion = object$mixture_proportion, 
                     trial_size = object$trial_size, 
                     spp_dispparam = NULL, 
                     spp_powerparam = NULL) 
    if(object$family$family[1] %in% c("Beta","nbinom2","Gamma","gaussian","tweedie"))
        use_model$spp_dispparam <- object$spp_nuisance$dispersion
    if(object$family$family[1] %in% c("tweedie"))
        use_model$spp_powerparam <- object$spp_nuisance$power
    
    if(!is.null(seed) || length(seed) > 0) 
        set.seed(seed[1])
    
    
    out <- replicate(nsim, 
                     create_samlife(family = use_model$family, 
                                    formula = use_model$formula, 
                                    data = data,
                                    betas = use_model$betas, 
                                    spp_intercepts = use_model$spp_intercepts, 
                                    spp_dispparam = use_model$spp_dispparam, 
                                    spp_powerparam = use_model$spp_powerparam, 
                                    mixture_proportion = use_model$mixture_proportion, 
                                    trial_size = use_model$trial_size))
    
    set.seed(NULL)
    return(out)
    }   


#' @rdname simulate.assam
#' @export 
simulate <- function(object, ...) {
        UseMethod("simulate")
    }

