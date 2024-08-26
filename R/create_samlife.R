#' @title Simulate data from a SAM
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Simulates multivariate abundance data based on a species archetype model and given the various parameter values as appropriate.

#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param formula An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). *Note there should be nothing on the left hand side of the "~". It should also include an intercept term.*
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param override_X Allows the user to directly supply the model matrix instead of having the function infer it from \code{formula} and \code{data}. This may be useful to speed things, especially if this function has to be applied repeatedly. If the argument is supplied, then anything supplied to \code{formula} and \code{data} is ignored.
#' @param offset A matrix of offset terms.  
#' @param betas A matrix of archetypal regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of archetypes.
#' @param spp_intercepts A vector of species-specific intercept.  
#' @param spp_dispparam A vector of species-specific dispersion parameters, to be used for distributions that require one.  
#' @param spp_powerparam A vector of species-specific power parameters, to be used for distributions that require one. 
#' @param mixture_proportion A vector of mixture proportions corresponding to the probability of belonging to each archetype.
#' @param trial_size Trial sizes to use for binomial distribution. This should equal to a scalar.
#' @param archetype_label If desired, the user can manually supply the archetype labels for each species. In this case, \code{mixture_proportion} must still be supplied but is subsequently ignored.

#' @details 
#' Simulates multivariate abundance data from a species archetype model (SAM). For the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_k,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_k} denotes the corresponding regression coefficients for archetype \eqn{k}. Based on the mean model given above, responses \eqn{y_{ij}} are then simulated from the assumed distribution, using the additional dispersion and power parameters as appropriate.
#' 
#' \subsection{Distributions}{
#' 
#' Currently the following response distributions are permitted: 
#' \describe{
#' \item{\code{Beta()}:}{Beta distribution using a logit link. The corresponding mean-variance relationship is given by \eqn{V = \mu(1-\mu)/(1+\phi)} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}
#' \item{\code{binomial()}:}{Binomial distribution. The corresponding mean-variance relationship is given by \eqn{V = N_{trial}\mu(1-\mu)} where \eqn{\mu} denotes the mean and \eqn{N_{trial}} is the trial size.}
#' \item{\code{Gamma()}:}{Gamma distribution, noting only the log link is permitted. The corresponding mean-variance relationship is given by \eqn{V = \phi\mu^2} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}
#' \item{\code{gaussian()}:}{Gaussian or normal distribution. The corresponding mean-variance relationship is given by \eqn{V = \phi^2}, where \eqn{\phi} is the standard deviation.}
#' \item{\code{poisson()}:}{Poisson distribution. The corresponding mean-variance relationship is given by \eqn{V = \mu} where \eqn{\mu} denotes the mean.}
#' \item{\code{nbinom2()}:}{Negative binomial distribution using a log link. The corresponding mean-variance relationship is given by \eqn{V = \mu + \mu^2/\phi} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}
#' \item{\code{tweedie()}:}{Tweedie distribution using a log link. The corresponding mean-variance relationship is given by \eqn{V = \phi\mu^{\rho}} where \eqn{\mu} denotes the mean, \eqn{\phi} is the dispersion parameter, and \eqn{\rho} is the power parameter.}
#' }
#' }
#' 
#' @return 
#' A list with the following components:
#' \describe{
#' \item{y:}{The simulated multivariate abundance response matrix.}
#' \item{archetype_label:}{A vector of archetype labels for each species.}
#' \item{linear_predictor:}{The matrix of linear predictors corresponding to the simulated multivariate abundance response matrix.}
#' }
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @examples
#' \dontrun{
#' #' ##----------------------
#' # Generate some multivariate abundance data from a SAM
#' ##----------------------
#' library(tidyverse)
#' library(mvtnorm)
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
#' mixture_proportion = true_mixprop)
#' }
#' 
#' @export
#' @import Matrix
#' @importFrom sdmTMB sdmTMB nbinom2 tweedie Beta
#' @importFrom stats as.formula binomial model.matrix rbeta rbinom rgamma rnorm rnbinom rpois plogis
#' @importFrom tweedie rtweedie
#' @md

## Data generation mechanism - intercept included by default
create_samlife <- function(family = binomial(), 
                           formula, 
                           data, 
                           override_X = NULL, 
                           offset = NULL, 
                           betas, 
                           spp_intercepts, 
                           spp_dispparam = NULL, 
                           spp_powerparam = NULL, 
                           mixture_proportion, 
                           trial_size = 1, 
                           archetype_label = NULL) {
    
    formula <- .check_X_formula(formula = formula, data = as.data.frame(data))          
    tmp_formula <- as.formula(paste("response", paste(as.character(formula),collapse = " ") ) )
    if(is.null(override_X)) {
        nullfit <- sdmTMB(tmp_formula, 
                           spatial = FALSE,
                           data = data.frame(data, response = rnorm(nrow(data))))
        X <- model.matrix(nullfit$formula[[1]], data = nullfit$data)[,-1] # Remove the intercept term
        }
    if(!is.null(override_X))
        X <- override_X
    
    
    if(nrow(betas) != length(mixture_proportion)) 
        stop("No. of mixing proportions should be equal to the number of rows in betas.")
    if(ncol(betas) != ncol(X)) 
        stop("No. of columns in betas should be equal to the number of columns in X. Note the latter should not contain an intercept term.")
    if(any(mixture_proportion < 0) || any(mixture_proportion > 1) || abs(sum(mixture_proportion) - 1) > 1e-12)
        stop("The mixture proportions mixture_proportion should be a vector with all elements between 0 and 1, and should sum to 1.")
    if(!(family$family %in% c("gaussian","Gamma","binomial","poisson","nbinom2","tweedie","Beta")))
        stop("family currently not supported. Sorry!")
    
    
    if(is.null(rownames(X)))
        rownames(X) <- paste0("unit", 1:nrow(X))
    if(is.null(colnames(X)))
        colnames(X) <- paste0("x", 1:ncol(X))
    num_units <- nrow(X)
    num_spp <- length(spp_intercepts)
    num_archetypes <- length(mixture_proportion)
    
    resp <- matrix(0, nrow = num_units, ncol = num_spp)
    rownames(resp) <- rownames(X)
    colnames(resp) <- paste0("spp", 1:num_spp)
    if(is.null(archetype_label))
        archetype_label <- sample(1:num_archetypes, size = num_spp, prob = mixture_proportion, replace = TRUE)
    
    spp_eta <- tcrossprod(cbind(1, X), cbind(spp_intercepts, betas[archetype_label,]))
    if(!is.null(offset))
        true_eta <- true_eta + offset
    
    for(j in 1:num_spp) {
        if(family$family == "Beta")
            resp[,j] <- rbeta(num_units, shape1 = spp_dispparam[j]*binomial()$linkinv(spp_eta[,j]), shape2 = spp_dispparam[j]*(1-binomial()$linkinv(spp_eta[,j])))
        if(family$family == "binomial")
            resp[,j] <- rbinom(num_units, size = trial_size, prob = family$linkinv(spp_eta[,j]))
        if(family$family == "gaussian")
            resp[,j] <- rnorm(num_units, mean = family$linkinv(spp_eta[,j]), sd = sqrt(spp_dispparam[j]))
        if(family$family == "Gamma")
            resp[,j] <- rgamma(num_units, scale = exp(spp_eta[,j])*spp_dispparam[j], shape = 1/spp_dispparam[j])
        if(family$family == "nbinom2")
            resp[,j] <- rnbinom(num_units, mu = exp(spp_eta[,j]), size = spp_dispparam[j])
        if(family$family == "poisson")
            resp[,j] <- rpois(num_units, lambda = family$linkinv(spp_eta[,j]))
        if(family$family == "tweedie") {
            resp[,j] <- tweedie::rtweedie(num_units, mu = exp(spp_eta[,j]), phi = spp_dispparam[j], power = spp_powerparam[j])
            }
        }
    
    return(list(y = resp, archetype_label = archetype_label, linear_predictor = spp_eta)) 
    }
