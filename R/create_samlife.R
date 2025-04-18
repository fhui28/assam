#' @title Simulate data from a SAM
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Simulates multivariate abundance data based on a species archetype model and given the various parameter values as appropriate.

#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param formula An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). *Note there should be nothing on the left hand side of the "~".*
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param which_spp_effects A vector identifying which columns of the model matrix induced by \code{formula} and \code{data} should be treated as species-specific effects. Default to 1, meaning only the first column i.e., the intercept, is species-specific.
#' @param offset A matrix of offset terms.  
#' @param betas A matrix of archetypal regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of archetypes.
#' @param spp_effects A matrix of species-specifics, where the number of rows defines the number of species in the simulated dataset, and the number of columns should equal the length of \code{which_spp_effects}.
#' @param spp_dispparam A vector of species-specific dispersion parameters, to be used for distributions that require one.  
#' @param spp_powerparam A vector of species-specific power parameters, to be used for distributions that require one. 
#' @param mesh Output from [sdmTMB::make_mesh()], from which species-specific spatial fields can be constructed and added to the linear predictor.
#' @param spp_spatial_sd A vector of standard deviations corresponding to the species-specific spatial fields. This corresponds to the \eqn{\sigma} parameters described in [sdmTMB's parametrization of Gaussian random fields](https://pbs-assess.github.io/sdmTMB/articles/model-description.html#gaussian-random-fields).
#' @param spp_spatial_range A vector of range parameters corresponding to the species-specific spatial fields. This is equal to \eqn{1/\kappa} where \eqn{\kappa} is described in [sdmTMB's parametrization of Gaussian random fields](https://pbs-assess.github.io/sdmTMB/articles/model-description.html#gaussian-random-fields).
#' @param mixing_proportion A vector of mixture proportions corresponding to the probability of belonging to each archetype.
#' @param trial_size Trial sizes to use for binomial distribution. This should equal to a scalar.
#' @param archetype_label If desired, the user can manually supply the archetype labels for each species. In this case, \code{mixing_proportion} must still be supplied but is subsequently ignored.
#' @param seed A seed that can be set for simulating datasets.
#'
#' @details 
#' Simulates multivariate abundance data from a species archetype model (SAM). For the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = u_i^\top\alpha_j + x_i^\top\beta_k,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{u_i^\top\alpha_j} corresponds to a component that is to kept species-specific e.g., species-specific intercept, \eqn{x_i^\top\beta_k}  denotes the component corresponding to effect of archetypal response \eqn{k}. Additionally, species-specific spatial fields can be included in the linear predictor e.g., to account for residual spatial correlation above and beyond that explained by the archetypal responses.
#' 
#' Based on the mean model given above, responses \eqn{y_{ij}} are then simulated from the assumed distribution, using the additional dispersion and power parameters as appropriate.
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
#' \item{spatial_fields:}{If applicable, a matrix of species-specific spatial fields.}
#' }
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
#'
#' set.seed(022025)
#' 
#' num_X <- 10
#' num_units <- 1000
#' num_spp <- 100
#' num_archetype <- 5
#' H <- outer(1:num_X, 1:num_X, "-")
#' H <- 0.5^abs(H)
#' covariate_dat <- data.frame(x = runif(num_units), y = runif(num_units))
#' covariate_dat <- bind_cols(covariate_dat,
#'     rmvnorm(num_units, sigma = H) %>% 
#'     as.data.frame %>% 
#'     rename_with(., .fn = function(x) paste0("covariate", x)))
#' rm(H)
#'
#' true_betas <- runif(num_archetype * num_X, -1, 1) %>% matrix(nrow = num_archetype)
#' true_species_effects <- matrix(runif(num_spp, -2, 0), ncol = 1)
#' true_dispparam <- 1/runif(num_spp, 0, 5) 
#' true_powerparam <- runif(num_spp, 1.4, 1.8)
#' true_mixprop <- c(0.2, 0.2, 0.3, 0.15, 0.15)
#' true_spatial_sd <- runif(num_spp, 0.2, 1)
#' true_spatial_range <- runif(num_spp, 0.1, 0.3)
#' 
#' 
#' simdat <- create_samlife(family = nbinom2(), 
#' formula = paste("~", paste0(colnames(covariate_dat)[-(1:2)], collapse = "+")) %>% as.formula, 
#' data = covariate_dat, 
#' betas = true_betas, 
#' spp_effects = true_species_effects, 
#' spp_dispparam = true_dispparam, 
#' spp_powerparam = true_powerparam, 
#' #mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 80),
#' #spp_spatial_sd = true_spatial_sd,
#' #spp_spatial_range = true_spatial_range,
#' mixing_proportion = true_mixprop,
#' seed = 022025)
#' }
#' 
#'
#' @export
#' @import Matrix
#' @importFrom sdmTMB make_mesh sdmTMB sdmTMB_simulate nbinom2 tweedie Beta
#' @importFrom stats as.formula binomial rbeta rbinom rgamma rnorm rnbinom rpois plogis
#' @importFrom tweedie rtweedie
#' @md

create_samlife <- function(family = binomial(), 
                           formula, 
                           data, 
                           which_spp_effects = 1,
                           offset = NULL, 
                           betas, 
                           spp_effects, 
                           spp_dispparam = NULL, 
                           spp_powerparam = NULL,
                           mesh = NULL,
                           spp_spatial_sd = NULL,
                           spp_spatial_range = NULL,
                           mixing_proportion,
                           trial_size = 1, 
                           archetype_label = NULL,
                           seed = NULL) {
    
    ##----------------
    #' # Checks and balances
    ##----------------
    if(!is.matrix(spp_effects))
        stop("spp_effects must be matrix, where the number of rows defines the number of species in the simulated dataset.")
    num_units <- nrow(data)
    num_spp <- nrow(spp_effects)
    num_archetypes <- length(mixing_proportion)
    formula <- .check_X_formula(formula = formula, data = as.data.frame(data))          
    
    if(nrow(betas) != length(mixing_proportion))
        stop("No. of mixing proportions should be equal to the number of rows in betas.")
    if(ncol(spp_effects) != length(which_spp_effects))
        stop("species_effects should be a matrix where the number of columns should equal the length of which_spp_effects.")
    if(any(mixing_proportion < 0) || any(mixing_proportion > 1) || abs(sum(mixing_proportion) - 1) > 1e-12)
        stop("The mixture proportions mixing_proportion should be a vector with all elements between 0 and 1, and should sum to 1.")
    if(!(family$family %in% c("gaussian","Gamma","binomial","poisson","nbinom2","tweedie","Beta")))
        stop("family currently not supported. Sorry!")
    
    add_spatial <- FALSE
    if(!is.null(mesh)) {
        if(!inherits(mesh, "sdmTMBmesh"))
            stop("If mesh is supplied for species-specific spatial fields, then the mesh argument must be an object class of \"sdmTMBmesh\".")
        add_spatial <- TRUE
        }
    
    ## Check and set spatial parameters
    .check_spp_spatial_parameters(spp_spatial_range = spp_spatial_range,
                                  spp_spatial_sd = spp_spatial_sd,
                                  add_spatial = add_spatial,
                                  spp_effects = spp_effects)
    
    ##----------------
    #' # Simulate data
    ##----------------
    get_spatial_fields <- NULL
    resp <- spp_eta <- matrix(0, nrow = num_units, ncol = num_spp)
    rownames(resp) <- rownames(spp_eta) <- rownames(data)
    colnames(resp) <- colnames(spp_eta) <- paste0("spp", 1:num_spp)
    
    set.seed(seed)
    if(is.null(archetype_label))
        archetype_label <- sample(1:num_archetypes, size = num_spp, prob = mixing_proportion, replace = TRUE)
    set.seed(NULL)
    
    
    #' ## For non-spatial models, skip using sdmTMB_simulate as it is *much* slower!
    if(!add_spatial) {
        tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ")))
        nullfit <- sdmTMB(tmp_formula,
                          spatial = FALSE,
                          data = data.frame(data, response = rnorm(nrow(data))))
        X <- model.matrix(nullfit$formula[[1]], data = nullfit$data)
        rm(nullfit)
        
        spp_eta <- tcrossprod(X[, which_spp_effects, drop = FALSE], spp_effects) + tcrossprod(X[, -which_spp_effects, drop = FALSE],  betas[archetype_label,])
        if(!is.null(offset))
            spp_eta <- spp_eta + offset
        }
    # get_spatial_fields <- NULL
    # if(!is.null(spp_spatial_sd)) {
    #     for(j in 1:num_spp) {
    #         get_spatial_fields <- cbind(get_spatial_fields,
    #                                     grf(grid = spatial_coordinates,
    #                                         nsim = 1,
    #                                         cov.model = "matern",
    #                                         cov.pars = c(spp_spatial_sd[j]^2, spp_spatial_range[j]), # grf parametrizes matern in terms of variance (sd^2) and spatial range phi, where the former is the same as tau and the latter is 1/kappa in (https://pbs-assess.github.io/sdmTMB/articles/model-description.html#gaussian-random-fields). 
    #                                         kappa = 1, # Smoothness set to 1 to exploit SPDE approach for estimation via sdmTMB later on
    #                                         messages = FALSE)$data
    #                                     ) 
    #         }
    #     
    #     spp_eta <- spp_eta + get_spatial_fields
    #     }
    
    for(j in 1:num_spp) {
        cw_seed <- NULL
        if(!is.null(seed))
            cw_seed <- seed + j
        
        if(!add_spatial) {
            set.seed(cw_seed)
            
            if(family$family[1] == "Beta")
                resp[,j] <- rbeta(num_units, shape1 = spp_dispparam[j]*binomial()$linkinv(spp_eta[,j]), shape2 = spp_dispparam[j]*(1-binomial()$linkinv(spp_eta[,j])))
            if(family$family[1] == "binomial")
                resp[,j] <- rbinom(num_units, size = trial_size, prob = family$linkinv(spp_eta[,j]))
            if(family$family[1] == "gaussian")
                resp[,j] <- rnorm(num_units, mean = family$linkinv(spp_eta[,j]), sd = spp_dispparam[j])
            if(family$family[1] == "Gamma")
                resp[,j] <- rgamma(num_units, scale = exp(spp_eta[,j])*spp_dispparam[j], shape = 1/spp_dispparam[j])
            if(family$family[1] == "nbinom2")
                resp[,j] <- rnbinom(num_units, mu = exp(spp_eta[,j]), size = spp_dispparam[j])
            if(family$family[1] == "poisson")
                resp[,j] <- rpois(num_units, lambda = family$linkinv(spp_eta[,j]))
            if(family$family[1] == "tweedie") {
                resp[,j] <- tweedie::rtweedie(num_units, mu = exp(spp_eta[,j]), phi = spp_dispparam[j], power = spp_powerparam[j])
                }
            
            set.seed(NULL)
            }
        
        if(add_spatial) {
            make_offset <- NULL
            if(!is.null(offset))
                make_offset <- offset[,j]
            
            make_B <- c(spp_effects[j,], betas[archetype_label[j],])
            make_B[which_spp_effects] <- spp_effects[j,]
            make_B[-which_spp_effects] <- betas[archetype_label[j],]
            sim_resp <- sdmTMB::sdmTMB_simulate(formula = formula,
                                                data = data,
                                                mesh = mesh,
                                                family = family,
                                                B = make_B,
                                                range = sqrt(8) * spp_spatial_range[j],
                                                sigma_O =  spp_spatial_sd[j],
                                                tweedie_p = spp_powerparam[j],
                                                phi = spp_dispparam[j],
                                                seed = cw_seed,
                                                offset = make_offset)
            
            resp[,j] <- sim_resp$observed
            spp_eta[,j] <- sim_resp$eta
            rm(make_B)
            if(add_spatial)
                get_spatial_fields <- cbind(get_spatial_fields, sim_resp$omega_s)
            }
        }
    
    if(add_spatial) {
        rownames(get_spatial_fields) <- rownames(data)
        colnames(get_spatial_fields) <- paste0("spp", 1:num_spp)
        }
    
    return(list(y = resp, 
                archetype_label = archetype_label, 
                linear_predictor = spp_eta, 
                spatial_fields = get_spatial_fields)) 
    }
