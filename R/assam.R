#' @title Approximate and scalable species archetype models (asSAMs)
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Fits approximate and scalable species archetype modeling (asSAMs) for model-based clustering of species based on their environmental response, into a small number of so-called archetypal responses. The basic idea is take the log-likelihood function of SAM, and then construct an approximation of this which (hopefully) is more scalable in both the number of sites and species.
#' 
#' @param y A multivariate abundance response matrix.
#' @param formula An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). *Note there should be nothing on the left hand side of the "~". It should also include an intercept term.*
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param offset A matrix of offset terms, of the same dimension as \code{y}.
#' @param trial_size Trial sizes to use for binomial distribution. This should equal to a scalar.
#' @param num_archetypes Number of archetypes (clusters) to assume in the asSAM.
#' @param mesh Output from [sdmTMB::make_mesh()], used for adding species-specific spatial fields to the linear predictor.
#' @param do_parallel Should parallel computing be used to fit the asSAM. Defaults to \code{FALSE}.
#' @param num_cores If \code{do_parallel = TRUE}, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to \code{parallel::detectCores() - 2}.
#' @param uncertainty_quantification Should uncertainly intervals be computed via parametric bootstrap?
#' @param control A list containing the following elements:
#' \itemize{
#' \item{max_iter:}{the maximum number of iterations in the EM algorithm. Usually convergence is quite quick e.g., less than 20 iterations.}
#' \item{tol:}{the convergence criterion; the difference in the log-likelihood value of the asSAM from successive iterations must be smaller than this value.}
#' \item{temper_prob:}{In the iteration of the EM algorithm, posterior probabilities from the E-step are "tempered" or push away from the 0/1 boundary. This is often useful to get the EM algorithm moving initially.}
#' #' \item{trace:}{controls if messages are printed as part of the estimation process to reflect progress.}
#' }
#' @param bootstrap_control A list containing the following elements to control the parametric bootstrap for uncertainty quantification:
#' \itemize{
#' \item{num_boot:}{the number of bootstrapped iterations to do. Defaults to 100, which can already take a long time but should be enough in a lot of settings for uncertainty quantification.}
#' \item{ci_alpha:}{the type-1 level for confidence interval construction. \code{100 * (1 - ci_alpha)} percent Confidence intervals are constructed.}
#' \item{seed:}{a seed that can be set for bootstrapping the datasets.}
#' #' \item{ci_type:}{type of confidence intervals to construct. Two options are currently available: 1) "percentile" confidence intervals based directly on the empirical quantiles of the bootstrap samples; 2) "expanded" percentile confidence intervals, which are typically slightly wider intervals that attempt to correct for a so-called "narrowness bias". Defaults to "percentile".}
#' }

#' @details 
#' or the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_k,}
#' 
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_k} denotes the corresponding regression coefficients for archetype \eqn{k}. Based on the mean model given above, responses \eqn{y_{ij}} are then simulated from the assumed distribution, using the additional dispersion and power parameters as appropriate. We refer the reader to [Dunstan et al., (2011)](https://doi.org/10.1016/j.ecolmodel.2010.11.030), [Hui et al., (2013)](https://doi.org/10.1890/12-1322.1), [Dunstan et al., (2013)](https://doi.org/10.1007/s13253-013-0146-x), and [Skipton Woolley's ecomix package](https://github.com/skiptoniam/ecomix) for more details about the formulations of SAMs.
#' 
#' The broad goal of this package is to construct a way of fitting SAMs that are, although approximate, more scalable in the number of sites and species (though not necessarily faster), hence the name asSAMs. We prefer to the corresponding manuscript (in preparation) for details, but to summarize, asSAMs are formed by constructing an approximate likelihood function for a SAM based on using ingredients (i.e., point estimates and the associated observed information matrix) from stacked species distribution models (which are fitted initially in parallel), and then building what is essentially a finite mixture of multivariate Gaussian distributions from this. This is then maximized using an EM algorithm, which can be done scalably and very quickly.
#' 
#' For uncertainty quantification, a parametric bootstrap approach is taken along the lines of Section 2.16.2 in [McLachlan and Peel (2004)](https://www.wiley.com/en-us/Finite+Mixture+Models-p-9780471654063). While computationally slow, it is simply in design and often does a decent job for mixture modeling! 
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
#' @section A note on parallelization: 
#' The scalability of asSAMs relies on being able to deploy [sdmTMB::sdmTMB()] is an efficient and faily optimized manner. Along these lines, please see [using sdmTMB in parallel](https://github.com/pbs-assess/sdmTMB/issues/368) and [using OpenBLAS](https://gist.github.com/seananderson/08a51e296a854f227a908ddd365fb9c1) and references therein for some tips to ensure things on your machine are more optimized. Thanks to Sean Anderson for this advice!
#' 
#' 
#' @return An object of class \code{assam} with the following elements (as appropriate, and not necessarily in the order below):
#' \item{call:}{The function call.}
#' \item{formula:}{Same as input argument.}
#' \item{family:}{Same as input argument.}
#' \item{num_nuisance_perspp:}{The number of "nuisance" parameters per species. For example, if \code{family = binomial()} and species-specific spatial fields are not included, then there are no nuisance parameters. If \code{family = nbinom2()} and species-specific spatial fields are included, say, then there are three nuisance parameters per species.}
#' \item{trial_size:}{Same as input argument.}
#' \item{offset:}{Same as input argument.}
#' \item{mesh:}{Same as input argument.}
#' \item{num_archetypes:}{Same as input argument.}
#' \item{uncertainty_quantification:}{Same as input argument.}
#' \item{spp_intercepts:}{Estimated species-specific intercepts.}
#' \item{betas:}{Estimated matrix of archetypal regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of archetypes.}
#' \item{mixture_proportion:}{Estimated vector of mixture proportions corresponding to the probability of belonging to each archetype.}
#' \item{spp_nuisance:}{Estimated matrix of species-specific nuisance parameters e.g., the dispersion parameter in the negative binomial distribution, and dispersion and power parameters in the Tweedie distribution, and so on.}
#' \item{posterior_probability:}{Estimated matrix of posterior probabilities for each species belong to each archetype. The number of rows in \code{posterior_probability} is equal to the number of species.}
#' \item{linear_predictor:}{Estimated array of archetype-specific linear predictors. The last dimension of the array corresponds to the number of archetypes.}
#' \item{spatial_fields:}{Predicted matrix species-specific spatial fields, if included.}
#' \item{logL:}{Estimated log-likelihood value of the asSAM i.e. the value of the approximated log-likelihood function at convergence.}
#' \item{df:}{Number of the estimated (freely varying) parameters in the asSAM.}
#' \item{linear_predictor:}{Estimated array of archetype-specific linear predictors. The last dimension of the array corresponds to the number of archetypes.}
#' \item{control:}{Same as input argument.}
#' \item{bootstrap_control:}{Same as input argument.}
#' \item{confidence_intervals:}{A list of estimated parametric bootstrap confidence intervals for all the parameters in the asSAM, with each element in the list corresponding to one of the parameters e.g., the species-specific intercepts, the mixture proportions, and so on.}
#' \item{bootstrap_paramters:}{A matrix of estimated (untransformed) parameters from the parametric bootstrap. *This output can be safely ignored by most users*, but those curious, it is used to uncertainty quantification for downstream predictions, say.}
#' \item{bootstrap_posterior_probability:}{An array of estimated posterior probabilities from the parametric bootstrap. *This output can be safely ignored by most users*, but those curious, it is used to uncertainty quantification for downstream predictions, say.}
#' \item{sdmTMB_fits}{A list of the set of stacked species distribution models fitted using [sdmTMB::sdmTMB()]. *This output can be safely ignored by most users*.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
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
#' true_intercepts <- runif(num_spp, -2, 0)  
#' true_dispparam <- 1/runif(num_spp, 1, 5) 
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
#' samfit <- assam(y = simdat$y,
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' family = nbinom2(),
#' uncertainty_quantification = TRUE,
#' num_archetypes = num_archetype,
#' num_cores = 8)
#' 
#'  
#' plot(true_intercepts, samfit$spp_intercepts); abline(0,1)
#' plot(true_dispparam, samfit$spp_nuisance$dispersion, log = "xy"); abline(0,1)
#' #' Note estimates for the archetypal responses and mixture proportions from (as)SAMs should be 
#' #' close to the corresponding true values, *up to a reordering* of the mixture component
#' #' s/archetypes (since the order is essentially arbitrary)
#' rbind(true_betas, samfit$betas) %>% 
#' t %>% 
#' as.data.frame %>%
#' ggpairs
#' table(simdat$archetype_label, apply(samfit$posterior_probability, 1, which.max))
#' 
#'  
#' ##----------------------
#' # Demonstrating basic use of functions for asSAM 
#' ##----------------------
#' samfit
#' summary(samfit)
#' 
#' fitted(samfit)
#'  
#' simulate(samfit)
#' 
#' residuals(samfit, type = "dunnsmyth")
#'  
#' #' Basic residual analysis
#' plot(samfit, y = simdat$y, transform_fitted_values = TRUE)
#'  
#' #' Archetype-level predictions
#' predict(samfit, newdata = covariate_dat, type = "archetype", se_fit = TRUE) 
#' 
#' #' Species-level predictions
#' predict(samfit,  newdata = covariate_dat, type = "species_max", num_cores = 8, se_fit = TRUE) 
#'  
#' }
#' 
#' 
#' 
#' @export assam
#'
#' @importFrom abind abind
#' @importFrom cluster pam
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom sdmTMB sdmTMB make_mesh nbinom2 tweedie Beta
#' @importFrom label.switching pra
#' @importFrom methods as
#' @importFrom parallel detectCores
#' @importFrom stats as.formula nlminb model.matrix qt plogis qlogis 
#' @importFrom TMB MakeADFun
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import Matrix
#' @md

assam <- function(y, 
                  formula, 
                  data, 
                  family, 
                  offset = NULL, 
                  trial_size = 1, 
                  num_archetypes,
                  mesh = NULL, 
                  do_parallel = TRUE, 
                  num_cores = NULL, 
                  uncertainty_quantification = TRUE, 
                  control = list(max_iter = 500, tol = 1e-5, temper_prob = 0.8, trace = FALSE),
                  bootstrap_control = list(num_boot = 100, ci_alpha = 0.05, seed = NULL, ci_type = "percentile")) {

    
    ##----------------
    # Checks and balances
    ##----------------
    if(!(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "poisson", "binomial", "tweedie")))
        stop("family currently not supported. Sorry!")

    if(length(trial_size) > 1)
        stop("Currently the size argument only supports a scalar. Sorry!")

    if(is.null(colnames(y)))
        colnames(y) <- paste0("spp", 1:ncol(y))
    if(is.null(rownames(y)))
        rownames(y) <- paste0("unit", 1:nrow(y))
          
    if(is.null(num_cores))
        registerDoParallel(cores = detectCores() - 2)
    if(!is.null(num_cores))
        registerDoParallel(cores = num_cores)
     
    .check_offset(offset = offset, y = y) 
    if(is.null(offset))
        offset <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    formula <- .check_X_formula(formula = formula, data = as.data.frame(data))     
    tmp_formula <- as.formula(paste("response", paste(as.character(formula),collapse = " ") ) )
    nullfit <- sdmTMB(tmp_formula, #' Currently not even sure you need this...at least it is not needed for spatial models!
                      spatial = FALSE, 
                      data = data.frame(data, response = rnorm(nrow(data))))
    X <- model.matrix(nullfit$formula[[1]], data = nullfit$data)[,-1] # Remove the intercept term
    num_X <- ncol(X)
    rm(tmp_formula, nullfit) 
    
    if(!is.null(mesh)) {
        if(class(mesh) != "sdmTMBmesh")
            stop("If mesh is supplied for species-specific spatial fields, then the mesh argument must be an object class of \"sdmTMBmesh\".")
        }
    
    control <- .fill_control(control = control)
    bootstrap_control <- .fill_bootstrap_control(control = bootstrap_control)

    num_unit <- nrow(y)
    num_spp <- ncol(y)
          

    ##----------------
    # Construct quadratic approximations for each species, and set up relevant quantities
    ##----------------
    message("Commencing fitting...")
    if(control$trace)
        message("Forming local quadratic approximation...")
     
    get_qa <- .quadapprox2_fn(family = family,
                             formula = formula, 
                             resp = y, 
                             data = data, 
                             mesh = mesh,
                             offset = offset,
                             trial_size = trial_size,
                             do_parallel = do_parallel) 
    get_qa$long_parameters <- apply(get_qa$parameters, 1, function(x) kronecker(rep(1,num_archetypes), x)) # Repeats species-specific estimates num_archetypes times; object has num_spp columns
    num_nuisance_perspp <- length(get_qa$parameters[1,]) - num_X - 1 # e.g., number of dispersion and power parameter per-species, along with parameters for species-specific spatial fields
    
    
    ### Make mapping matrix, maps psi to long_parameters
    #' Psi sequence: species-specific intercepts; archetypal regression coefficients by archetype; species-specific nuisance parameters, along with parameters for species-specific spatial fields
    mapping_mat <- Matrix(0, nrow = length(get_qa$long_parameters), ncol = num_spp + num_X * num_archetypes + num_spp*num_nuisance_perspp)
    makecolnames <- c(paste0("spp_intercept", 1:num_spp), 
                      paste0(rep(paste0("archetype", 1:num_archetypes), each = num_X), rep(paste0("_beta", 1:num_X), num_archetypes)))
    if(num_nuisance_perspp > 0)
        makecolnames <- c(makecolnames,
                          paste0(rep(paste0("spp_nuisance", 1:num_spp), each = num_nuisance_perspp), rep(colnames(get_qa$parameters)[num_X + 1 + 1:num_nuisance_perspp], num_spp))
                          )
    colnames(mapping_mat) <- makecolnames
    rm(makecolnames)
    
    makerownames_mappingmat_fn <- function(k, j) {
        out <- c(paste0("spp_intercept", j), 
                 paste0("archetype", k, "_beta", 1:num_X))
        if(num_nuisance_perspp > 0)
            out <- c(out, paste0("spp_nuisance", j, colnames(get_qa$parameters)[num_X + 1 + 1:num_nuisance_perspp]))
        
        return(out)
        }
    rownames(mapping_mat) <- as.vector(sapply(1:num_spp, function(j) as.vector(sapply(1:num_archetypes, makerownames_mappingmat_fn, j = j))))
    rm(makerownames_mappingmat_fn)
    for(k0 in 1:ncol(mapping_mat)) {
        mapping_mat[which(rownames(mapping_mat) == colnames(mapping_mat)[k0]), k0] <- 1
        }
    gc()
          
          
    ##----------------
    #' # Run EM algorithm
    ##----------------
    em_fn <- function(qa_object) {
         counter <- 0
         diff <- 10
         cw_logL <- -Inf
         cw_params <- numeric(ncol(mapping_mat))
         post_prob <- matrix(NA, nrow = num_spp, ncol = num_archetypes)
         rownames(post_prob) <- colnames(y)
         colnames(post_prob) <- paste0("archetype", 1:num_archetypes)
         logL_spp <- NULL         
         basic_bigW <- bdiag(lapply(1:num_spp, function(j) kronecker(Diagonal(n = num_archetypes), qa_object$hessian[[j]]))) # Needed as part of constructing the big weight matrix in the M-step
         
         
         while(diff > control$tol & counter < control$max_iter) {
             if(counter == 0) {
                 do_kmeans <- cluster::pam(qa_object$parameters[, 1 + 1:num_X, drop = FALSE], k = num_archetypes)  
                 cw_betas <- do_kmeans$medoids
                 cw_spp_intercept <- qa_object$parameters[, 1]
                 cw_nuisance <- NULL
                 if(num_nuisance_perspp > 0)
                    cw_nuisance <- qa_object$parameters[, num_X + 1 + 1:num_nuisance_perspp, drop = FALSE]
                 cw_mixprop <- as.vector(table(do_kmeans$clustering)) / num_spp
                 rm(do_kmeans)
                 }
             
             ##-------------------
             #' ## E-step
             #' Ignore the normalizing constant of the normal distribution as that does not vary as a function of archetype anyway
             ##-------------------
             for(j in 1:num_spp) { 
                 cw_Quad <- sapply(1:num_archetypes, function(k) {
                     cw_v <- matrix(c(cw_spp_intercept[j], cw_betas[k,], cw_nuisance[j,]) - qa_object$parameters[j,], ncol = 1)
                     return(-0.5 * (crossprod(cw_v, qa_object$hessian[[j]]) %*% cw_v))
                     })
                 eps <- max(cw_Quad)
                 # cw_Quad <- sapply(1:num_archetypes, function(k) {
                 #     out <- dmvnorm(c(cw_spp_intercept[j], cw_betas[k,], cw_nuisance[j,]), mean = qa_object$parameters[j,], sigma = solve(qa_object$hessian[[j]]), log = TRUE)
                 #     return(out)
                 #     })
                 logL_spp[j] <- log(sum(cw_mixprop * exp(cw_Quad - eps))) + eps
                 post_prob[j,] <- exp((log(cw_mixprop) + cw_Quad) - logL_spp[j])
                 rm(eps, cw_Quad)
                 }
               
             new_logL <- sum(logL_spp)
             if(counter == 0) {
                 ## Temper the classification in the initial E-step
                 alpha_temper <- (1 - control$temper_prob * num_archetypes) / (control$temper_prob * (2-num_archetypes) - 1)          
                 for(j in 1:num_spp)
                     post_prob[j,] <- (2*alpha_temper*post_prob[j,]-alpha_temper+1)/(2*alpha_temper - alpha_temper*num_archetypes + num_archetypes)
                 new_logL <- -1e8               
                 }
                    
             ##-------------------
             #' ## M-step 
             ##-------------------
             new_mixprop <- colMeans(post_prob)
             bigW <- Diagonal(x = rep(as.vector(t(sqrt(post_prob))), each = nrow(qa_object$hessian[[j]]))) %*% basic_bigW %*% Diagonal(x = rep(as.vector(t(sqrt(post_prob))), each = nrow(qa_object$hessian[[j]])))
             MtWM <- forceSymmetric(crossprod(mapping_mat, bigW) %*% mapping_mat)
             new_params <- solve(MtWM, crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters))
             new_params <- as.vector(new_params)
             names(new_params) <- colnames(mapping_mat)
             new_spp_intercept <- as.vector(new_params[grep("spp_intercept", names(new_params))])
             new_betas <- matrix(new_params[grep("archetype", names(new_params))], nrow = num_archetypes, byrow = TRUE)
             new_nuisance <- NULL
             if(num_nuisance_perspp > 0) {
                 new_nuisance <- matrix(new_params[grep("spp_nuisance", names(new_params))], nrow = num_spp, byrow = TRUE)
                 colnames(new_nuisance) <- colnames(get_qa$parameters)[num_X + 1 + 1:num_nuisance_perspp]  
                 }
             rm(bigW, MtWM)

             ##-------------------
             #' ## Finish iteration
             ##-------------------
             diff <- new_logL - cw_logL 
             if(control$trace == TRUE) {
                 message("Iteration: ", counter, "\t New Log-likelihood:", round(new_logL, 4), "\t Difference in Log-likelihood: ", round(new_logL - cw_logL, 4))
                 #, "\t Norm difference in parameters: ", round(sum((new_params - cw_params)^2), 5))
                 }
                    
             cw_logL <- new_logL
             cw_mixprop <- new_mixprop
             cw_spp_intercept <- new_spp_intercept
             cw_betas <- new_betas
             cw_nuisance <- new_nuisance
             cw_params <- new_params
             counter <- counter + 1          
             }
         
         names(new_spp_intercept) <- colnames(y)
         rownames(new_betas) <- names(new_mixprop)
         colnames(new_betas) <- colnames(X)
         
         tmp_nuisance <- NULL
         if(num_nuisance_perspp > 0) {
             if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
                 tmp_nuisance$dispersion <- exp(new_nuisance[,grep("ln_phi", colnames(new_nuisance))])
             if(family$family[1] == "tweedie")
                 tmp_nuisance$power <- plogis(new_nuisance[,grep("thetaf", colnames(new_nuisance))]) + 1
             if(!is.null(mesh)) { 
                 tmp_nuisance$spatial_range <- 1/exp(new_nuisance[,grep("ln_kappa", colnames(new_nuisance))])
                 est_tau <- exp(new_nuisance[,grep("ln_tau_O", colnames(new_nuisance))])
                 tmp_nuisance$spatial_SD <- 1/(sqrt(4*pi) * est_tau * exp(new_nuisance[,grep("ln_kappa", colnames(new_nuisance))]))
                 }
                 
             tmp_nuisance <- as.data.frame(tmp_nuisance)
             rownames(new_nuisance) <- rownames(tmp_nuisance) <- colnames(y)            
             }
         
         return(list(new_logL = new_logL, 
                     new_mixprop = new_mixprop, 
                     new_spp_intercept = new_spp_intercept, 
                     new_betas = new_betas, 
                     new_nuisance = new_nuisance, 
                     new_transformed_nuisance = tmp_nuisance,
                     new_params = new_params, 
                     counter = counter, 
                     post_prob = post_prob))          
         }

    if(control$trace)
        message("Commencing EM algorithm...")
    do_em <- em_fn(qa_object = get_qa)

    try_counter <- 0
    if(any(do_em$new_mixprop < 1e-3) & try_counter < 20) {
        message("Mixture component is being emptied...altering initial temp probability and restarting EM-algorithm to try and fix this.")
        control$temper_prob <- control$temper_prob + 0.025
          
        do_em <- em_fn(qa_object = get_qa)
        try_counter <- try_counter + 1
        }

 
    ##----------------
    #' # Format output
    ##----------------
    message("Fitting completed, applying finishing touches...")
    
    out_assam <- list(call = match.call(), 
                      formula = formula,
                      family = family,
                      num_nuisance_perspp = num_nuisance_perspp,
                      trial_size = trial_size,
                      offset = as(offset, "sparseMatrix"),
                      mesh = mesh,
                      num_archetypes = num_archetypes,
                      uncertainty_quantification = uncertainty_quantification,
                      spp_intercepts = do_em$new_spp_intercept,
                      betas = do_em$new_betas,
                      mixture_proportion = do_em$new_mixprop,
                      spp_nuisance = as.data.frame(do_em$new_transformed_nuisance),
                      posterior_probability = do_em$post_prob)

    #' ## Predict species-specific field in an ad-hoc but scalable manner based on the most likely archetype that the species belong.
    get_spatial_fields <- foreach(l = 1:num_spp, .combine = "cbind") %dopar% .predict_spatial_fields(l = l, 
                                                                                                     mesh = mesh,
                                                                                                     qa_object = get_qa, 
                                                                                                     em_object = do_em, 
                                                                                                     family = family)
    if(!is.null(mesh)) {
        rownames(get_spatial_fields) <- rownames(y)
        colnames(get_spatial_fields) <- colnames(y)
        }
        
    get_eta <- tcrossprod(X, out_assam$betas)
    out_assam$linear_predictor <- array(NA, dim = c(num_unit, num_spp, num_archetypes),
                                        dimnames = list(units = rownames(y), spp = colnames(y), archetype = names(out_assam$mixture_proportion)))
    for(k0 in 1:num_archetypes) {
        out_assam$linear_predictor[,,k0] <- matrix(out_assam$spp_intercepts, nrow = num_unit, ncol = num_spp, byrow = TRUE) + matrix(get_eta[,k0], nrow = num_unit, ncol = num_spp, byrow = FALSE) + offset
        if(!is.null(mesh))
            out_assam$linear_predictor[,,k0] <- out_assam$linear_predictor[,,k0] + get_spatial_fields
        }
    rm(get_eta)
     
    out_assam$spatial_fields <- get_spatial_fields
    out_assam$logL <- do_em$new_logL
    out_assam$df <- num_spp*(num_nuisance_perspp + 1) + prod(dim(out_assam$betas)) + (num_archetypes - 1)
    out_assam$control <- control
    out_assam$bootstrap_control <- bootstrap_control
    gc()
    
    
    ##----------------
    #' # Standard Error using parametric bootstrap 
    #' Computation time not great on this at the moment unless you have a HPC!!! Starts from estimated parameters to give a little speed on for the quadratic approximations
    ##----------------
    if(uncertainty_quantification) {
        message("Performing parametric bootstrap to obtain uncertainty quantification...this will take a while so go a brew a cup of tea (or two)!")
        bootstrap_control$ci_type <- match.arg(bootstrap_control$ci_type, choices = c("percentile", "expanded")) 
        
        #' ## Bootstrap datasets -- This is *very* slow due to to sdmTMB_simulate
        class(out_assam) <- "assam"
        out_assam$sdmTMB_fits <- get_qa$sdmTMB_fits
        bootresp <- simulate.assam(out_assam,
                                   nsim = bootstrap_control$num_boot,
                                   do_parallel = TRUE,
                                   num_cores = num_cores,
                                   seed = bootstrap_control$seed)
        gc()
        
        #' ## Fit ASSAM to each bootstrapped dataset
        bootcov_fn <- function(b0, control) { 
            ##----------------
            #' ## Construct quadratic approximations for each species in bootstrap dataset, and set up relevant quantities
            ##----------------
            get_boot_qa <- try(.quadapprox2_fn(family = family, 
                                               formula = formula, 
                                               resp = bootresp[,b0]$y, 
                                               data = data, 
                                               mesh = mesh,
                                               offset = offset,
                                               trial_size = trial_size,
                                               do_parallel = do_parallel,
                                               return_fits = FALSE),
                               silent = TRUE)
            gc()
            if(inherits(get_boot_qa, "try-error"))
                return(get_boot_qa)
            
            get_boot_qa$long_parameters <- apply(get_boot_qa$parameters, 1, function(x) kronecker(rep(1,num_archetypes), x))
               
            ##----------------
            #' # Run EM algorithm for bootstrap datasets
            ##----------------
            do_em <- try(em_fn(qa_object = get_boot_qa), silent = TRUE)
            
            setTxtProgressBar(pb, b0)
            return(do_em)
            }          
        
        pb <- txtProgressBar(min = 0, max = bootstrap_control$num_boot, style = 3)
        bootrun <- lapply(1:bootstrap_control$num_boot, bootcov_fn, control = control) 
        close(pb) 
        rm(pb)

        find_errors <- which(sapply(bootrun, function(x) inherits(x, "try-error")))
        if(length(find_errors) > 0) {
            warning("Bootstrapped datasets ", find_errors, "had problems during fitting, and subsequently ignored...\nIf the number of datasets with fitting problems is large, they may point deeper issues with the asSAM")
            bootrun <- bootrun[-find_errors]
            }
        rm(find_errors)
        gc()
        
        
        #' ## Account for potential label-switching across the bootstrapped datasets, and also transform nuisance parameters as appropriate
        boot_params <- lapply(bootrun, function(x) cbind(x$new_betas, x$new_mixprop))
        boot_params <- abind::abind(boot_params, along = 0)
        switch_labels <- label.switching::pra(mcmc.pars = boot_params, 
                                              pivot = cbind(out_assam$betas, out_assam$mixture_proportion))$permutations
        bootrun <- lapply(1:length(bootrun), function(k0) {
            bootrun[[k0]]$new_mixprop <- bootrun[[k0]]$new_mixprop[switch_labels[k0,]]
            bootrun[[k0]]$new_betas <- bootrun[[k0]]$new_betas[switch_labels[k0,],]
            #bootrun[[k0]]$post_prob <- bootrun[[k0]]$post_prob[,switch_labels[k0,]] #' Not actually used later on so omit!

            #' A vector of all parameters ordered in the same way as the columns of the mapping matrix, plus the mixture proportions
            bootrun[[k0]]$boot_params <- c(bootrun[[k0]]$new_spp_intercept, 
                                           as.vector(t(bootrun[[k0]]$new_betas)), 
                                           as.vector(unlist(t(bootrun[[k0]]$new_nuisance))), # Note untransformed parameters passed here as this is what comes out of EM alg
                                           bootrun[[k0]]$new_mixprop)
            names(bootrun[[k0]]$boot_params) <- c(colnames(mapping_mat), paste0(names(out_assam$mixture_proportion), "_", "mixture_proportion"))
            
            bootrun[[k0]]$new_logL <- bootrun[[k0]]$post_prob <- bootrun[[k0]]$new_params <- NULL
            return(bootrun[[k0]])
            })
        rm(boot_params, switch_labels, bootresp)
        
        
        #' ## Form the intervals
        form_cis <- list()
        modified_alpha <- bootstrap_control$ci_alpha
        if(bootstrap_control$ci_type == "expanded")
            modified_alpha <- 2*pnorm(sqrt(num_spp/(num_spp-1)) * qt(bootstrap_control$ci_alpha/2, df = num_spp - 1))
        form_cis$spp_intercepts <- data.frame( t(apply(sapply(bootrun, function(x) x$new_spp_intercept), 1, quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)))
        form_cis$betas <- apply(abind::abind(lapply(bootrun, function(x) x$new_betas), along = 0), c(2,3), quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)
        form_cis$betas <- list(lower = form_cis$betas[1,,], upper = form_cis$betas[2,,])
        if(num_nuisance_perspp > 0) {
            form_cis$spp_nuisance <- apply(abind::abind(lapply(bootrun, function(x) x$new_transformed_nuisance), along = 0), c(2,3), quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)
            form_cis$spp_nuisance <- list(lower = form_cis$spp_nuisance[1,,,drop = FALSE], upper = form_cis$spp_nuisance[2,,,drop = FALSE])
            }
        form_cis$mixture_proportion <- data.frame( t(apply(sapply(bootrun, function(x) x$new_mixprop), 1, quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)))
        
        colnames(form_cis$spp_intercepts) <- colnames(form_cis$mixture_proportion) <- c("lower", "upper")
        rownames(form_cis$betas$lower) <- rownames(form_cis$betas$lower) <- rownames(out_assam$betas)
        if(num_nuisance_perspp > 0) 
            dimnames(form_cis$spp_nuisance$lower)[[3]] <- dimnames(form_cis$spp_nuisance$upper)[[3]] <- colnames(out_assam$spp_nuisance)
        
        out_assam$confidence_intervals <- form_cis
        out_assam$bootstrap_parameters <- t(sapply(bootrun, function(x) x$boot_params))
        
        
        #' ## Calculate bootstrapped posterior probabilities of *original species data* belong to each archetype
        .calc_posterior_prob <- function(cw_bootstrap_spp_intercept, cw_bootstrap_betas, cw_bootstrap_nuisance, cw_bootstrap_mixprop, qa_object) {
            post_prob <- matrix(NA, nrow = num_spp, ncol = num_archetypes)
            rownames(post_prob) <- colnames(y)
            colnames(post_prob) <- paste0("archetype", 1:num_archetypes)
            logL_spp <- NULL         
            
            for(j in 1:num_spp) { 
                cw_Quad <- sapply(1:num_archetypes, function(k) {
                    cw_v <- matrix(c(cw_bootstrap_spp_intercept[j], cw_bootstrap_betas[k,], cw_bootstrap_nuisance[j,]) - qa_object$parameters[j,], ncol = 1)
                    return(-0.5 * (crossprod(cw_v, qa_object$hessian[[j]]) %*% cw_v))
                    })
                eps <- max(cw_Quad)
                
                logL_spp[j] <- log(sum(cw_bootstrap_mixprop * exp(cw_Quad - eps))) + eps
                post_prob[j,] <- exp((log(cw_bootstrap_mixprop) + cw_Quad) - logL_spp[j])
                rm(eps, cw_Quad)
                }
            
            return(post_prob)
            }
        
        out_assam$bootstrap_posterior_probability <- abind::abind(lapply(1:length(bootrun), function(k0) 
            .calc_posterior_prob(cw_bootstrap_spp_intercept = bootrun[[k0]]$new_spp_intercept, 
                                 cw_bootstrap_betas = bootrun[[k0]]$new_betas, 
                                 cw_bootstrap_nuisance = bootrun[[k0]]$new_nuisance, 
                                 cw_bootstrap_mixprop = bootrun[[k0]]$new_mixprop, 
                                 qa_object = get_qa)),
            along = 3)
        } 
          
    ##----------------
    #' # Done! Final touches if any
    ##----------------
    class(out_assam) <- "assam"
    return(out_assam)
    }
     

# Hidden function to predict species-specific spatial fields from an assam fit. 
# An ad-hoc but computationally scalable approach is adopted where the field is predicted based on the most likely archetype the species is in.
#' @noMd
#' @noRd
.predict_spatial_fields <- function(l, mesh, qa_object, em_object, family) {
    if(is.null(mesh))
        return(NULL)
    
    if(!is.null(mesh)) {
        #' Recompile the TMB object from the original sdmTMB fits, evaluated at the assam parameter values
        #' Note an nlminb is required to update the environment parameter values 
        use_pars <- qa_object$sdmTMB_fits[[l]]$tmb_obj$env$parList()
        use_pars[["b_j"]] <- c(em_object$new_spp_intercept[l], em_object$new_betas[which.max(em_object$post_prob[l,]), ])
        if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
            use_pars[["ln_phi"]] <- em_object$new_nuisance[l, grep("ln_phi", colnames(em_object$new_nuisance))]
        if(family$family[1] == "tweedie") 
            use_pars[["thetaf"]] <- em_object$new_nuisance[l, grep("thetaf", colnames(em_object$new_nuisance))]
        use_pars[["ln_tau_O"]] <- em_object$new_nuisance[l, grep("ln_tau_O", colnames(em_object$new_nuisance))]
        use_pars[["ln_kappa"]] <- matrix(em_object$new_nuisance[l, grep("ln_kappa", colnames(em_object$new_nuisance))], nrow = 2, ncol = 1) # This has two rows as set up as sdmTMB
        
        #' Truncate spatial parameters that are very large in magnitude 
        if(use_pars[["ln_tau_O"]] < -30) use_pars[["ln_tau_O"]] <- -30
        if(use_pars[["ln_tau_O"]] > 30) use_pars[["ln_tau_O"]] <- 30
        if(any(use_pars[["ln_kappa"]] < -30)) use_pars[["ln_kappa"]] <- matrix(-30, nrow = 2, ncol = 1)
        if(any(use_pars[["ln_kappa"]] > 30)) use_pars[["ln_kappa"]] <- matrix(30, nrow = 2, ncol = 1)
        
        use_map <- list(b_j = as.factor(rep(NA, length(use_pars[["b_j"]])))) 
        if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie") > 0) 
            use_map$ln_phi <- as.factor(NA) 
        if(family$family[1] == "tweedie") 
            use_map$thetaf <- as.factor(NA)
        use_map$ln_tau_O <- as.factor(NA)
        use_map$ln_kappa <- as.factor(matrix(NA, 2, 1)) 

        new_tmb_obj <- TMB::MakeADFun(data = qa_object$sdmTMB_fits[[l]]$tmb_data,
                                      profile = qa_object$sdmTMB_fits[[l]]$control$profile,
                                      parameters = qa_object$sdmTMB_fits[[l]]$tmb_params,
                                      map = use_map,
                                      random = qa_object$sdmTMB_fits[[l]]$tmb_random,
                                      DLL = "sdmTMB",
                                      silent = TRUE)

        new_fit0 <- nlminb(start = new_tmb_obj$par,
                           objective = new_tmb_obj$fn,
                           gradient = new_tmb_obj$gr,
                           control = qa_object$sdmTMB_fits[[l]]$nlminb_control)
        

        #' Use TMB's report function to return the estimated field of assam
        r <- new_tmb_obj$report(new_tmb_obj$env$last.par.best)
        return(as.vector(r$omega_s_A))
        }
    }




