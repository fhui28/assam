#' @title Approximate and scalable species archetype models (asSAMs)
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Fits approximate and scalable species archetype modeling (asSAMs) for model-based clustering of species based on their environmental response, into a small number of so-called archetypal responses. The basic idea is take the log-likelihood function of SAM, and then construct an approximation of this which (hopefully) is more scalable in both the number of sites and species.
#' 
#' @param y A multivariate abundance response matrix.
#' @param formula An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). *Note there should be nothing on the left hand side of the "~".* Currently, smooth terms are not permitted.
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param which_spp_effects A vector identifying which columns of the model matrix induced by \code{formula} and \code{data} should be treated as species-specific effects. Default to 1, meaning only the first column i.e., the intercept, is species-specific.
#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param offset A matrix of offset terms, of the same dimension as \code{y}.
#' @param trial_size Trial sizes to use for binomial distribution. This should equal to a scalar.
#' @param num_archetypes Number of archetypes (clusters) to assume in the asSAM.
#' @param mesh Output from [sdmTMB::make_mesh()], used for adding species-specific spatial fields to the linear predictor.
#' @param do_parallel Should parallel computing be used to fit the asSAM. Defaults to \code{FALSE}.
#' @param num_cores If \code{do_parallel = TRUE}, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to \code{parallel::detectCores() - 2}.
#' @param beta_selection Should variable selection be performed on the archetypal regression coefficients via the broken adaptive ridge [BAR](https://www.sciencedirect.com/science/article/pii/S0047259X17305067) penalty? Defaults to \code{FALSE}.
#' @param uncertainty_quantification Should uncertainly intervals be computed via parametric bootstrap?
#' @param supply_quadapprox An object of class \code{assam_quadapprox}, which is (mostly likely) obtained as a consequence of running an initial fit with \code{do_assam_fit = FALSE}. Supplying this can be useful when multiple asSAMs e.g., with different values of \code{num_archetypes} are needed; see [https://github.com/fhui28/assam/issues/8](https://github.com/fhui28/assam/issues/8) for an example of its usage. 
#' @param do_assam_fit If \code{FALSE}, then the ingredients needed to construct the approximation to the log-likelihood are returned, *without fitting the asSAM itself*. This can be useful when multiple asSAMs e.g., with different values of \code{num_archetypes} are needed. Otherwise, this function should be kept at its default value of \code{TRUE}. 
#' @param control A list containing the following elements:
#' \describe{
#' \item{max_iter:}{the maximum number of iterations in the EM algorithm. Usually convergence is quite quick e.g., less than 20 iterations.}
#' \item{tol:}{the convergence criterion; the difference in the log-likelihood value of the asSAM from successive iterations must be smaller than this value.}
#' \item{temper_prob:}{in the iteration of the EM algorithm, posterior probabilities from the E-step are "tempered" or push away from the 0/1 boundary. This is often useful to get the EM algorithm moving initially.}
#' \item{trace:}{controls if messages are printed as part of the estimation process to reflect progress.}
#' \item{beta_lower:}{a vector that can be used to constrain the lower limit of the regression coefficients for each (and every) archetype, along with the species-specific effects. The length of this vector should be the same as the number of columns of the model matrix induced by \code{formula} and \code{data}. However, no checks are made on this vector, so *please ensure you get the length right to ensure the correct implementation!* Defaults to \code{NULL}, in which there is no constraint.}
#' \item{beta_upper:}{a vector that can be used to constrain the upper limit of the regression coefficients for each (and every) archetype, along with the species-specific effects. The length of this vector should be the same as the number of columns of the model matrix induced by \code{formula} and \code{data}. However, no checks are made on this vector, so *please ensure you get the length right to ensure the correct implementation!* Defaults to \code{NULL}, in which there is no constraint.}
#' }
#' @param beta_selection_control A list containing the following elements to control the broken adaptive ridge (BAR) penalty for variable selection on the archetypal regression coefficients:
#' \describe{
#' \item{lambda:}{the tuning parameter for the BAR penalty. Note the function only accepts a single value for this; if you wish to construct a regularization path for the archetypal regression coefficients, then please consider using the [passam()] function instead.}
#' \item{warm_start:}{a list containing a set of values to "warm start" the BAR optimization part of the EM algorithm. This can be useful to speed up the optimization process, but also in some cases can help to reduce inconsistencies between the penalized estimates obtained here versus as part of the regularization path using [passam()]. The list must contain the following elements: 1) \code{beta}, a matrix of the archetypal regression coefficients; 2) \code{spp_effects}, a vector or matrix of species-specific effects; 3) \code{spp_nuisance}, a vector or matrix of species-specific nuisance parameters; 4) \code{mixture_proportion}, a vector of mixture proportions; 5) \code{posterior_probability}, a matrix of posterior probabilities for each species belonging to each archetype.}. 
#' \item{max_iter:}{the maximum number of iterations in the BAR optimization part of the EM algorithm.}
#' \item{eps:}{the convergence criterion; the norm of the difference between all estimated parameters from successive iterations must be smaller than this value.}
#' \item{round_eps:}{a tolerance to round values to zero. The technically not needed as the BAR penalty will produce exactly zero estimates up to machine error, but is included anyway, but is included anyway.}
#' }
#' @param bootstrap_control A list containing the following elements to control the parametric bootstrap for uncertainty quantification:
#' \describe{
#' \item{method:}{method of parametric bootstrap to use. Two options are currently available: 1) "full_bootstrap", which is a full parametric bootstrap where new multivariate abundance responses are simulated and a new approximate likelihood function is formed each time; 2) "fast_bootstrap" which bootstraps directly/only off the approximate likelihood. 
#' Method 1 should be more accurate, but is computationally slower and less scalable. Defaults to "full_bootstrap".}
#' \item{num_boot:}{the number of bootstrapped iterations to do. Defaults to 100, which can already take a long time but should be enough in a lot of settings for uncertainty quantification.}
#' \item{ci_alpha:}{the type-1 level for confidence interval construction. \code{100 * (1 - ci_alpha)} percent Confidence intervals are constructed.}
#' \item{seed:}{a seed that can be set for bootstrapping the datasets.}
#' \item{ci_type:}{type of confidence intervals to construct. Two options are currently available: 1) "percentile" confidence intervals based directly on the empirical quantiles of the bootstrap samples; 2) "expanded" percentile confidence intervals, which are typically slightly wider intervals that attempt to correct for a so-called "narrowness bias". Defaults to "percentile".}
#' }

#' @details 
#' For the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = u_i^\top\alpha_j + x_i^\top\beta_k,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{u_i^\top\alpha_j} corresponds to a component that is to kept species-specific e.g., species-specific intercept, \eqn{x_i^\top\beta_k}  denotes the component corresponding to effect of archetypal response \eqn{k}. Additionally, species-specific spatial fields can be included in the linear predictor e.g., to account for residual spatial correlation above and beyond that explained by the archetypal responses. Conditional on the mean model above, the \eqn{y_{ij}} are assumed to be come from some response distribution using the additional dispersion and power parameters as appropriate. We refer the reader to [Dunstan et al., (2011)](https://doi.org/10.1016/j.ecolmodel.2010.11.030), [Hui et al., (2013)](https://doi.org/10.1890/12-1322.1), [Dunstan et al., (2013)](https://doi.org/10.1007/s13253-013-0146-x), and [Skipton Woolley's ecomix package](https://github.com/skiptoniam/ecomix) for more details about the formulations of SAMs. 
#' 
#' The broad goal of this package is to construct a way of fitting SAMs that are, although approximate, more scalable in the number of sites and species (though not necessarily faster), hence the name asSAMs. We refer to the corresponding manuscript (in preparation) for details, but to summarize, asSAMs are formed by constructing an approximate likelihood function for a SAM based on using ingredients (i.e., point estimates and the associated observed information matrix) from stacked species distribution models (which are fitted initially in parallel), and then building what is essentially a finite mixture of multivariate Gaussian distributions from this. This is then maximized using an EM algorithm, which can be done scalably and very quickly.
#' 
#' For uncertainty quantification, two forms of parametric bootstrap are available along the lines of Section 2.16.2 in [McLachlan and Peel (2004)](https://www.wiley.com/en-us/Finite+Mixture+Models-p-9780471654063). 
#' 
#' Common lower and upper limit constraints can be placed on the \eqn{\beta_k} for each (and every) archetype. 
#' 
#' Variable selection on the elements of \eqn{\beta_k} can be performed via the broken adaptive ridge [BAR](https://www.sciencedirect.com/science/article/pii/S0047259X17305067) penalty, via the \code{beta_selection} and \code{beta_selection_control} arguments. The BAR penalty can be interpreted as a kind of approximation to the \eqn{L_0} penalty, and encourages sparsity in the archetypal regression coefficients e.g., to uncover what covariates are informative for each of the archetypal responses. Note the function is set up to only accept a single value for the tuning parameter: to construct a full regularization path for the archetypal regression coefficients given the BAR penalty and/or to perform selection on the mixture proportions and thus select the number of archetypes to include in a asSAM, please consider using the [passam()] function instead. Note, in the case model selection is performed, bootstrap uncertainty quantification is performed *conditional on selected model*. 

#' 
#' 
#' \subsection{Distributions}{
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
#' The scalability of asSAMs relies on being able to deploy [sdmTMB::sdmTMB()] is an efficient and fully optimized manner. Along these lines, please see [using sdmTMB in parallel](https://github.com/pbs-assess/sdmTMB/issues/368) and [using OpenBLAS](https://gist.github.com/seananderson/08a51e296a854f227a908ddd365fb9c1) and references therein for some tips to ensure things on your machine are more optimized. Thanks to Sean Anderson for this advice!
#' 
#' 
#' @return If \code{do_assam_fit = FALSE}, then a object of class \code{assam_quadapprox} is returned, which contains the ingredients needed to fit asSAMs. 
#' 
#' Otherwise, if \code{do_assam_fit = TRUE} as is the default, then an object of class \code{assam} with the following elements (as appropriate, and not necessarily in the order below):
#' \item{call:}{The function call.}
#' \item{formula:}{Same as input argument.}
#' \item{family:}{Same as input argument.}
#' \item{num_nuisance_perspp:}{The number of "nuisance" parameters per species. For example, if \code{family = binomial()} and species-specific spatial fields are not included, then there are no nuisance parameters. If \code{family = nbinom2()} and species-specific spatial fields are included, say, then there are three nuisance parameters per species.}
#' \item{trial_size:}{Same as input argument.}
#' \item{offset:}{Same as input argument.}
#' \item{add_spatial:}{Were species-specific spatial fields included?}
#' \item{mesh:}{Same as input argument.}
#' \item{num_archetypes:}{Same as input argument.}
#' \item{uncertainty_quantification:}{Same as input argument.}
#' \item{spp_effects:}{Estimated species-specific effects i.e., \eqn{alpha_j}.}
#' \item{betas:}{Estimated matrix of archetypal regression coefficients corresponding to the model matrix created i.e., \eqn{beta_k}. The number of rows in \code{betas} is equal to the number of archetypes.}
#' \item{mixture_proportion:}{Estimated vector of mixture proportions corresponding to the probability of belonging to each archetype.}
#' \item{spp_nuisance:}{Estimated matrix of species-specific nuisance parameters e.g., the dispersion parameter in the negative binomial distribution, and dispersion and power parameters in the Tweedie distribution, and so on.}
#' \item{posterior_probability:}{Estimated matrix of posterior probabilities for each species belong to each archetype. The number of rows in \code{posterior_probability} is equal to the number of species.}
#' \item{linear_predictor:}{Estimated array of archetype-specific linear predictors for all species. The last dimension of the array corresponds to the number of archetypes.}
#' \item{spatial_fields:}{Predicted matrix species-specific spatial fields, if included. The spatial field prediction is developed based on the most likely archetype that the species belongs to, as judged by the posterior probabilities.}
#' \item{logL:}{Estimated log-likelihood value of the asSAM i.e. the value of the approximated log-likelihood function at convergence.}
#' \item{df:}{Number of the estimated (freely varying) parameters in the asSAM.}
#' \item{linear_predictor:}{Estimated array of archetype-specific linear predictors. The last dimension of the array corresponds to the number of archetypes.}
#' \item{control:}{Same as input argument.}
#' \item{bootstrap_control:}{Same as input argument.}
#' \item{beta_selection_control:}{Same as input argument.}
#' \item{confidence_intervals:}{A list of estimated parametric bootstrap confidence intervals for all the parameters in the asSAM, with each element in the list corresponding to one of the parameters e.g., the species-specific effects, the mixture proportions, and so on.}
#' \item{bootstrap_paramters:}{A matrix of estimated (untransformed) parameters from the parametric bootstrap. *This output can be safely ignored by most users*, but those curious, it is used to uncertainty quantification for downstream predictions, say.}
#' \item{bootstrap_posterior_probability:}{An array of estimated posterior probabilities from the parametric bootstrap. *This output can be safely ignored by most users*, but those curious, it is used to uncertainty quantification for downstream predictions, say.}
#' \item{sdmTMB_fits}{A list of the set of stacked species distribution models fitted using [sdmTMB::sdmTMB()]. *This output can be safely ignored by most users*.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @examples
#' \dontrun{
#' ##----------------------
#' # Example 1: Generate some multivariate abundance (count) data from a SAM
#' ##----------------------
#' library(tidyverse)
#' library(mvtnorm)
#' library(GGally)
#' library(doParallel)
#' 
#' set.seed(022025)
#' 
#' num_X <- 10
#' num_units <- 1000
#' num_spp <- 100
#' num_archetype <- 5
#' H <- outer(1:num_X, 1:num_X, "-")
#' H <- 0.5^abs(H)
#' covariate_dat <- rmvnorm(num_units, sigma = H) %>% 
#'     as.data.frame %>% 
#'     rename_with(., .fn = function(x) paste0("covariate", x))
#' rm(H)
#' 
#' true_betas <- runif(num_archetype * num_X, -1, 1) %>% matrix(nrow = num_archetype)
#' true_spp_effects <- matrix(runif(num_spp, -2, 0), ncol = 1)
#' true_dispparam <- 1/runif(num_spp, 1, 5) 
#' true_powerparam <- runif(num_spp, 1.4, 1.8)
#' true_mixprop <- c(0.2, 0.2, 0.3, 0.15, 0.15)
#'  
#' simdat <- create_samlife(family = nbinom2(), 
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula, 
#' data = covariate_dat, 
#' betas = true_betas, 
#' spp_effects = true_spp_effects, 
#' spp_dispparam = true_dispparam, 
#' spp_powerparam = true_powerparam, 
#' mixture_proportion = true_mixprop,
#' seed = 022025)
#'
#'  
#' ## Fit asSAM and assess results 
#' #' **Most users should start here**
#' samfit <- assam(y = simdat$y,
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' family = nbinom2(),
#' num_archetypes = num_archetype,
#' num_cores = detectCores() - 2)
#' 
#'  
#' plot(true_spp_effects, samfit$spp_effects); abline(0,1)
#' plot(true_dispparam, samfit$spp_nuisance$dispersion, log = "xy"); abline(0,1)
#' #' Note estimates for the archetypal responses and mixture proportions from (as)SAMs should be 
#' #' close to the corresponding true values, *up to a reordering* of the mixture component
#' #' s/archetypes (since the order is essentially arbitrary)
#' rbind(true_betas, samfit$betas) %>% 
#' t %>% 
#' as.data.frame %>%
#' GGally::ggpairs(.)
#' table(simdat$archetype_label, apply(samfit$posterior_probability, 1, which.max))
#' 
#'  
#' ## Demonstrating basic use of functions for asSAM 
#' samfit
#' summary(samfit)
#' 
#' fitted(samfit)
#'  
#' simulate(samfit, data = covariate_dat)
#' 
#' residuals(samfit, type = "dunnsmyth")
#'  
#' #' Basic residual analysis
#' plot(samfit, transform_fitted_values = TRUE, envelope = FALSE)
#'  
#' #' Archetype-level predictions
#' predict(samfit, newdata = covariate_dat, type = "archetype", se_fit = TRUE) 
#' 
#' #' Species-level predictions
#' predict(samfit, newdata = covariate_dat, type = "species_max", 
#' num_cores = detectCores() - 2, se_fit = TRUE) 
#'  
#'  
#'  
#' ##----------------------
#' # Example 2: Demonstrating variable selection on the archetypal regression coefficients
#' # Generate some multivariate abundance (non-negative continuous) data from a sparse SAM
#' # Note only a single tuning parameter is used below; please see the [passam()] for 
#' # constructing a proper regularization path, as well as to perform selection on the 
#' # mixing proportions i.e., choose the number of archetypes. 
#' ##----------------------
#' true_betas <- runif(num_archetype * num_X, -1, 1) %>% matrix(nrow = num_archetype)
#' true_betas[which(abs(true_betas) < 0.4)] <- 0 # Making archetypal coefficients sparse
#' true_betas
#' 
#' 
#' simdat <- create_samlife(family = tweedie(),
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' betas = true_betas,
#' spp_effects = true_spp_effects,
#' spp_dispparam = true_dispparam,
#' spp_powerparam = true_powerparam,
#' mixture_proportion = true_mixprop,
#' seed = 022025)
#' 
#' 
#' ## Fit asSAM and assess results 
#' #' **Most users should start here**
#' samfit_select <- assam(y = simdat$y,
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' family = tweedie(),
#' beta_selection = TRUE,
#' num_archetypes = num_archetype,
#' beta_selection_control = list(lambda = 0.1), # Note this an arbitrary choice!
#' bootstrap_control = list(num_boot = 10), 
#' num_cores = detectCores() - 2)
#' 
#' samfit_select
#' samfit_select$betas
#' true_betas
#' 
#' 
#' plot(true_spp_effects, samfit_select$spp_effects); abline(0,1)
#' plot(true_dispparam, samfit_select$spp_nuisance$dispersion, log = "xy"); abline(0,1)
#' plot(true_powerparam, samfit_select$spp_nuisance$power, log = "xy"); abline(0,1)
#' #' Note estimates for the archetypal responses and mixture proportions from (as)SAMs should be 
#' #' close to the corresponding true values, *up to a reordering* of the mixture component
#' #' s/archetypes (since the order is essentially arbitrary)
#' rbind(true_betas, samfit_select$betas) %>% 
#' t %>% 
#' as.data.frame %>%
#' GGally::ggpairs(.)
#' table(simdat$archetype_label, apply(samfit_select$posterior_probability, 1, which.max))
#' 
#' 
#' ## Demonstrating basic use of functions for asSAM 
#' summary(samfit_select)
#' 
#' fitted(samfit_select)
#'  
#' simulate(samfit_select, data = covariate_dat)
#' 
#' residuals(samfit, type = "dunnsmyth")
#'  
#' #' Basic residual analysis
#' plot(samfit_select, transform_fitted_values = TRUE, envelope = FALSE)
#'  
#' #' Archetype-level predictions
#' predict(samfit_select, newdata = covariate_dat, type = "archetype") 
#' 
#' #' Species-level predictions
#' predict(samfit_select, newdata = covariate_dat, type = "species_max") 
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
#' @importFrom sdmTMB sdmTMB sdmTMBcontrol make_mesh nbinom2 tweedie Beta
#' @importFrom label.switching pra
#' @importFrom methods as
#' @importFrom parallel detectCores
#' @importFrom quadprog solve.QP
#' @importFrom stats as.formula nlminb model.matrix qt plogis qlogis 
#' @importFrom TMB MakeADFun
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import Matrix
#' @md

assam <- function(y, 
                  formula, 
                  data, 
                  which_spp_effects = 1,
                  family, 
                  offset = NULL, 
                  trial_size = 1, 
                  num_archetypes,
                  mesh = NULL, 
                  do_parallel = TRUE, 
                  num_cores = NULL, 
                  beta_selection = FALSE,
                  uncertainty_quantification = FALSE,
                  supply_quadapprox = NULL,
                  do_assam_fit = TRUE,
                  control = list(max_iter = 500, tol = 1e-4, 
                                 temper_prob = 0.8, trace = FALSE, 
                                 beta_lower = NULL, beta_upper = NULL),
                  beta_selection_control = list(lambda = 1, warm_start = NULL, 
                                                max_iter = 100, eps = 1e-4, round_eps = 1e-5),
                  bootstrap_control = list(method = "full_bootstrap", 
                                           num_boot = 100, ci_alpha = 0.05, seed = NULL, 
                                           ci_type = "percentile")) {
    
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
    if(!is.null(num_cores)) {
        registerDoParallel(cores = num_cores)
        }
     
    .check_offset(offset = offset, y = y) 
    if(is.null(offset))
        offset <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    formula <- .check_X_formula(formula = formula, data = as.data.frame(data))     
    tmp_formula <- as.formula(paste("response", paste(as.character(formula),collapse = " ") ) )
    nullfit <- sdmTMB(tmp_formula,
                      spatial = FALSE,
                      data = data.frame(data, response = rnorm(nrow(data))))
    X <- model.matrix(nullfit$formula[[1]], data = nullfit$data)
    num_X <- ncol(X)
    rm(tmp_formula, nullfit) 
    
    if(is.null(mesh)) {
        add_spatial <- FALSE
        final_mesh <- NULL
        }
    if(!is.null(mesh)) {
        add_spatial <- TRUE
        if(!inherits(mesh, "sdmTMBmesh"))
            stop("If mesh is supplied for species-specific spatial fields, then the mesh argument must be an object class of \"sdmTMBmesh\".")
        final_mesh <- mesh
        }
    
    control <- .fill_control(control = control)
    beta_selection_control <- .fill_beta_selection_control(control = beta_selection_control)
    bootstrap_control <- .fill_bootstrap_control(control = bootstrap_control)
    bootstrap_control$method <- match.arg(bootstrap_control$method, choices = c("full_bootstrap", "fast_bootstrap"))
    bootstrap_control$ci_type <- match.arg(bootstrap_control$ci_type, choices = c("percentile", "expanded"))
    
    .check_beta_options(control = control, 
                        beta_selection_control = beta_selection_control, 
                        beta_selection = beta_selection, 
                        uncertainty_quantification = uncertainty_quantification)
    
    num_unit <- nrow(y)
    num_spp <- ncol(y)
          

    ##----------------
    # Construct quadratic approximations for each species, and set up relevant quantities
    ##----------------
    message("Commencing fitting...")
    if(!is.null(supply_quadapprox)) {
        if(!inherits(supply_quadapprox, "assam_quadapprox"))
            stop("If supply_quadapprox is supplied, then it must be an object of class \"assam_quadapprox\".")
        
        message("Species-specific quadratic approximations supplied by practitioner...thanks!")
        
        get_qa <- supply_quadapprox 
        }
    
    if(is.null(supply_quadapprox)) {
        if(control$trace)
            message("Forming species-specific quadratic approximations...")
        
        get_qa <- .quadapprox2_fn(family = family,
                                  formula = formula, 
                                  resp = y, 
                                  data = data,
                                  add_spatial = add_spatial,
                                  mesh = final_mesh,
                                  offset = offset,
                                  trial_size = trial_size,
                                  do_parallel = do_parallel,
                                  control = control) 
        }
    
    if(!do_assam_fit) {
        message("Returning species-specific quadratic approximations; no assam fitted...have a lovely day!")
        return(get_qa)
        }   
    
    ### Make mapping matrix, maps psi to long_parameters
    #' Psi sequence: species-specific effects; archetypal regression coefficients by archetype; species-specific nuisance parameters, along with parameters for species-specific spatial fields
    get_qa$long_parameters <- apply(get_qa$parameters, 1, function(x) kronecker(rep(1,num_archetypes), x)) # Repeats species-specific estimates num_archetypes times; object has num_spp columns
    num_nuisance_perspp <- length(get_qa$parameters[1,]) - num_X # e.g., number of dispersion and power parameter per-species, along with parameters for species-specific spatial fields
    qa_parameters_colnames <- colnames(get_qa$parameters)
    
    mapping_mat <- Matrix(0, nrow = length(get_qa$long_parameters), ncol = num_spp * length(which_spp_effects) + (num_X - length(which_spp_effects)) * num_archetypes + num_spp*num_nuisance_perspp)
    makecolnames <- c(paste0(rep(paste0("spp_effects", 1:num_spp), each = length(which_spp_effects)), rep(paste0("_alpha", 1:length(which_spp_effects)), num_spp)), 
                      paste0(rep(paste0("archetype", 1:num_archetypes), each = num_X - length(which_spp_effects)), rep(paste0("_beta", 1:(num_X - length(which_spp_effects))), num_archetypes)))
    if(num_nuisance_perspp > 0)
        makecolnames <- c(makecolnames,
                          paste0(rep(paste0("spp_nuisance", 1:num_spp), each = num_nuisance_perspp), rep(colnames(get_qa$parameters)[num_X + 1:num_nuisance_perspp], num_spp))
                          )
    colnames(mapping_mat) <- makecolnames
    rm(makecolnames)
    
    makerownames_mappingmat_fn <- function(k, j) {
        out <- character(length = num_X)
        out[which_spp_effects] <- paste0("spp_effects", j, "_alpha", 1:length(which_spp_effects))
        out[-which_spp_effects] <- paste0("archetype", k, "_beta", 1:(num_X - length(which_spp_effects)))
        if(num_nuisance_perspp > 0)
            out <- c(out, paste0("spp_nuisance", j, colnames(get_qa$parameters)[num_X + 1:num_nuisance_perspp]))
        
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
    em_fn <- function(qa_object, 
                      dobar_penalty,
                      warm_start = NULL,
                      betamatrix_selection = NULL) {
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
                 if(is.null(warm_start)) {
                     do_kmeans <- cluster::pam(qa_object$parameters[, (1:num_X)[-which_spp_effects], drop = FALSE], k = num_archetypes)  
                     if(is.null(betamatrix_selection)) {
                         cw_betas <- do_kmeans$medoids
                         cw_mixprop <- as.vector(table(do_kmeans$clustering)) / num_spp
                         }
                     if(!is.null(betamatrix_selection)) {
                         #' Reordering to try and get to a starting value that best respects the desired sparsity pattern. No guarantees this actually works that well in practice, especially if a lot of non-zero coefficients are weak!
                         switch_labels <- label.switching::pra(mcmc.pars = abind::abind(do_kmeans$medoids, along = 0), 
                                                             pivot = betamatrix_selection)$permutations
                         cw_betas <- do_kmeans$medoids[switch_labels[1,], ]
                         cw_mixprop <- as.vector(table(do_kmeans$clustering))[switch_labels[1,]] / num_spp
                         rm(switch_labels)
                         }
                     cw_spp_effects <- qa_object$parameters[, which_spp_effects, drop = FALSE]
                     cw_nuisance <- NULL
                     if(num_nuisance_perspp > 0)
                        cw_nuisance <- qa_object$parameters[, num_X + 1:num_nuisance_perspp, drop = FALSE]
                     rm(do_kmeans)
                     }
                 if(!is.null(warm_start)) {
                     cw_betas <- warm_start$betas
                     cw_spp_effects <- warm_start$spp_effects
                     cw_nuisance <- warm_start$spp_nuisance
                     cw_mixprop <- warm_start$mixing_proportion
                     }
                 } 
             
             
             ##-------------------
             #' ## E-step
             #' Ignore the normalizing constant of the normal distribution as that does not vary as a function of archetype anyway
             ##-------------------
             for(j in 1:num_spp) { 
                 cw_Quad <- sapply(1:num_archetypes, function(k) {
                     cw_params <- c(cw_spp_effects[j,], cw_betas[k,], cw_nuisance[j,])
                     cw_params[which_spp_effects] <- cw_spp_effects[j,]
                     cw_params[(1:num_X)[-which_spp_effects]] <- cw_betas[k,]
                     
                     cw_v <- matrix(cw_params - qa_object$parameters[j,], ncol = 1)
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
                 if(num_archetypes > 1) {
                     ## Temper the classification in the initial E-step if num_archetypes > 1
                     alpha_temper <- (1 - control$temper_prob * num_archetypes) / (control$temper_prob * (2-num_archetypes) - 1)          
                     for(j in 1:num_spp)
                         post_prob[j,] <- (2*alpha_temper*post_prob[j,]-alpha_temper+1)/(2*alpha_temper - alpha_temper*num_archetypes + num_archetypes)
                     new_logL <- -1e8               
                     }
                 if(!is.null(warm_start) & dobar_penalty) {
                     post_prob <- warm_start$posterior_probability
                     #new_logL <- warm_start$new_logL
                    }
                }
                    
             ##-------------------
             #' ## M-step 
             ##-------------------
             new_mixprop <- colMeans(post_prob)
             bigW <- Diagonal(x = rep(as.vector(t(sqrt(post_prob))), each = nrow(qa_object$hessian[[j]]))) %*% basic_bigW %*% Diagonal(x = rep(as.vector(t(sqrt(post_prob))), each = nrow(qa_object$hessian[[j]])))
             MtWM <- forceSymmetric(crossprod(mapping_mat, bigW) %*% mapping_mat)
             new_params <- solve(MtWM, crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters))
             new_params <- as.vector(new_params)
             
             
             #' ### For selection on the archetypal coefficients via the broken adaptive ridge (BAR) penalty 
             if(dobar_penalty) {
                 beta_s_err <- Inf
                 beta_s_counter <- 0
                 cw_params <- new_params
                 Dbar <- Diagonal(n = length(new_params))
                 diag(Dbar)[-grep("archetype", colnames(mapping_mat))] <- 0
                 lambda_used <- num_spp *  beta_selection_control$lambda 
                 
                 while(beta_s_counter < beta_selection_control$max_iter & beta_s_err > beta_selection_control$eps) {
                     GammaMatrix_params <- Diagonal(n = length(new_params))
                     diag(GammaMatrix_params)[grep("archetype", colnames(mapping_mat))] <- cw_params[grep("archetype", colnames(mapping_mat))]
                     
                     new_params <- GammaMatrix_params %*% solve(GammaMatrix_params %*% MtWM %*% GammaMatrix_params + lambda_used * Dbar) %*% GammaMatrix_params %*% (crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters))
                     new_params <- as.vector(new_params)
                     names(new_params) <- colnames(mapping_mat)
                     
                     beta_s_err <- sum((new_params - cw_params)^2)
                     cw_params <- new_params
                     beta_s_counter <- beta_s_counter + 1
                     }
                 #' Note the convergence properties of the algorithm does ensure exact zeros can be achieved (up to machine error); see [https://doi.org/10.1155/2016/3456153]. However once can also round to zero if interested.
                 new_params[abs(new_params) < beta_selection_control$round_eps] <- 0
                 rm(beta_s_counter, beta_s_err, GammaMatrix_params, Dbar)
                 }
             
             
             #' ### For constrained solution
             if(!is.null(control$beta_lower) | !is.null(control$beta_upper)) {
                 if(!is.null(control$beta_lower)) {
                     makeAmat_lower <- diag(x = 0, nrow = ncol(mapping_mat))
                     rownames(makeAmat_lower) <- colnames(makeAmat_lower) <- colnames(mapping_mat)
                     diag(makeAmat_lower)[grep("spp_effects", rownames(makeAmat_lower))] <- 1
                     diag(makeAmat_lower)[grep("_beta", rownames(makeAmat_lower))] <- 1
                        
                     makebvec_lower <- numeric(ncol(mapping_mat))
                     names(makebvec_lower) <- colnames(mapping_mat)
                     makebvec_lower[grep("spp_effects", names(makebvec_lower))] <- control$beta_lower[which_spp_effects]
                     makebvec_lower[grep("_beta", names(makebvec_lower))] <- rep(control$beta_lower[-which_spp_effects], num_archetypes)
                     
                     final_makeAmat <- makeAmat_lower
                     final_makebvec <- makebvec_lower
                     }
                 if(!is.null(control$beta_upper)) {
                     makeAmat_upper <- diag(x = 0, nrow = ncol(mapping_mat))
                     rownames(makeAmat_upper) <- colnames(makeAmat_upper) <- colnames(mapping_mat)
                     diag(makeAmat_upper)[grep("spp_effects", rownames(makeAmat_upper))] <- -1
                     diag(makeAmat_upper)[grep("_beta", rownames(makeAmat_upper))] <- -1
                     
                     makebvec_upper <- numeric(ncol(mapping_mat))
                     names(makebvec_upper) <- colnames(mapping_mat)
                     makebvec_upper[grep("spp_effects", names(makebvec_upper))] <- control$beta_upper[which_spp_effects]
                     makebvec_upper[grep("_beta", names(makebvec_upper))] <- rep(control$beta_upper[-which_spp_effects], num_archetypes)
                     
                     final_makeAmat <- makeAmat_upper
                     final_makebvec <- makebvec_upper
                     }
                 if(!is.null(control$beta_lower) & !is.null(control$beta_upper)) {
                     final_makeAmat <- t(rbind(makeAmat_lower, makeAmat_upper))
                     final_makebvec <- c(makebvec_lower, makebvec_upper)
                     
                     rm(makeAmat_lower, makeAmat_upper, makebvec_lower, makebvec_upper)
                     }                 
                 
                 
                 #' Refine Amat and bvec arguments to things which only have constraints
                 remove_unconstrained_elements <- which(diag(final_makeAmat) == 0) 
                 if(length(remove_unconstrained_elements) > 0) {
                     final_makeAmat <- final_makeAmat[-remove_unconstrained_elements,]
                     final_makebvec <- final_makebvec[-remove_unconstrained_elements]
                    }
                 remove_unconstrained_elements <- which(!is.finite(final_makebvec))
                 if(length(remove_unconstrained_elements) > 0) {
                     final_makeAmat <- final_makeAmat[-remove_unconstrained_elements,]
                     final_makebvec <- final_makebvec[-remove_unconstrained_elements]
                     }
                 rm(remove_unconstrained_elements)
                 
                 Mstep_update <- quadprog::solve.QP(Dmat = MtWM, 
                                                    dvec = crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters),
                                                    Amat = t(final_makeAmat),
                                                    bvec = final_makebvec)
                 
                 new_params <- Mstep_update$solution
                 new_params[abs(new_params) < .Machine$double.eps] <- 0
                 rm(Mstep_update, final_makeAmat, final_makebvec)
                 }
             
             #' ### For subset solution
             if(!is.null(betamatrix_selection)) {
                 A <- Matrix(0, nrow = length(new_params), ncol = sum(betamatrix_selection == 0))
                 find_zero_indices <- matrix(grep("archetype", colnames(mapping_mat)), nrow = num_archetypes, byrow = TRUE) * (betamatrix_selection == 0)
                 find_zero_indices <- t(find_zero_indices)[t(find_zero_indices) != 0]
                 for(l0 in 1:ncol(A))
                     A[find_zero_indices[l0], l0] <- 1
                 rm(find_zero_indices)

                 Mstep_update <- quadprog::solve.QP(Dmat = MtWM,
                                                    dvec = crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters),
                                                    Amat = A,
                                                    meq = ncol(A),
                                                    factorized = FALSE)
                 new_params <- Mstep_update$solution
                 new_params[abs(new_params) < .Machine$double.eps] <- 0
                 rm(Mstep_update, A)
                 }

             names(new_params) <- colnames(mapping_mat)
             new_spp_effects <- matrix(new_params[grep("spp_effects", names(new_params))], nrow = num_spp, byrow = TRUE)
             new_betas <- matrix(new_params[grep("archetype", names(new_params))], nrow = num_archetypes, byrow = TRUE)
             new_nuisance <- NULL
             if(num_nuisance_perspp > 0) {
                 new_nuisance <- matrix(new_params[grep("spp_nuisance", names(new_params))], nrow = num_spp, byrow = TRUE)
                 colnames(new_nuisance) <- qa_parameters_colnames[num_X + 1:num_nuisance_perspp]  
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
             cw_spp_effects <- new_spp_effects
             cw_betas <- new_betas
             cw_nuisance <- new_nuisance
             cw_params <- new_params
             counter <- counter + 1          
             }
         
         rownames(new_spp_effects) <- colnames(y)
         rownames(new_betas) <- names(new_mixprop)
         colnames(new_spp_effects) <- colnames(X)[which_spp_effects]
         colnames(new_betas) <- colnames(X)[-which_spp_effects]
         
         tmp_nuisance <- NULL
         if(num_nuisance_perspp > 0) {
             if(family$family[1] %in% c("Beta", "gaussian", "Gamma", "nbinom2", "tweedie")) 
                 tmp_nuisance$dispersion <- exp(new_nuisance[,grep("ln_phi", colnames(new_nuisance))])
             if(family$family[1] == "tweedie")
                 tmp_nuisance$power <- plogis(new_nuisance[,grep("thetaf", colnames(new_nuisance))]) + 1
             if(add_spatial) { 
                 tmp_nuisance$spatial_range <- 1/exp(new_nuisance[,grep("ln_kappa", colnames(new_nuisance))])
                 est_tau <- exp(new_nuisance[,grep("ln_tau_O", colnames(new_nuisance))])
                 tmp_nuisance$spatial_SD <- 1/(sqrt(4*pi) * est_tau * exp(new_nuisance[,grep("ln_kappa", colnames(new_nuisance))]))
                 }
                 
             tmp_nuisance <- as.data.frame(tmp_nuisance)
             rownames(new_nuisance) <- rownames(tmp_nuisance) <- colnames(y)            
             }
         
         return(list(new_logL = new_logL, 
                     new_mixprop = new_mixprop, 
                     new_spp_effects = new_spp_effects, 
                     new_betas = new_betas, 
                     new_nuisance = new_nuisance, 
                     new_transformed_nuisance = tmp_nuisance,
                     new_params = new_params, 
                     counter = counter, 
                     post_prob = post_prob))          
         }

    if(control$trace)
        message("Commencing EM algorithm...")
    make_warm_start <- NULL
    if(!is.null(beta_selection_control$warm_start)) {
        make_warm_start <- beta_selection_control$warm_start
        if(!is.null(make_warm_start$spp_effects))
            make_warm_start$spp_effects <- matrix(make_warm_start$spp_effects, nrow = num_spp)
        if(!is.null(make_warm_start$spp_nuisance))
            make_warm_start$spp_nuisance <- matrix(make_warm_start$spp_nuisance, nrow = num_spp)
        make_warm_start$posterior_probability <- matrix(make_warm_start$posterior_probability, nrow = num_spp)
        }
    do_em <- em_fn(qa_object = get_qa,
                   dobar_penalty = beta_selection,
                   warm_start = make_warm_start)

    
    try_counter <- 0
    while(any(do_em$new_mixprop < 1e-3) & try_counter < 20) {
        message("Mixture component is being emptied...altering initial temp probability and restarting EM-algorithm to try and fix this.")
        control$temper_prob <- control$temper_prob + 0.025
          
        do_em <- em_fn(qa_object = get_qa,
                       dobar_penalty = beta_selection)
        try_counter <- try_counter + 1
        }

 
    ##----------------
    #' # Format output
    ##----------------
    message("Fitting completed, applying finishing touches...")
    
    out_assam <- list(call = match.call(), 
                      formula = formula,
                      family = family,
                      which_spp_effects = which_spp_effects,
                      num_nuisance_perspp = num_nuisance_perspp,
                      trial_size = trial_size,
                      offset = as(offset, "sparseMatrix"),
                      beta_selection = beta_selection,
                      add_spatial = add_spatial,
                      mesh = mesh,
                      num_archetypes = num_archetypes,
                      uncertainty_quantification = uncertainty_quantification,
                      spp_effects = do_em$new_spp_effects,
                      betas = do_em$new_betas,
                      mixture_proportion = do_em$new_mixprop,
                      spp_nuisance = as.data.frame(do_em$new_transformed_nuisance),
                      posterior_probability = do_em$post_prob)

    #' ## Predict species-specific field in an ad-hoc but scalable manner based on the most likely archetype that the species belong.
    get_spatial_fields <- foreach(l = 1:num_spp, .combine = "cbind") %dopar% .predict_spatial_fields(l = l,
                                                                                                     add_spatial = add_spatial,
                                                                                                     mesh = final_mesh,
                                                                                                     qa_object = get_qa, 
                                                                                                     em_object = do_em, 
                                                                                                     family = family,
                                                                                                     which_spp_effects = which_spp_effects)
    if(add_spatial) {
        rownames(get_spatial_fields) <- rownames(y)
        colnames(get_spatial_fields) <- colnames(y)
        }
        
    get_eta <- tcrossprod(X[, -which_spp_effects, drop = FALSE], out_assam$betas)
    out_assam$linear_predictor <- array(NA, dim = c(num_unit, num_spp, num_archetypes),
                                        dimnames = list(units = rownames(y), spp = colnames(y), archetype = names(out_assam$mixture_proportion)))
    for(k0 in 1:num_archetypes) {
        out_assam$linear_predictor[,,k0] <- tcrossprod(X[, which_spp_effects, drop = FALSE], out_assam$spp_effects) + matrix(get_eta[,k0], nrow = num_unit, ncol = num_spp, byrow = FALSE) + offset
        if(add_spatial)
            out_assam$linear_predictor[,,k0] <- out_assam$linear_predictor[,,k0] + get_spatial_fields
        }
    rm(get_eta)
     
    out_assam$spatial_fields <- get_spatial_fields
    out_assam$logL <- do_em$new_logL
    out_assam$df <- num_spp*(num_nuisance_perspp + length(which_spp_effects)) + prod(dim(out_assam$betas)) + (num_archetypes - 1)
    out_assam$control <- control
    out_assam$bootstrap_control <- bootstrap_control
    out_assam$beta_selection_control <- beta_selection_control
    if(uncertainty_quantification & bootstrap_control$method == "full_bootstrap") {
        save(get_qa, file = file.path(tempdir(), "allsdmTMBfits.RData")) # Save sdmTMB_fits into a temporary directory and remove it before running bootstrap...saves a lot of memory which in turn speeds up bootstrapping by a lot!
        rm(get_qa)
        }
    if(uncertainty_quantification & bootstrap_control$method == "fast_bootstrap") {
        out_assam$sdmTMB_fits <- get_qa$sdmTMB_fits
        get_qa$sdmTMB_fits <- NULL
        }
    if(!uncertainty_quantification) {
        out_assam$sdmTMB_fits <- get_qa$sdmTMB_fits
        get_qa$sdmTMB_fits <- NULL
        }
    
    
    ##----------------
    #' # Standard Error using full or fast but crude parametric bootstrap approach
    #' Computation time for full parametric bootstrap is crap at the moment unless you have a HPC!!! Starts from estimated parameters to give a little speed on for the quadratic approximations, so not used...
    #' In the case model selection is performed i.e., beta_selection = TRUE, **bootstrap is performed conditional on selected model** 
    ##----------------
    if(uncertainty_quantification & bootstrap_control$method == "fast_bootstrap") {
        message("Performing a fast but crude parametric bootstrap approach to obtain uncertainty quantification. Please take the results with a grain of salt!")
        
        #' ## Bootstrap datasets
        bootresp <- .fastsimulate_assam(qa_object = get_qa,
                                        em_object = do_em,
                                        num_X = num_X,
                                        which_spp_effects = which_spp_effects,
                                        nsim = bootstrap_control$num_boot,
                                        do_parallel = do_parallel,
                                        num_cores = num_cores,
                                        seed = bootstrap_control$seed)
        
        
        #' ## Fit asSAM to each bootstrapped dataset
        bootcov_fast_fn <- function(b0) { 
            get_boot_qa <- get_qa
            get_boot_qa$parameters <- bootresp[[b0]]$bootstrap_parameters
            get_boot_qa$long_parameters <- apply(get_boot_qa$parameters, 1, function(x) kronecker(rep(1, num_archetypes), x)) # Repeats species-specific estimates num_archetypes times; object has num_spp columns
            
            ##----------------
            #' # Run EM algorithm for bootstrap datasets
            ##----------------
            if(!beta_selection)
                do_boot_em <- try(em_fn(qa_object = get_boot_qa,
                                        dobar_penalty = FALSE), silent = TRUE)
            if(beta_selection)
                do_boot_em <- try(em_fn(qa_object = get_boot_qa, 
                                        dobar_penalty = FALSE,
                                        betamatrix_selection = out_assam$betas), silent = TRUE)
            
            setTxtProgressBar(pb, b0)
            return(do_boot_em)
            }          
        
        pb <- txtProgressBar(min = 0, max = bootstrap_control$num_boot, style = 3)
        bootrun <- lapply(1:bootstrap_control$num_boot, bootcov_fast_fn)
        close(pb) 
        rm(pb)
        }
    
    if(uncertainty_quantification & bootstrap_control$method == "full_bootstrap") {
        message("Performing full parametric bootstrap to obtain uncertainty quantification...this will take a while so go a brew a cup of tea (or two)!")
        
        #' ## Bootstrap datasets
        class(out_assam) <- "assam"
        bootresp <- simulate.assam(out_assam,
                                   data = data,
                                   nsim = bootstrap_control$num_boot,
                                   do_parallel = do_parallel,
                                   num_cores = num_cores,
                                   seed = bootstrap_control$seed)
        

        #' ## Fit asSAM to each bootstrapped dataset
        bootcov_fn <- function(b0, control) { 
            ##----------------
            #' ## Construct quadratic approximations for each species in bootstrap dataset, and set up relevant quantities
            ##----------------
            get_boot_qa <- try(.quadapprox2_fn(family = family, 
                                               formula = formula, 
                                               resp = bootresp[,b0]$y, 
                                               data = data, 
                                               add_spatial = add_spatial,
                                               mesh = final_mesh,
                                               offset = offset,
                                               trial_size = trial_size,
                                               do_parallel = do_parallel,
                                               control = control,
                                               return_fits = FALSE),
                               silent = TRUE)
            if(inherits(get_boot_qa, "try-error"))
                return(get_boot_qa)
            
            get_boot_qa$long_parameters <- apply(get_boot_qa$parameters, 1, function(x) kronecker(rep(1,num_archetypes), x))
               
            ##----------------
            #' # Run EM algorithm for bootstrap datasets
            ##----------------
            if(!beta_selection)
                do_em <- try(em_fn(qa_object = get_boot_qa,
                                   dobar_penalty = FALSE), silent = TRUE)
            if(beta_selection)
                do_em <- try(em_fn(qa_object = get_boot_qa,
                                   dobar_penalty = FALSE,
                                   betamatrix_selection = out_assam$betas), silent = TRUE)
            
            setTxtProgressBar(pb, b0)
            return(do_em)
            }          
        
        pb <- txtProgressBar(min = 0, max = bootstrap_control$num_boot, style = 3)
        bootrun <- lapply(1:bootstrap_control$num_boot, bootcov_fn, control = control)
        close(pb) 
        rm(pb)

        load(file = file.path(tempdir(), "allsdmTMBfits.RData"))
        out_assam$sdmTMB_fits <- get_qa$sdmTMB_fits
        get_qa$sdmTMB_fits <- NULL
        file.remove(file.path(tempdir(), "allsdmTMBfits.RData"))
        }
        
    if(uncertainty_quantification) {
        find_errors <- which(sapply(bootrun, function(x) inherits(x, "try-error")))
        if(length(find_errors) > 0) {
            warning("Bootstrap datasets #", find_errors, " encountered problems during fitting, and are subsequently ignored...\nIf the number of datasets with fitting problems is large, it may point to deeper issues with the asSAM itself.")
            bootrun <- bootrun[-find_errors]
            }
        rm(find_errors)
        gc()
        
        
        #' ## Account for potential label-switching across the bootstrapped datasets, and also transform nuisance parameters as appropriate
        if(num_archetypes > 1) {
            boot_params <- lapply(bootrun, function(x) cbind(x$new_betas, x$new_mixprop))
            boot_params <- abind::abind(boot_params, along = 0)
            switch_labels <- label.switching::pra(mcmc.pars = boot_params, 
                                                  pivot = cbind(out_assam$betas, out_assam$mixture_proportion))$permutations #' This should technically be avoided when model selection has been performed and a particular sparsity pattern is sought. But in practice hopefully this makes little difference!
            }
        if(num_archetypes == 1) {
            boot_params <- NULL
            switch_labels <- NULL
            }
        
        bootrun <- lapply(1:length(bootrun), function(k0) {
            if(num_archetypes > 1) {
                bootrun[[k0]]$new_mixprop <- bootrun[[k0]]$new_mixprop[switch_labels[k0,]]
                bootrun[[k0]]$new_betas <- bootrun[[k0]]$new_betas[switch_labels[k0,],]
                #bootrun[[k0]]$post_prob <- bootrun[[k0]]$post_prob[,switch_labels[k0,]] #' Not actually used later on so omit! 
                }
            
            #' A vector of all parameters ordered in the same way as the columns of the mapping matrix, plus the mixture proportions
            if(num_nuisance_perspp == 0)
                bootrun[[k0]]$boot_params <- c(as.vector(t(bootrun[[k0]]$new_spp_effects)), 
                                               as.vector(t(bootrun[[k0]]$new_betas)), 
                                               bootrun[[k0]]$new_mixprop)
            if(num_nuisance_perspp > 0)
                bootrun[[k0]]$boot_params <- c(as.vector(t(bootrun[[k0]]$new_spp_effects)), 
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
        form_cis$spp_effects <- apply(abind::abind(lapply(bootrun, function(x) x$new_spp_effects), along = 0), c(2,3), quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)
        form_cis$spp_effects <- list(lower = form_cis$spp_effects[1,,,drop=FALSE], upper = form_cis$spp_effects[2,,,drop=FALSE])
        if(length(which_spp_effects) == 1) {
            form_cis$spp_effects$lower <- matrix(form_cis$spp_effects$lower, nrow = 1)
            form_cis$spp_effects$upper <- matrix(form_cis$spp_effects$upper, nrow = 1)
            }
        
        form_cis$betas <- apply(abind::abind(lapply(bootrun, function(x) x$new_betas), along = 0), c(2,3), quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)
        form_cis$betas <- list(lower = form_cis$betas[1,,], upper = form_cis$betas[2,,])
        if(num_archetypes == 1) {
            form_cis$betas$lower <- matrix(form_cis$betas$lower, nrow = 1)
            form_cis$betas$upper <- matrix(form_cis$betas$upper, nrow = 1)
            }
        
        if(num_nuisance_perspp > 0) {
            form_cis$spp_nuisance <- apply(abind::abind(lapply(bootrun, function(x) x$new_transformed_nuisance), along = 0), c(2,3), quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)
            form_cis$spp_nuisance <- list(lower = form_cis$spp_nuisance[1,,,drop = FALSE], upper = form_cis$spp_nuisance[2,,,drop = FALSE])
            
            }
        if(num_archetypes > 1)
            form_cis$mixture_proportion <- data.frame( t(apply(sapply(bootrun, function(x) x$new_mixprop), 1, quantile, prob = c(modified_alpha/2, 1 - modified_alpha/2), na.rm = TRUE)))
        
        rownames(form_cis$spp_effects$lower) <- rownames(form_cis$spp_effects$upper) <- colnames(out_assam$spp_effects)
        if(num_archetypes > 1) {
            colnames(form_cis$mixture_proportion) <- c("lower", "upper")
            }
        rownames(form_cis$betas$lower) <- rownames(form_cis$betas$lower) <- rownames(out_assam$betas)
        if(num_nuisance_perspp > 0) 
            dimnames(form_cis$spp_nuisance$lower)[[3]] <- dimnames(form_cis$spp_nuisance$upper)[[3]] <- colnames(out_assam$spp_nuisance)
        
        out_assam$confidence_intervals <- form_cis
        out_assam$bootstrap_parameters <- t(sapply(bootrun, function(x) x$boot_params))
        
        
        #' ## Calculate bootstrapped posterior probabilities of *original species data* belong to each archetype
        .calc_posterior_prob <- function(cw_bootstrap_spp_effects, cw_bootstrap_betas, cw_bootstrap_nuisance, cw_bootstrap_mixprop, qa_object) {
            post_prob <- matrix(NA, nrow = num_spp, ncol = num_archetypes)
            rownames(post_prob) <- colnames(y)
            colnames(post_prob) <- paste0("archetype", 1:num_archetypes)
            logL_spp <- NULL         
            
            for(j in 1:num_spp) { 
                cw_Quad <- sapply(1:num_archetypes, function(k) {
                    cw_params <- c(cw_bootstrap_spp_effects[j,], cw_bootstrap_betas[k,], cw_bootstrap_nuisance[j,])
                    cw_params[which_spp_effects] <- cw_bootstrap_spp_effects[j,]
                    cw_params[(1:num_X)[-which_spp_effects]] <- cw_bootstrap_betas[k,]
                    
                    cw_v <- matrix(cw_params - qa_object$parameters[j,], ncol = 1)
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
            .calc_posterior_prob(cw_bootstrap_spp_effects = bootrun[[k0]]$new_spp_effects, 
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
.predict_spatial_fields <- function(l, 
                                    add_spatial, 
                                    mesh, 
                                    qa_object, 
                                    em_object, 
                                    family,
                                    which_spp_effects) {
    if(!add_spatial)
        return(NULL)
    
    if(add_spatial) {
        #' Recompile the TMB object from the original sdmTMB fits, evaluated at the assam parameter values
        #' Note an nlminb is required to update the environment parameter values 
        use_pars <- qa_object$sdmTMB_fits[[l]]$tmb_obj$env$parList()
        cw_b_j <- c(em_object$new_spp_effects[l,], em_object$new_betas[which.max(em_object$post_prob[l,]), ])
        cw_b_j[which_spp_effects] <- em_object$new_spp_effects[l,]
        cw_b_j[-which_spp_effects] <- em_object$new_betas[which.max(em_object$post_prob[l,]), ] 
        use_pars[["b_j"]] <- cw_b_j
        rm(cw_b_j)
        
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
                                      parameters = use_pars,
                                      map = use_map,
                                      random = qa_object$sdmTMB_fits[[l]]$tmb_random,
                                      DLL = "sdmTMB",
                                      silent = TRUE)

        #' Method 1 -- Based on predict.sdmTMB but slower due to the use of sdreport to extract random effects
        # new_tmb_obj$fn(new_tmb_obj$par) # need to initialize the new TMB object once
        # new_tmb_sdreport <- TMB::sdreport(new_tmb_obj, par.fixed = new_tmb_obj$par) # Update random effects
        # r <- new_tmb_obj$report(new_tmb_obj$env$last.par) # last.par taken since it is the newest set of parameters
        
        #' Method 2 -- Own approach which is faster but requires running a new, single optimization
        new_fit0 <- nlminb(start = new_tmb_obj$par,
                           objective = new_tmb_obj$fn,
                           gradient = new_tmb_obj$gr,
                           control = qa_object$sdmTMB_fits[[l]]$nlminb_control)

        #' Use TMB's report function to return the estimated field of assam
        r <- new_tmb_obj$report(new_tmb_obj$env$last.par.best)
        return(as.vector(r$omega_s_A))
        }
    }



# Hidden function taken directly from sdmTMB utils.R
#' @noMd
#' @noRd
.get_pars <- function(object) {
    # based on glmmTMB:
    ee <- object$tmb_obj$env
    x <- ee$last.par.best
    if (length(ee$random) > 0) x <- x[-ee$random]
    p <- ee$parList(x = x)
    p
    }




#' score_fn <- function(j) {
#'     #' Psi sequence: species-specific effects; archetypal regression coefficients by archetype; species-specific nuisance parameters, along with parameters for species-specific spatial fields
#'     some_scores <- function(k0) {
#'         cw_params <- c(do_em$new_spp_effects[j,], do_em$new_betas[k0,], do_em$new_nuisance[j,])
#'         cw_params[which_spp_effects] <- do_em$new_spp_effects[j,]
#'         cw_params[(1:num_X)[-which_spp_effects]] <- do_em$new_nuisance[k0,]
#'         
#'         Infoproddiff <- -get_qa$hessian[[j]] %*% (cw_params - get_qa$parameters[j,])
#'         cw_spp_effects <- Infoproddiff[which_spp_effects]
#'         cw_beta <- Infoproddiff[(1:num_X)[-which_spp_effects]].
#'         cw_nuisance <- Infoproddiff[num_X + 1:num_nuisance_perspp]
#'         
#'         return(list(spp_effects = matrix(cw_spp_effects, ncol = 1), 
#'                     beta = matrix(cw_beta, ncol = 1), 
#'                     spp_nuisance = matrix(cw_nuisance, ncol = 1)))
#'     }
#'     all_scores <- lapply(1:num_archetypes, some_scores)
#'     
#'     score_spp_effects <- rowSums(matrix(out_assam$posterior_probability[j,], nrow = length(which_spp_effects), ncol = num_archetypes, byrow = TRUE) * matrix(sapply(all_scores, function(x) x$spp_effects), ncol = num_archetypes))
#'     score_beta <- matrix(out_assam$posterior_probability[j,], nrow = length((1:num_X)[-which_spp_effects]), ncol = num_archetypes, byrow = TRUE) * matrix(sapply(all_scores, function(x) x$beta), ncol = num_archetypes)
#'     score_beta <- as.vector(score_beta)
#'     score_nuisance <- rowSums(matrix(out_assam$posterior_probability[j,], nrow = num_nuisance_perspp, ncol = num_archetypes, byrow = TRUE) * matrix(sapply(all_scores, function(x) x$spp_nuisance), ncol = num_archetypes))
#'     
#'     score_psi <- numeric(ncol(mapping_mat))
#'     names(score_psi) <- colnames(mapping_mat)
#'     score_psi[grep(paste0("spp_effects",j,"_"), names(score_psi))] <- score_spp_effects
#'     score_psi[grep("archetype", names(score_psi))] <- score_beta
#'     score_psi[grep(paste0("spp_nuisance",j,"[:alpha:]"), names(score_psi))] <- score_nuisance
#'     rm(all_scores, score_spp_effects, score_beta, score_nuisance)            
#'     
#'     #' Note is omits the mixing proportion of the last archetype due to sum to one constraint
#'     score_mixprop <- out_assam$posterior_probability[j,-num_archetypes]/out_assam$mixture_proportion[-num_archetypes] - out_assam$posterior_probability[j,num_archetypes]/out_assam$mixture_proportion[num_archetypes] 
#'     
#'     out <- c(score_mixprop, score_psi)
#'     return(out)
#' }
#' 
#' all_scores <- foreach(j = 1:num_spp, .combine = "cbind") %dopar% score_fn(j)        
#' 
#' empirical_info <- matrix(rowSums(apply(all_scores, 2, tcrossprod)), nrow = nrow(all_scores)) 
#' rownames(empirical_info) <- colnames(empirical_info) <- c(paste0("mixture_proportion", 1:(num_archetypes-1)), colnames(mapping_mat))