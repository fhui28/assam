#' @title Penalized approximate and scalable species archetype models (asSAMs)
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Fits approximate and scalable species archetype modeling (asSAMs) for model-based clustering of species based on their environmental response, into a small number of so-called archetypal responses. This is a modification of the main [assam()] function where a penalty is augmented to the approximate log-likelihood function to encourage sparsity in the archetypal regression coefficients.
#' 
#' 
#' @param y A multivariate abundance response matrix.
#' @param formula An object of class "formula", which represents a symbolic description of the full model matrix to be created (based on using this argument along with the \code{data} argument). *Note there should be nothing on the left hand side of the "~".* Currently, smooth terms are not permitted.
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param which_spp_effects A vector identifying which columns of the model matrix induced by \code{formula} and \code{data} should be treated as species-specific effects. Default to 1, meaning only the first column i.e., the intercept, is species-specific.
#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param offset A matrix of offset terms, of the same dimension as \code{y}.
#' @param trial_size Trial sizes to use for binomial distribution. This should equal to a scalar.
#' @param num_archetypes Number of archetypes (clusters) to assume in the asSAM.
#' @param mesh Output from [sdmTMB::make_mesh()], used for adding species-specific spatial fields to the linear predictor.
#' @param do_parallel Should parallel computing be used to fit the asSAM. Defaults to \code{FALSE}.
#' @param num_cores If \code{do_parallel = TRUE}, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to \code{parallel::detectCores() - 2}.
#' @param supply_quadapprox An object of class \code{assam_quadapprox}, which is (mostly likely) obtained as a consequence of running an initial fit using [assam()] with \code{do_assam_fit = FALSE}. 
#' @param nlambda The number of tuning parameters values to consider when forming the regularization path.
#' @param lambda_min_ratio The smallest value for the tuning parameter lambda, as a fraction of the maximum internally-derived lambda value.
#' @param lambda A user-supplied tuning parameter sequence. Note this is usually not supplied, as it is standard to have the function itself compute its own lambda sequence based on \code{nlambda} and \code{lambda_min_ratio}.
#' @param min_df The minimum number of non-zero archetypal regression coefficients allowed in the asSAM. This is useful to supply if, when the function tries to find a appropriate lambda sequence, the user wants the maximum value of lambda to not necessarily shrink to all the penalized estimates to zero. Defaults to zero.
#' @param control A list containing the following elements:
#' \describe{
#' \item{max_iter:}{the maximum number of iterations in the EM algorithm. Usually convergence is quite quick e.g., less than 20 iterations.}
#' \item{tol:}{the convergence criterion; the difference in the log-likelihood value of the asSAM from successive iterations must be smaller than this value.}
#' \item{temper_prob:}{in the iteration of the EM algorithm, posterior probabilities from the E-step are "tempered" or push away from the 0/1 boundary. This is often useful to get the EM algorithm moving initially.}
#' \item{trace:}{controls if messages printed as part of the estimation process to reflect progress.}
#' }
#' @param beta_selection_control A list containing the following elements to control the broken adaptive ridge (BAR) penalty for variable selection on the archetypal regression coefficients:
#' \describe{
#' \item{max_iter:}{the maximum number of iterations in the BAR optimization part of the EM algorithm.}
#' \item{eps:}{the convergence criterion; the norm of the difference between all estimated parameters from successive iterations must be smaller than this value.}
#' \item{round_eps:}{a tolerance to round values to zero. The technically not needed as the BAR penalty will produce exactly zero estimates up to machine error, but is included anyway, but is included anyway.}
#' }
#' 
#' 
#' @details 
#' For the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = u_i^\top\alpha_j + x_i^\top\beta_k,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{u_i^\top\alpha_j} corresponds to a component that is to kept species-specific e.g., species-specific intercept, \eqn{x_i^\top\beta_k}  denotes the component corresponding to effect of archetypal response \eqn{k}. Additionally, species-specific spatial fields can be included in the linear predictor e.g., to account for residual spatial correlation above and beyond that explained by the archetypal responses. Conditional on the mean model above, the \eqn{y_{ij}} are assumed to be come from some response distribution using the additional dispersion and power parameters as appropriate. We refer the reader to the main [assam()] function for more details and references for SAMs.
#' 
#' This function is specifically designed for performing variable selection, and building an associated regularization path, on the elements of \eqn{\beta_k} via the broken adaptive ridge [BAR](https://www.sciencedirect.com/science/article/pii/S0047259X17305067) penalty. The BAR penalty can be interpreted as a kind of approximation to the \eqn{L_0} penalty, and encourages sparsity in the archetypal regression coefficients e.g., to uncover what covariates are informative for each of the archetypal responses. After building a regularization path using this function, one can examine the path and select the tuning parameter minimizing one of the provided information criterion, say, to pass to [assam()] function for fitting the final, sparse asSAM.
#' 
#' @section Distributions and parallelization:
#' Please see the [assam()] function for more details on the distributions currently supported by the package, and how to properly parallelize the fitting process.
#' 
#' @return An object of class \code{passam} with the following elements (as appropriate, and not necessarily in the order below):
#' \item{call:}{The function call.}
#' \item{formula:}{Same as input argument.}
#' \item{family:}{Same as input argument.}
#' \item{num_nuisance_perspp:}{The number of "nuisance" parameters per species. For example, if \code{family = binomial()} and species-specific spatial fields are not included, then there are no nuisance parameters. If \code{family = nbinom2()} and species-specific spatial fields are included, say, then there are three nuisance parameters per species.}
#' \item{trial_size:}{Same as input argument.}
#' \item{offset:}{Same as input argument.}
#' \item{add_spatial:}{Were species-specific spatial fields included?}
#' \item{mesh:}{Same as input argument.}
#' \item{num_archetypes:}{Same as input argument.}
#' \item{lambda:}{The actual sequence of tuning parameter, lambda, values used.}
#' \item{betas_df:}{The number of non-zero archetypal regression coefficients at each value of lambda.}
#' \item{betas_path:}{The estimated archetypal regression coefficients at each value of lambda. This takes the form of a three-dimensional array, where the third dimension corresponds to the lambda. So example \code{betas_path[,,1]} corresponds to the estimated archetypal regression coefficients at the first value of lambda.}
#' \item{logL:}{Estimated log-likelihood value of the asSAM i.e. the value of the approximated log-likelihood function, at each value of lambda.}
#' \item{AIC:}{The Akaike Information Criterion at each value of lambda. The AIC is defined as \eqn{-2\log(L) + 2df}, where \eqn{L} is the log-likelihood value and \eqn{df} is the number of non-zero archetypal regression coefficients.}
#' \item{BIC:}{The Bayesian Information Criterion at each value of lambda. The BIC is defined as \eqn{-2\log(L) + \log(M)df}, where \eqn{M} is the number of species.}
#' \item{BIC2:}{The Bayesian Information Criterion, version 2, at each value of lambda. The BIC is defined as \eqn{-2\log(L) + \log(MN)df}, where \eqn{M} is the number of species and \eqn{N} is the number of observational units. Compared to BIC, BIC2 uses a more severe model complexity penalty and thus is more likely to selection a sparse asSAM.}
#' \item{regularization_frame:}{A data frame combining some of the above information for easy "digestion".}
#' \item{control:}{Same as input argument.}
#' \item{beta_selection_control:}{Same as input argument.}
#' \item{quad_approx:}{An object of class \code{assam_quadapprox}, which contains the ingredients needed to fit asSAMs. This is useful for supplying to downstream applications of the [assam()] function e.g., see the example below.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @examples
#' \dontrun{
#' ##----------------------
#' # Generate some multivariate abundance (count) data from a sparse SAM
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
#' true_betas[which(abs(true_betas) < 0.4)] <- 0 # Making archetypal coefficients sparse
#' true_spp_effects <- matrix(runif(num_spp, -2, 0), ncol = 1)
#' true_dispparam <- 1/runif(num_spp, 1, 5) 
#' true_powerparam <- runif(num_spp, 1.4, 1.8)
#' true_mixprop <- c(0.2, 0.2, 0.3, 0.15, 0.15)
#' true_betas
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
#' ## First construct regularization path for asSAMs  
#' #' **Most users should start here**
#' # Note this can take a bit of time...apologies!
#' samfit_select <- passam(y = simdat$y,
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' family = nbinom2(),
#' min_df = 5, #' Minimum tuning parameter is such that there are at least five non-zero coefficients
#' num_archetypes = num_archetype,
#' control = list(trace = 1),
#' num_cores = detectCores() - 2)
#' 
#' samfit_select$regularization_frame
#' use_BIC2 <- which.min(samfit_select$BIC2)
#' samfit_select$betas_path[,,use_BIC2]
#' 
#' 
#' ## Now fit the final asSAMs given a chosen value of the tuning parameter
#' # Note the supply of the quadratic approximations to speed up [assam()]!
#' samfit <- assam(y = simdat$y,
#' formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#' data = covariate_dat,
#' family = nbinom2(),
#' beta_selection = TRUE,
#' uncertainty_quantification = FALSE,
#' num_archetypes = num_archetype,
#' supply_quadapprox = samfit_select$quad_approx,
#' beta_selection_control = list(lambda = samfit_select$lambda[use_BIC2]),
#' num_cores = detectCores() - 2)
#' 
#' 
#' samfit
#' samfit$betas
#' true_betas
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
#' }


#' @export passam
#'
#' @importFrom abind abind
#' @importFrom cluster pam
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom sdmTMB sdmTMB sdmTMBcontrol make_mesh nbinom2 tweedie Beta
#' @importFrom label.switching pra
#' @importFrom methods as
#' @importFrom parallel detectCores
#' @importFrom stats as.formula nlminb model.matrix qt plogis qlogis 
#' @importFrom TMB MakeADFun
#' @import Matrix
#' @md

passam <- function(y, 
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
                  supply_quadapprox = NULL,
                  nlambda = 100, 
                  lambda_min_ratio = 1e-6,
                  lambda = NULL, 
                  min_df = 0,
                  control = list(max_iter = 500, tol = 1e-5, 
                                 temper_prob = 0.8, trace = FALSE),
                  beta_selection_control = list(max_iter = 100, eps = 1e-5, round_eps = 1e-6)) {
    
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
    beta_selection_control <- .fill_beta_selection_control(control = beta_selection_control, lambda_fill = FALSE)

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
    #' # Set up penalized EM algorithm function
    ##----------------
    pem_fn <- function(qa_object, lambda, beta_selection_control) {
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
                do_kmeans <- cluster::pam(qa_object$parameters[, (1:num_X)[-which_spp_effects], drop = FALSE], k = num_archetypes)  
                cw_betas <- do_kmeans$medoids
                cw_spp_effects <- qa_object$parameters[, which_spp_effects, drop = FALSE]
                cw_nuisance <- NULL
                if(num_nuisance_perspp > 0)
                    cw_nuisance <- qa_object$parameters[, num_X + 1:num_nuisance_perspp, drop = FALSE]
                cw_mixprop <- as.vector(table(do_kmeans$clustering)) / num_spp
                rm(do_kmeans)
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
                }
            
            ##-------------------
            #' ## M-step
            #' Make the unpenalized estimate at the current M-step first (more for formatting reasons, and it has relatively little overhead, before penalizing)
            ##-------------------
            new_mixprop <- colMeans(post_prob)
            bigW <- Diagonal(x = rep(as.vector(t(sqrt(post_prob))), each = nrow(qa_object$hessian[[j]]))) %*% basic_bigW %*% Diagonal(x = rep(as.vector(t(sqrt(post_prob))), each = nrow(qa_object$hessian[[j]])))
            MtWM <- forceSymmetric(crossprod(mapping_mat, bigW) %*% mapping_mat)
            new_params <- solve(MtWM, crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters))
            new_params <- as.vector(new_params)
            
            
            beta_s_err <- Inf
            beta_s_counter <- 0
            cw_params <- new_params
            Dbar <- Diagonal(n = length(new_params))
            diag(Dbar)[-grep("archetype", colnames(mapping_mat))] <- 0
            lambda_used <- num_spp * lambda 
            
            while(beta_s_counter < beta_selection_control$max_iter & beta_s_err > beta_selection_control$eps) {
                GammaMatrix_params <- Diagonal(n = length(new_params))
                diag(GammaMatrix_params)[grep("archetype", colnames(mapping_mat))] <- cw_params[grep("archetype", colnames(mapping_mat))]
                
                new_params <- GammaMatrix_params %*% solve(GammaMatrix_params %*% MtWM %*% GammaMatrix_params + lambda_used * Dbar) %*% GammaMatrix_params %*% crossprod(mapping_mat, bigW) %*% as.vector(qa_object$long_parameters)
                new_params <- as.vector(new_params)
                names(new_params) <- colnames(mapping_mat)
                
                beta_s_err <- sum((new_params - cw_params)^2)
                cw_params <- new_params
                beta_s_counter <- beta_s_counter + 1
                }
            #' Note the convergence properties of the algorithm does ensure exact zeros can be achieved (up to machine error); see [https://doi.org/10.1155/2016/3456153]. However once can also round to zero if interested.
            new_params[abs(new_params) < beta_selection_control$round_eps] <- 0
            rm(beta_s_counter, beta_s_err, GammaMatrix_params, Dbar)
            
            
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
            if(control$trace) {
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
    ##----------------
    ## Attempt to determine a minimum and maximum lambda, if required, for BAR penalty
    ##----------------
    if(!is.null(lambda))
        lambdaseq <- lambda
    if(is.null(lambda)) {
        message("Finding an appropriate sequence of lambda values...")
        lambda_max <- 1
        fast_beta_selection_control <- beta_selection_control
        fast_beta_selection_control$max_iter <- 20 
        fit_togetmaxlambda <- pem_fn(qa_object = get_qa, 
                                     lambda = lambda_max,
                                     beta_selection_control = fast_beta_selection_control)
        
        while(sum(fit_togetmaxlambda$new_betas != 0) > min_df) {
            lambda_max <- lambda_max*2
            fit_togetmaxlambda <- pem_fn(qa_object = get_qa, 
                                         lambda = lambda_max,
                                         beta_selection_control = fast_beta_selection_control)
            }
        
        lambdaseq <- lseq(lambda_max, lambda_max*lambda_min_ratio, length = nlambda, decreasing = TRUE)
        rm(fit_togetmaxlambda, fast_beta_selection_control)
        }
    
    ##----------------
    #' # Run penalized EM at each point on regularization path
    ##----------------
    tic <- proc.time()
    message("Constructing regularization path...")

    cwfit_fn <- function(l) {
        if(control$trace)
            message("Commencing penalized EM algorithm at lambda = ", round(lambdaseq[l], 4))
        
        cwfit <- pem_fn(qa_object = get_qa, 
                        lambda = lambdaseq[l], 
                        beta_selection_control = beta_selection_control)
        try_counter <- 0
        while(any(cwfit$new_mixprop < 1e-3) & try_counter < 20) {
            if(control$trace)
                message("Mixture component is being emptied...altering initial temp probability and restarting EM-algorithm to try and fix this.")
            control$temper_prob <- control$temper_prob + 0.025
            
            cwfit <- pem_fn(qa_object = get_qa, 
                            lambda = lambdaseq[l], 
                            beta_selection_control = beta_selection_control)
            try_counter <- try_counter + 1
            }
        
        return(cwfit)
        }
    
    allfits <- foreach(l = 1:length(lambdaseq)) %dopar% cwfit_fn(l = l)
    toc <- proc.time()
    gc()
    
    coefficients_path <- abind::abind(lapply(allfits, function(x) x$new_betas), along = 3)
    df_path <- sapply(allfits, function(x) sum(x$new_betas != 0))
    passam_logL_path <- sapply(allfits, function(x) x$new_logL)
    
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
                      add_spatial = add_spatial,
                      mesh = mesh,
                      num_archetypes = num_archetypes)
    
    out_assam$lambda <- lambdaseq
    out_assam$betas_df <- df_path
    out_assam$betas_path <- coefficients_path
    out_assam$logL <- passam_logL_path
    out_assam$AIC <- -2*passam_logL_path + 2*df_path
    out_assam$BIC <- -2*passam_logL_path + log(num_spp)*df_path
    out_assam$BIC2 <- -2*passam_logL_path + log(num_spp*num_unit)*df_path
    out_assam$regularization_frame <- data.frame(lambda = lambdaseq, 
                                                 df = df_path, 
                                                 logL = passam_logL_path, 
                                                 AIC = out_assam$AIC, 
                                                 BIC = out_assam$BIC, 
                                                 BIC2 = out_assam$BIC2)
    out_assam$control <- control
    out_assam$beta_selection_control <- beta_selection_control
    out_assam$quad_approx <- get_qa

    ##----------------
    #' # Done! Final touches if any
    ##----------------
    class(out_assam) <- "passam"
    return(out_assam)
    }
     


