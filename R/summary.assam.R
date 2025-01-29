#' @title Summary for an approximate and scalable SAM object
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{assam} object and produces some basic summary information.
#'
#' @param object An object of class \code{assam}.
#' @param digits The number of significant figures to use when printing.
#' @param ... Not used.
#' 
#' @details 
#' Currently, the function returns estimated mixture proportions, species-specific intercepts, and archetypal regression coefficients. If the \code{object$uncertainty_quantification == TRUE}, then summary tables are also produced for archetypal regression coefficients containing parametric bootstrap confidence intervals and statistical significance based off this. Note the length of the confidence intervals was determined as part of the fitting the asSAM; see \code{object$bootstrap_control$ci_alpha}.
#' 
#' @return An object of class "summary.assam" which includes the following components, not necessarily in the order below (and as appropriate):
#' \describe{
#' \item{call: }{The matched function call of \code{object}.}
#' \item{mixture_proportion:}{Estimated vector of mixture proportions corresponding to the probability of belonging to each archetype.}
#' \item{spp_intercepts:}{Estimated species-specific intercepts.}
#' \item{betas:}{Estimated matrix of archetypal regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of archetypes.}
#' \item{betas_results: }{If the \code{object$uncertainty_quantification = TRUE}, a summary table corresponding to archetypal regression coefficients.}
#' }
#' @details # Warning
#' - The current summary function is pretty basic (apologies!), and in the future we may add some more useful information to the summary output.
#' 
#' - Inferential results on the archetypal regression coefficients should be interpreted with a large grain of salt if lower/upper constraints were placed on these during the fitting pro 
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
#' set.seed(022025)
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
#'  simdat <- create_samlife(family = nbinom2(), 
#'  formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula, 
#'  data = covariate_dat, 
#'  betas = true_betas, 
#'  spp_intercept = true_intercepts, 
#'  spp_dispparam = true_dispparam, 
#'  spp_powerparam = true_powerparam, 
#'  mixture_proportion = true_mixprop,
#' seed = 022025)
#'  
#'  
#'  ##----------------------
#'  # Fit asSAM and assess results 
#'  #' **Most users should start here**
#'  ##----------------------
#'  samfit <- assam(y = simdat$y,
#'  formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
#'  data = covariate_dat,
#'  family = nbinom2(),
#'  num_archetypes = num_archetype,
#'  num_cores = 8)
#'  
#'  plot(true_intercepts, samfit$spp_intercepts); abline(0,1)
#'  plot(true_dispparam, samfit$spp_nuisance$dispersion, log = "xy"); abline(0,1)
#'  #' Note estimates for the archetypal responses and mixture proportions from (as)SAMs should be 
#'  #' close to the corresponding true values, *up to a reordering* of the mixture component
#'  #' s/archetypes (since the order is essentiually arbtirary)
#'  rbind(true_betas, samfit$betas) %>% 
#'  t %>% 
#'  as.data.frame %>%
#'  ggpairs
#'  table(simdat$archetype_label, apply(samfit$posterior_probability, 1, which.max))
#'  
#'  
#'  samfit
#'  summary(samfit)
#' }
#' 
#' 
#' @method summary assam
#' @export
#' @md

summary.assam <- function(object, digits = 4, ...) {
    summary_output <- list(call = object$call, 
                           mixture_proportion = round(object$mixture_proportion, digits),
                           spp_intercepts = round(object$spp_intercepts, digits),
                           betas = round(object$betas, digits))
    
    if(object$uncertainty_quantification) {          
        betas_resultstab <- data.frame(as.data.frame.table(object$betas),
                                       lower = as.data.frame.table(round(object$confidence_intervals$betas$lower, digits))[,3], 
                                       upper = as.data.frame.table(round(object$confidence_intervals$betas$upper, digits))[,3])
        betas_resultstab$excludeszero <- as.numeric(betas_resultstab$lower > 0 | betas_resultstab$upper < 0)
        colnames(betas_resultstab)[1:3] <- c("Species", "Covariate", "Estimate")
        
        summary_output$betas_results <- betas_resultstab
        }
    
    
    class(summary_output) <- "summary.assam"
    summary_output 
    }	

