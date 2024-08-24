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
#' The current summary function is pretty basic (apologies!), and in the future we may add some more useful information to the summary output. 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
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

