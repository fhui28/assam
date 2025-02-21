#' @title Print a (approximate and scalable) SAM object
#' 
#' @description 
#' `r lifecycle::badge("stable")
#'
#' The default print method for a \code{assam} object.
#' 
#' @param x An object of class \code{assam}.
#' @param ... Not used.
#'
#' @details
#' Print out details such as the function call, assumed family/response distribution, number of sites and species, and archetypal response-environment relationship fitted.
#'
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' @examples
#' \dontrun{
#' #' Please see the help file for assam for example.
#' }
#' 
#' @method print assam 
#' @export
#' @md



print.assam <- function(x, ...) {
    message("Call:")
    print(x$call)
    message()
    
    message("--- Approximate and scalable species archetype model (asSAM) ---") 
    message("Family: ", x$family$family[1], "\nNo. of sites: ", dim(x$linear_predictor)[1], "\nNo. of species: ", dim(x$linear_predictor)[2]) 
    message("No. of archetypes: ", nrow(x$betas)) 
    message("Archetypal response-environment relationship: ", x$formula) 
    message("Were constraints imposed on the archetypal regression coefficients? (TRUE/FALSE): ", (!is.null(x$control$beta_lower) | !is.null(x$control$beta_upper))) 
    message("Was variable selection performed on the regression coefficients? (TRUE/FALSE): ", x$beta_selection) 
    message("Were confident intervals constructed? (TRUE/FALSE): ", x$uncertainty_quantification) 
    message("Were species-specific spatial fields included? (TRUE/FALSE): ", x$add_spatial)
    message("---") 
    }
