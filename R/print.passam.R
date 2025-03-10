#' @title Print a approximate and scalable SAM (pasSAM) object
#' 
#' @description 
#' `r lifecycle::badge("stable")
#'
#' The default print method for a \code{passam} object.
#' 
#' @param x An object of class \code{passam}.
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
#' #' Please see the help file for assam and passam for example.
#' }
#' 
#' @method print passam 
#' @export
#' @md



print.passam <- function(x, ...) {
    message("Call:")
    print(x$call)
    message()
    
    message("--- Penalized approximate and scalable species archetype model (asSAM) ---") 
    message("Family: ", x$family$family[1], "\nNo. of sites: ", dim(x$quad_approx$sdmTMB_fits[[1]]$data)[1], "\nNo. of species: ", dim(x$quad_approx$parameters)[2]) 
    message("(Maximumz) No. of archetypes: ", x$num_archetypes) 
    message("Archetypal response-environment relationship: ", x$formula) 
    message("Were species-specific spatial fields included? (TRUE/FALSE): ", x$add_spatial)
    message("")
    message("Please see regularization_frame for information about the estimated regularization path.")
    message("---") 
    }
