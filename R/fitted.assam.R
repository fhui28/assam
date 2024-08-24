#' @title Extract fitted values from an approximate and scalable SAM object
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Extracts fitted mean values from a \code{assam} object.
#' 
#' @param object An object of class \code{assam}.
#' @param type The type of fitted value to obtain. Currently one of three options are possible: \code{all}, which for each species returns the fitted values for all archetypes; \code{max}, which for each species returns the fitted values corresponding to the most likely archetype the species belongs to (as based on the posterior probabilities); \code{mean}, which for each species returns a weighted mean of the fitted values from each archetype, where the weights are given by the posterior probabilities.
#' @param ... Not used.
#' 
#' @details 
#' To clarify, the returned fitted values are on the response scale i.e., a matrix of the estimated means \eqn{\hat{\mu}_{ij}} after model fitting.
#' 
#' @return A matrix of fitted values.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>

#' @export
#' @md

fitted.assam <- function(object, type = "all", ...) {
    if(!inherits(object, "assam")) 
        stop("`object' is not of class \"assam\"")
    
    type <- match.arg(type, choices = c("all", "max", "mean"))
    
    if(type == "all")
        return(object$family$linkinv(object$linear_predictor))
    
    if(type == "max") {
        num_spp <- nrow(object$posterior_probability)
        get_eta <- sapply(1:num_spp, function(j) object$linear_predictor[,j,which.max(object$posterior_probability[j,])])
        return(object$family$linkinv(get_eta))
        }
    
    if(type == "mean") {
        num_spp <- nrow(object$posterior_probability)
        num_archetype <- ncol(object$posterior_probability)
        num_units <- dim(object$linear_predictor)[1]
        get_eta <- object$family$linkinv(object$linear_predictor)
        get_eta <- sapply(1:num_spp, function(j) rowSums(matrix(object$posterior_probability[j,], nrow = num_units, ncol = num_archetype, byrow = TRUE)*get_eta[,j,]))
        return(get_eta)
        }
    
    }


