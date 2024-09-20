#' @title Extract residuals from an approximate and scalable SAM object
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Calculate either probability integral transform (PIT) residuals and Dunn-Smyth residuals residuals from a fitted \code{assam} object.
#' 
#' @param object An object of class \code{assam}.
#' @param type The type of residuals which should be returned. The options currently available are "PIT" and "dunnsmyth".
#' @param seed This can be used set the seed when constructing the PIT and Dunn-Smyth residuals, which for discrete responses involve some degree of jittering.  
#' @param ... Not used.
#' 
#' @details 
#' Of the two options available \code{type = "PIT"} returns the probability integral transform residuals that are also used in [DHARMa::simulateResiduals()] and [mpcmp::rPIT()], among other packages. If the (estimated) model is correct, then these residuals should behave as random variables from a standard uniform distribution [Dunn and Smyth, 1996](https://doi.org/10.2307/1390802). Note there is a level of jittering used in producing the PIT residuals. On the other hand, \code{type = "dunnsmyth"} returns the Dunn-Smyth residuals that are also used in [boral::ds.residuals()] and[gllvm::residuals.gllvm()], among other packages. If the (estimated) model is correct, then these residuals should behave as random variables from a standard normal distribution (Dunn and Smyth, 1996). Note there is a level of jittering used in producing the Dunn-Smyth residuals.
#' 
#' Both PIT and Dunn-Smyth residuals were adapted by [Dunstan et al., (2013)](https://doi.org/10.1007/s13253-013-0146-x) to species archetype models, and we follow their developments in building this function.
#' 
#' @return A matrix of residuals.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' 
#' @examples
#' \dontrun{
#' #' Please see the help file for assam for example.
#' }
#' 
#' @export
#' @importFrom stats runif qnorm pbeta pbinom pgamma plogis pnorm ppois pnbinom
#' @importFrom tweedie ptweedie 
#' @md


residuals.assam <- function(object, type = "dunnsmyth", seed = NULL, ...) {
    if(!inherits(object, "assam")) 
        stop("`object' is not of class \"assam\"")
    
    y <- sapply(object$sdmTMB_fits, function(x) x$data$response)
    
    type <- match.arg(type, choices = c("dunnsmyth", "PIT"))
    get_fitted <- fitted.assam(object, type = "all")
    num_archetype <- ncol(object$posterior_probability)
    out <- matrix(0, nrow = nrow(y), ncol = ncol(y)) 
    
    if(object$family$family[1] %in% c("binomial", "poisson", "nbinom2", "tweedie"))
        a <- b <- matrix(0, nrow = nrow(y), ncol = ncol(y)) 
    
    set.seed(seed)
    
    for(k0 in 1:num_archetype) {
        if(object$family$family[1] %in% c("Beta")) 
            out <- out + object$mixture_proportion[k0]*pbeta(y, 
                                                             shape1 = matrix(object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE) * get_fitted[,,k0],
                                                             shape2 = matrix(object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE) * (1-get_fitted[,,k0]))
        
        if(object$family$family[1] %in% c("gaussian")) 
            out <- out + object$mixture_proportion[k0]*pnorm(y, mean = get_fitted[,,k0], sd = matrix(object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
        
        if(object$family$family[1] %in% c("Gamma")) 
            out <- out + object$mixture_proportion[k0]*pgamma(y, 
                                                              scale = matrix(object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE) * get_fitted[,,k0], 
                                                              shape = matrix(1/object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
        
        if(object$family$family[1] %in% c("poisson")) {
            a <- a + object$mixture_proportion[k0]*ppois(y-1, lambda = get_fitted[,,k0])
            b <- b + object$mixture_proportion[k0]*ppois(y, lambda = get_fitted[,,k0])
            }
        
        if(object$family$family[1] %in% c("binomial")) {
            a <- a + object$mixture_proportion[k0]*pbinom(y-1, size = object$trial_size, prob = get_fitted[,,k0])
            b <- b + object$mixture_proportion[k0]*pbinom(y, size = object$trial_size, prob = get_fitted[,,k0])
            }
        
        if(object$family$family[1] %in% c("nbinom2")) {
            a <- a + object$mixture_proportion[k0]*pnbinom(y-1, mu = get_fitted[,,k0], size = matrix(object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
            b <- b + object$mixture_proportion[k0]*pnbinom(y, mu = get_fitted[,,k0], size = matrix(object$spp_nuisance$dispersion, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
            }
        
        if(object$family$family[1] %in% c("tweedie")) {
            for(j in 1:ncol(y)) {
                a[,j] <- a[,j] + object$mixture_proportion[k0]*tweedie::ptweedie(y[,j], 
                                                                                 mu = get_fitted[,j,k0],
                                                                                 phi = object$spp_nuisance$dispersion[j],
                                                                                 power = object$spp_nuisance$power[j])
                b[,j] <- b[,j] + object$mixture_proportion[k0]*tweedie::ptweedie(y[,j], 
                                                                                 mu = get_fitted[,j,k0],
                                                                                 phi = object$spp_nuisance$dispersion[j],
                                                                                 power = object$spp_nuisance$power[j])
            }
            a[y == 0] <- 0
            }
        
        }

    if(object$family$family[1] %in% c("binomial","poisson","nbinom2","tweedie")) {
        out <- matrix(runif(length(y), min = a, max = b), nrow = nrow(y), ncol = ncol(y))
        out[out < 1e-12] <- 1e-12        
        out[out > (1-1e-12)] <- 1-1e-12        
        }
    
    if(type == "dunnsmyth")
        out <- qnorm(out)
    
    set.seed(NULL)
    return(out)
    }



