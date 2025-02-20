#' @noRd
#' @noMd
.fill_control <- function(control) {
    if(is.null(control$max_iter))
        control$max_iter <- 500
    if(is.null(control$tol))
        control$tol <- 1e-5
    if(is.null(control$temper_prob))
        control$temper_prob <- 0.8
    if(is.null(control$trace))
        control$trace <- FALSE
    if(is.null(control$beta_lower))
        control$beta_lower <- NULL
    if(is.null(control$beta_upper))
        control$beta_upper <- NULL
    if(is.null(control$betamatrix_select))
        control$betamatrix_select <- NULL
    
    return(control)
    }


#' @noRd
#' @noMd
.fill_bootstrap_control <- function(control) {
    if(is.null(control$num_boot))
        control$num_boot <- 100
    if(is.null(control$method))
        control$method <- "full_bootstrap"
    if(is.null(control$ci_alpha))
        control$ci_alpha <- 0.05
    if(is.null(control$seed))
        control$seed <- NULL
    if(is.null(control$ci_type))
        control$ci_type <- "percentile"
    
    return(control)
    }


#' @noRd
#' @noMd
.check_betamatrix_selection <- function(control, num_archetypes, num_X) {
    if(!is.null(control$betamatrix_selection)) {
        if(!is.null(control$beta_lower) | !is.null(control$beta_upper))
            stop("If control$betamatrix_selection is supplied, then control$beta_lower and control$beta_upper can not also be supplied.")
    }
    
    if(!is.null(control$betamatrix_selection)) {
        if(!is.matrix(control$betamatrix_selection))
            stop("control$betamatrix_selection must be a matrix.")
        if(nrow(control$betamatrix_selection) != num_archetypes)
            stop("The number of rows in control$betamatrix_selection must be equal to num_archetypes.")
        if(ncol(control$betamatrix_selection) != (num_X-1))
            stop("The number of columns in control$betamatrix_selection must be equal to the number of columns in the model matrix induced by formula and data, minus one.")
        if(!all(control$betamatrix_selection %in% 0:1))
            stop("control$betamatrix_selection must be a binary matrix.")
        }
    }



#' @noRd
#' @noMd
.check_offset <- function(offset = NULL, y) {
    if(!is.null(offset)) { 
        if(!is.matrix(offset)) 
            stop("offset should be a matrix with the same dimensions as y.")
        if(nrow(offset) != nrow(y)) 
            stop("offset should be a matrix with the same dimensions as y.")
        if(ncol(offset) != ncol(y)) 
            stop("offset should be a matrix with the same dimensions as y.")
        } 
    }


#' @noRd
#' @noMd
.check_X_formula <- function(formula, data) {
    formulaX <- as.formula(formula)
    
    termsinformula <- as.character(formula)
    if(length(termsinformula) == 3)
        termsinformula <- termsinformula[-2]
    formula <- formula(paste(termsinformula, collapse = " "))
    
    return(formula)
    }     

#' @noRd
#' @noMd
.check_spp_spatial_parameters <- function(spp_spatial_range, spp_spatial_sd, add_spatial, spp_effects) {
    if(!is.null(spp_spatial_range)) {
        if(is.null(spp_spatial_sd) | !add_spatial)
            stop("If species-specific spatial fields are to be included, then mesh, spp_spatial_range, spp_spatial_sd must all be supplied.")
        if(length(spp_spatial_range) != nrow(spp_effects))
            stop("Both spp_spatial_range and spp_spatial_sd must be vectors of length equal to nrow(spp_effects).")
        if(!all(spp_spatial_range > 0))
            stop("spp_spatial_range must be vectors of strictly positive elements.")
        }
    
    if(!is.null(spp_spatial_sd)) {
        if(is.null(spp_spatial_range) | !add_spatial)
            stop("If species-specific spatial fields are to be included, then mesh, spp_spatial_range, spp_spatial_sd must all be supplied.")
        if(length(spp_spatial_sd) != nrow(spp_effects))
            stop("Both spp_spatial_range and spp_spatial_sd must be vectors of length equal to now(spp_effects).")
        if(!all(spp_spatial_sd > 0))
            stop("spp_spatial_sd must be vectors of strictly positive elements.")
        }
    
    }


