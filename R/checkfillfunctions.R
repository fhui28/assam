.fill_control <- function(control) {
    if(is.null(control$max.iter))
        control$max.iter <- 500
    if(is.null(control$tol))
        control$tol <- 1e-5
    if(is.null(control$trace))
        control$trace <- FALSE
    
    return(control)
    }


.fill_bootstrap_control <- function(control) {
    if(is.null(control$num_boot))
        control$num_boot <- 100
    if(is.null(control$ci_alpha))
        control$ci_alpha <- 0.05
    if(is.null(control$seed))
        control$seed <- NULL
    if(is.null(control$ci_type))
        control$ci_type <- "percentile"
    
    return(control)
    }


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


.check_X_formula <- function(formula, data) {
    formulaX <- as.formula(formula)
    
    termsinformula <- as.character(formula)
    if(length(termsinformula) == 3)
        termsinformula <- termsinformula[-2]
    formula <- formula(paste(termsinformula, collapse = " "))
    
    return(formula)
}     

