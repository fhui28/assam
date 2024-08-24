function() {
    formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula
    family <- tweedielogfam()
    resp = simdat$y
    data = covariate_dat
    offset = NULL
    trial_size = 1
    do_parallel = TRUE
    }
  	

.quadapprox_fn <- function(family, formula, resp, data, offset = NULL, 
                          trial_size = 1, 
                          do_parallel = TRUE,
                          num_nuisance_perspp = NULL,
                          start_archetype_labels = NULL,
                          start = NULL) {     
    
    num_spp <- ncol(resp)
     
    spp_qa_fn <- function(j) {
        tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
        
        if(is.null(start))
            use_start <- NULL
        if(!is.null(start))  { #' Assumes start_archetype_labels has also being supplied!
            use_start <- list(beta = c(start$spp_intercepts[j], start$betas[start_archetype_labels[j],]))
            if(num_nuisance_perspp > 0)
                use_start$betad <- log(start$spp_nuisance$dispersion[j])
            if(family$family[1] == "tweedie")
                use_start$psi <- qlogis(start$spp_nuisance$power[j]-1)
            }
        
        if(family$family %in% c("binomial")) {
            tmp_formula <- as.formula(paste("cbind(response, trial_size - response)", paste(as.character(formula), collapse = " ") ) )
            }
        new_offset <- offset[,j]
        #Hmat <- diag(control$ridge+1e-15, nrow = num_X)
         
        if(family$family %in% c("gaussian", "poisson", "Gamma", "binomial", "tweedie")) {
            fit0 <- glmmTMB(tmp_formula,
                            data = data.frame(response = resp[,j], data, trial_size = trial_size), 
                            #knots = knots, select = select, gamma = full_gamma[j]  
                            offset = new_offset, family = family, 
                            start = use_start)
            }
        if(family$family %in% c("negative.binomial")) {
            fit0 <- glmmTMB(tmp_formula, 
                            data = data.frame(response = resp[,j], data), 
                            #knots = knots, select = select, gamma = full_gamma[j]
                            offset = new_offset, family = "nbinom2",
                            start = use_start)
            }
        if(family$family %in% c("Beta")) {
            fit0 <- glmmTMB(tmp_formula, 
                            data = data.frame(response = resp[,j], data), 
                            #knots = knots, select = select, gamma = full_gamma[j]
                            offset = new_offset, family = "beta",
                            start = use_start)
            }
        fit0$parameters <- fit0$obj$env$last.par.best
        fit0$hessian <- fit0$obj$he(fit0$obj$env$last.par.best)
        return(fit0)
        }


    if(do_parallel)
        all_quadapprox <- foreach(j = 1:num_spp) %dopar% spp_qa_fn(j = j)          
    if(!do_parallel)
        all_quadapprox <- foreach(j = 1:num_spp) %do% spp_qa_fn(j = j)          

    out <- list(parameters = t(sapply(1:num_spp, function(k) all_quadapprox[[k]]$parameters)), 
                hessian = lapply(1:num_spp, function(k) all_quadapprox[[k]]$hessian))
     
    #rownames(out$parameters) <- colnames(resp)
    #names(out$hessian) <- colnames(resp)
    return(out)
    }

          