.quadapprox2_fn <- function(family, 
                            formula, 
                            resp, 
                            data, 
                            mesh = NULL, 
                            offset = NULL, 
                            trial_size = 1, 
                            do_parallel = TRUE,
                            return_fits = TRUE) {     
    
    num_spp <- ncol(resp)
    add_spatial <- FALSE
    if(!is.null(mesh))
        add_spatial <- TRUE
    if(!add_spatial)
        mesh <- sdmTMB::make_mesh(data = data, xy_cols = colnames(data)[1:2], n_knots = 5) # Make any old mesh as sdmTMB must supplied this
    
    spp_qa_fn <- function(l) {
        tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
        
        if(family$family %in% c("binomial")) {
            tmp_formula <- as.formula(paste("cbind(response, trial_size - response)", paste(as.character(formula), collapse = " ") ) )
            }
        new_offset <- offset[,l]
        #Hmat <- diag(control$ridge+1e-15, nrow = num_X)
        
        if(family$family %in% c("gaussian", "poisson", "Gamma", "binomial")) {
            fit0 <- sdmTMB(tmp_formula,
                           data = data.frame(response = resp[,l], data, trial_size = trial_size), 
                           spatial = add_spatial,
                           mesh =  mesh,
                           offset = new_offset, 
                           family = family)
                }
        if(family$family %in% c("tweedie")) {
            fit0 <- sdmTMB(tmp_formula,
                           data = data.frame(response = resp[,l], data, trial_size = trial_size), 
                           spatial = add_spatial,
                           mesh =  mesh,
                           offset = new_offset, 
                           family = sdmTMB::tweedie())
            }
        if(family$family %in% c("nbinom2")) {
            fit0 <- sdmTMB(tmp_formula, 
                           data = data.frame(response = resp[,l], data), 
                           spatial = add_spatial,
                           mesh =  mesh,
                           offset = new_offset, 
                           family = sdmTMB::nbinom2())
            }
        if(family$family %in% c("Beta")) {
            fit0 <- sdmTMB(tmp_formula, 
                           data = data.frame(response = resp[,l], data), 
                           spatial = add_spatial,
                           mesh =  mesh,
                           offset = new_offset, 
                           family = sdmTMB::Beta())
            }
        
        fit0$parameters <- fit0$sd_report$par.fixed
        fit0$hessian <- solve(fit0$sd_report$cov.fixed)
        return(fit0)
        }
    
    if(do_parallel)
        all_quadapprox <- foreach(l = 1:num_spp) %dopar% spp_qa_fn(l = l) 
    if(!do_parallel)
        all_quadapprox <- foreach(l = 1:num_spp) %do% spp_qa_fn(l = l)          

    out <- list(parameters = t(sapply(1:num_spp, function(k) all_quadapprox[[k]]$parameters)), 
                hessian = lapply(1:num_spp, function(k) all_quadapprox[[k]]$hessian))
    if(return_fits)
        out$sdmTMB_fits <- all_quadapprox
    
    return(out)
    }


#' Old function using glmmTMB
    # .quadapprox_fn <- function(family, formula, resp, data, offset = NULL, 
    #                       trial_size = 1, 
    #                       do_parallel = TRUE,
    #                       num_nuisance_perspp = NULL,
    #                       start_archetype_labels = NULL,
    #                       start = NULL) {     
    # 
    # num_spp <- ncol(resp)
    #  
    # spp_qa_fn <- function(l) {
    #     tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
    #     
    #     if(is.null(start))
    #         use_start <- NULL
    #     if(!is.null(start))  { #' Assumes start_archetype_labels has also being supplied!
    #         use_start <- list(beta = c(start$spp_intercepts[l], start$betas[start_archetype_labels[l],]))
    #         if(num_nuisance_perspp > 0)
    #             use_start$betad <- log(start$spp_nuisance$dispersion[l])
    #         if(family$family[1] == "tweedie")
    #             use_start$psi <- qlogis(start$spp_nuisance$power[l]-1)
    #         }
    #     
    #     if(family$family %in% c("binomial")) {
    #         tmp_formula <- as.formula(paste("cbind(response, trial_size - response)", paste(as.character(formula), collapse = " ") ) )
    #         }
    #     new_offset <- offset[,l]
    #     #Hmat <- diag(control$ridge+1e-15, nrow = num_X)
    #      
    #     if(family$family %in% c("gaussian", "poisson", "Gamma", "binomial", "tweedie")) {
    #         fit0 <- glmmTMB(tmp_formula,
    #                         data = data.frame(response = resp[,l], data, trial_size = trial_size), 
    #                         #knots = knots, select = select, gamma = full_gamma[j]  
    #                         offset = new_offset, family = family, 
    #                         start = use_start)
    #         }
    #     if(family$family %in% c("nbinom2")) {
    #         fit0 <- glmmTMB(tmp_formula, 
    #                         data = data.frame(response = resp[,l], data), 
    #                         #knots = knots, select = select, gamma = full_gamma[j]
    #                         offset = new_offset, family = glmmTMB::nbinom2,
    #                         start = use_start)
    #         }
    #     if(family$family %in% c("Beta")) {
    #         fit0 <- glmmTMB(tmp_formula, 
    #                         data = data.frame(response = resp[,l], data), 
    #                         #knots = knots, select = select, gamma = full_gamma[j]
    #                         offset = new_offset, family = glmmTMB::beta_family,
    #                         start = use_start)
    #         }
    #     fit0$parameters <- fit0$obj$env$last.par.best
    #     fit0$hessian <- fit0$obj$he(fit0$obj$env$last.par.best)
    #     return(fit0)
    #     }
    # 
    # 
    # if(do_parallel)
    #     all_quadapprox <- foreach(l = 1:num_spp) %dopar% spp_qa_fn(l = l)          
    # if(!do_parallel)
    #     all_quadapprox <- foreach(l = 1:num_spp) %do% spp_qa_fn(l = l)          
    # 
    # out <- list(parameters = t(sapply(1:num_spp, function(k) all_quadapprox[[k]]$parameters)), 
    #             hessian = lapply(1:num_spp, function(k) all_quadapprox[[k]]$hessian))
    #  
    # #rownames(out$parameters) <- colnames(resp)
    # #names(out$hessian) <- colnames(resp)
    # return(out)
    # }

