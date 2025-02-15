.quadapprox2_fn <- function(family, 
                            formula, 
                            resp, 
                            data, 
                            add_spatial = FALSE,
                            mesh = NULL, 
                            offset = NULL, 
                            trial_size = 1, 
                            do_parallel = TRUE,
                            control = NULL,
                            return_fits = TRUE) {     
    
    num_spp <- ncol(resp)
    if(add_spatial) {
        if(inherits(mesh, "sdmTMBmesh"))
            stop("If mesh is supplied for species-specific spatial fields, then the mesh argument must be an object class of \"sdmTMBmesh\".")
        }

    make_sdmTMBcontrol <- sdmTMB::sdmTMBcontrol(get_joint_precision = FALSE)
    if(!is.null(control$beta_lower))
        make_sdmTMBcontrol$lower <- list(b_j = control$beta_lower)
    if(!is.null(control$beta_upper))
        make_sdmTMBcontrol$upper <- list(b_j = control$beta_upper)
    if(!is.null(control$beta_lower) | !is.null(control$beta_upper))
        make_sdmTMBcontrol$newton_loops <- 0 # Needed to ensure constraints are actually respected; see (https://github.com/pbs-assess/sdmTMB/issues/394#issuecomment-2619995494)


    spp_qa_fn <- function(l) {
        tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
        
        if(family$family %in% c("binomial")) {
            tmp_formula <- as.formula(paste("cbind(response, trial_size - response)", paste(as.character(formula), collapse = " ") ) )
            }
        new_offset <- offset[,l]

        if(family$family %in% c("gaussian", "poisson", "Gamma", "binomial", "nbinom2", "tweedie")) {
            if(!add_spatial) {
                fit0 <- sdmTMB(tmp_formula,
                               data = data.frame(response = resp[,l], data, trial_size = trial_size), 
                               spatial = add_spatial,
                               offset = new_offset, 
                               family = family,
                               control = make_sdmTMBcontrol)
                }
            if(add_spatial) {
                fit0 <- sdmTMB(tmp_formula,
                               data = data.frame(response = resp[,l], data, trial_size = trial_size), 
                               spatial = add_spatial,
                               mesh =  mesh,
                               offset = new_offset, 
                               family = family,
                               control = make_sdmTMBcontrol)
                }
            }
        if(family$family %in% c("Beta")) {
            if(!add_spatial) {
                fit0 <- sdmTMB(tmp_formula, 
                               data = data.frame(response = resp[,l], data), 
                               spatial = add_spatial,
                               offset = new_offset, 
                               family = sdmTMB::Beta(),
                               control = make_sdmTMBcontrol)
                }
            if(add_spatial) {
                fit0 <- sdmTMB(tmp_formula, 
                               data = data.frame(response = resp[,l], data), 
                               spatial = add_spatial,
                               mesh =  mesh,
                               offset = new_offset, 
                               family = sdmTMB::Beta(),
                               control = make_sdmTMBcontrol)
                }
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
    
    class(out) <- "assam_quadapprox"
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

