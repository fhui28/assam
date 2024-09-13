function() {
    
    reprex({
    rm(list = ls())
    library(foreach)
    library(doParallel)
    registerDoParallel(cores = 8)
    
    #' # Simulate some multivariate abundance data
    set.seed(092024)
    covariate_dat <- data.frame(x = runif(1000), 
                                y = runif(1000),
                                cov1 = rnorm(1000),
                                cov2 = rnorm(1000),
                                cov3 = rnorm(1000),
                                cov4 = rbinom(1000, size = 1, prob = 0.5),
                                cov5 = runif(1000))
    
    resp_dat <- matrix(rnorm(1000*50), nrow = 1000)
    colnames(resp_dat) <- paste("species", 1:50)
    
    
    #' # Fit stacked non-spatial models
    #' Everything seems to work as is with %dopar% actually working in parallel and taking less time
    testfn <- function(l) {
        fit0_orig <- sdmTMB::sdmTMB(response ~ cov1 + cov2 + cov3 + cov4 + cov5,
                                    data = data.frame(response = resp_dat[,l], covariate_dat), 
                                    spatial = FALSE,
                                    mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 25),
                                    family = gaussian())
        }
    
    system.time(foreach(l = 1:50) %do% testfn(l))
    system.time(foreach(l = 1:50) %dopar% testfn(l))
    
    
    #' # Fit stacked spatial models 
    #' Now %dopar% actually takes longer and when I have a look over at cores being like it looks like both are using multiple cores.
    testfn <- function(l) {
        fit0_orig <- sdmTMB::sdmTMB(response ~ cov1 + cov2 + cov3 + cov4 + cov5,
                                    data = data.frame(response = resp_dat[,l], covariate_dat), 
                                    spatial = TRUE,
                                    mesh = sdmTMB::make_mesh(covariate_dat, xy_cols = c("x", "y"), n_knots = 25),
                                    family = gaussian())
        }
    
    system.time(foreach(l = 1:50) %do% testfn(l))
    system.time(foreach(l = 1:50) %dopar% testfn(l))
    
    
    sessioninfo::session_info()
    })
}

