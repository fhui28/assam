rm(list = ls())
library(tidyverse)
library(Matrix)
library(foreach)
library(abind)
library(label.switching)
library(GGally)
library(doParallel)
library(glmmTMB)
library(cluster)
library(mvtnorm)
here::i_am("testv0.R")
library(here)

source(here("assam.R"))
source(here("betalogitfam.R"))
source(here("create_samlife.R"))
source(here("checkfillfunctions.R"))
source(here("fitted.assam.R"))
source(here("nb2.R"))
source(here("plot.assam.R"))
source(here("predict.assam.R"))
source(here("print.assam.R"))
source(here("quadapprox_fn.R"))
source(here("residuals.assam.R"))
source(here("simulate.assam.R"))
source(here("summary.assam.R"))
source(here("tweedielogfam.R"))

##----------------------
## Generate some multivariate abundance data from a SAM
##----------------------
set.seed(092024)

num_X <- 10
num_units <- 1000
num_spp <- 80
num_archetype <- 5
H <- outer(1:num_X, 1:num_X, "-")
H <- 0.5^abs(H)
covariate_dat <- rmvnorm(num_units, sigma = H) %>% 
    as.data.frame %>% 
    rename_with(., .fn = function(x) paste0("covariate", x))
rm(H)

#true_betas <- rnorm(num_archetype*num_X,0,sd=0.5) %>% 
#     matrix(nrow = num_archetype)
true_betas <- runif(num_archetype * num_X, -1, 1) %>% 
    matrix(nrow = num_archetype)
true_intercepts <- runif(num_spp, -3, 0)  
true_dispparam <- 1/runif(num_spp, 0, 5) 
true_powerparam <- runif(num_spp, 1.4, 1.8)
true_mixprop <- c(0.2, 0.25, 0.3, 0.1, 0.15)


simdat <- create_samlife(family = nb2(), 
                         formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula, 
                         data = covariate_dat, 
                         betas = true_betas, 
                         spp_intercept = true_intercepts, 
                         spp_dispparam = true_dispparam, 
                         spp_powerparam = true_powerparam, 
                         mixture_proportion = true_mixprop)


##----------------------
## Fit asSAM and assess results and functions
##----------------------
tic <- proc.time()
testfit <- assam(y = simdat$y,
                 formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                 data = covariate_dat,
                 family = nb2(),
                 num_archetypes = num_archetype,
                 uncertainty_quantification = FALSE,
                 num_cores = 8)
toc <- proc.time()
toc-tic
# user  system elapsed 
# 25.101   1.479   4.576 


tic <- proc.time()
testfit <- assam(y = simdat$y,
                 formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) %>% as.formula,
                 data = covariate_dat,
                 family = nb2(),
                 num_archetypes = num_archetype,
                 uncertainty_quantification = TRUE,
                 num_cores = 8)
toc <- proc.time()
toc-tic
# user   system  elapsed 
# 2991.192  541.550  537.708 

plot(true_intercepts, testfit$spp_intercepts); abline(0,1)
plot(true_dispparam, testfit$spp_nuisance$dispersion, log = "xy"); abline(0,1)
plot(true_powerparam, testfit$spp_nuisance$power); abline(0,1)
rbind(true_betas, testfit$betas) %>% 
    t %>% 
    as.data.frame %>%
    ggpairs
table(simdat$archetype_label, 
      apply(testfit$posterior_probability, 1, which.max))


testfit
summary(testfit)

fitted(testfit)
simulate(testfit, data = covariate_dat)
residuals(testfit, y = simdat$y, type = "dunnsmyth")
plot(testfit, y = simdat$y, transform_fitted_values = TRUE)
plot(simdat$linear_predictor, fitted(testfit, type = "mean") %>% log)

predict(testfit, 
        newdata = covariate_dat, 
        type = "species_mean") - fitted(testfit, type = "mean")

predict(testfit, 
        newdata = covariate_dat, 
        type = "archetype",
        se_fit = TRUE) 

predict(testfit, 
        newdata = covariate_dat, 
        type = "species_max",
        num_cores = 8,
        se_fit = TRUE) 



#----------------------
## Compare to Skip's ecomix package
#' There is a known issue with dispersion parameters in negative binomial SAMs: see (https://github.com/skiptoniam/ecomix/issues/23)
##----------------------
#devtools::install_github('skiptoniam/ecomix')
library(ecomix)

archetype_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:num_spp), collapse = ','),
                                          ")", paste("~ ", paste0(colnames(covariate_dat), collapse = "+"))))

tic <- proc.time()
testfit_ecomix <- species_mix(archetype_formula = archetype_form,
                              species_formula = ~ 1,
                              data = data.frame(simdat$y, const = 1, covariate_dat),
                              family = 'bernoulli', 
                              power = 1.6,
                              nArchetypes = 5)
toc <- proc.time()
# user  system elapsed 
# 13.526   0.041  13.568 


ggpairs(data.frame(true = true_intercepts, ecomix = testfit_ecomix$alpha, assam = testfit$spp_intercepts))
ggpairs(data.frame(true = true_dispparam, ecomix = 1/testfit_ecomix$theta, assam = testfit$spp_nuisance$dispersion)) +
    scale_x_continuous(trans = "log1p") +
    scale_y_continuous(trans = "log1p")
data.frame(true = true_mixprop, ecomix = testfit_ecomix$pi, assam = testfit$mixture_proportion)

rbind(true_betas, coef(testfit_ecomix)$beta) %>% 
    t %>% 
    as.data.frame %>%
    ggpairs
table(simdat$archetype_label, 
      apply(testfit_ecomix$tau, 1, which.max))

#' ## Reorder rows of estimated archetypal coefficient matrices to get it as close as possible to the order of true coefficients
norm(true_betas - coef(testfit_ecomix)$beta[table(simdat$archetype_label, 
                                             apply(testfit_ecomix$tau, 1, which.max)) %>% 
                                           apply(., 1, which.max) %>% 
                                           as.integer(),])
norm(true_betas - testfit$betas[table(simdat$archetype_label, 
                                      apply(testfit$posterior_probability, 1, which.max)) %>% 
                                    apply(., 1, which.max) %>% 
                                    as.integer(),])

tic <- proc.time()
testfit_ecomix_boot <- vcov(testfit_ecomix,
                            nboot = 100,
                            method = "BayesBoot",
                            mc.cores = 8)
toc <- proc.time()
toc-tic
# user   system  elapsed 
# 1098.582    4.297  189.438 

testfit_ecomix$vcov <- testfit_ecomix_boot

summary(testfit_ecomix)
