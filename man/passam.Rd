% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/passam.R
\name{passam}
\alias{passam}
\title{Penalized approximate and scalable species archetype models (pasSAMs)}
\usage{
passam(
  y,
  formula,
  data,
  which_spp_effects = 1,
  family,
  offset = NULL,
  trial_size = 1,
  num_archetypes,
  mesh = NULL,
  do_parallel = TRUE,
  num_cores = NULL,
  selection_on = "betas",
  supply_quadapprox = NULL,
  nlambda = 100,
  lambda_min_ratio = 0.001,
  lambda = NULL,
  control = list(max_iter = 500, tol = 1e-04, temper_prob = 0.7, trace = FALSE),
  beta_selection_control = list(min_df = 0, max_iter = 100, eps = 1e-04, round_eps =
    1e-05)
)
}
\arguments{
\item{y}{A multivariate abundance response matrix.}

\item{formula}{An object of class "formula", which represents a symbolic description of the full model matrix to be created (based on using this argument along with the \code{data} argument). \emph{Note there should be nothing on the left hand side of the "~".} Currently, smooth terms are not permitted.}

\item{data}{A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument).}

\item{which_spp_effects}{A vector identifying which columns of the model matrix induced by \code{formula} and \code{data} should be treated as species-specific effects. Default to 1, meaning only the first column i.e., the intercept, is species-specific.}

\item{family}{a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.}

\item{offset}{A matrix of offset terms, of the same dimension as \code{y}.}

\item{trial_size}{Trial sizes to use for binomial distribution. This should equal to a scalar.}

\item{num_archetypes}{Number of archetypes (clusters) to assume in the asSAM. Note when \code{selection_on = "mixing_proportions"}, then this should be an upper bound i.e., a maximum number of archetypes the practitioner wants in their asSAM.}

\item{mesh}{Output from \code{\link[sdmTMB:make_mesh]{sdmTMB::make_mesh()}}, used for adding species-specific spatial fields to the linear predictor.}

\item{do_parallel}{Should parallel computing be used to fit the asSAM. Defaults to \code{FALSE}.}

\item{num_cores}{If \code{do_parallel = TRUE}, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to \code{parallel::detectCores() - 2}.}

\item{selection_on}{Should the penalty be applied, and thus sparsity be encouraged, on the archetypal regression coefficients or the mixing proportions? Current choices are \code{"betas"}, which is former, or \code{"mixing_proportions"}, which is the latter. Defaults to \code{"betas"}.}

\item{supply_quadapprox}{An object of class \code{assam_quadapprox}, which is (mostly likely) obtained as a consequence of running an initial fit using \code{\link[=assam]{assam()}} with \code{do_assam_fit = FALSE}.}

\item{nlambda}{The number of tuning parameters values to consider when forming the regularization path.}

\item{lambda_min_ratio}{The smallest value for the tuning parameter lambda, as a fraction of the maximum internally-derived lambda value.}

\item{lambda}{A user-supplied tuning parameter sequence. Note this is usually not supplied, as it is standard to have the function itself compute its own lambda sequence based on \code{nlambda} and \code{lambda_min_ratio}.}

\item{control}{A list containing the following elements:
\describe{
\item{max_iter:}{the maximum number of iterations in the EM algorithm. Usually convergence is quite quick e.g., less than 20 iterations.}
\item{tol:}{the convergence criterion; the difference in the log-likelihood value of the asSAM from successive iterations must be smaller than this value.}
\item{temper_prob:}{in the iteration of the EM algorithm, posterior probabilities from the E-step are "tempered" or push away from the 0/1 boundary. This is often useful to get the EM algorithm moving initially.}
\item{trace:}{controls if messages printed as part of the estimation process to reflect progress.}
}}

\item{beta_selection_control}{A list containing the following elements to control the broken adaptive ridge (BAR) penalty for variable selection on the archetypal regression coefficients:
\describe{
\item{min_df:}{The minimum number of non-zero archetypal regression coefficients allowed in the asSAM. This is useful to supply if, when the function tries to find an appropriate lambda sequence, the largest value of lambda is determined such that the number of archetypal regression coefficients estimated as non-zero is equal to or below this number. Defaults to zero.}
\item{max_iter:}{the maximum number of iterations in the BAR optimization part of the EM algorithm.}
\item{eps:}{the convergence criterion; the norm of the difference between all estimated parameters from successive iterations must be smaller than this value.}
\item{round_eps:}{a tolerance to round values to zero. The technically not needed as the BAR penalty will produce exactly zero estimates up to machine error, but is included anyway, but is included anyway.}
}}
}
\value{
An object of class \code{passam} with the following elements (as appropriate, and not necessarily in the order below):
\item{call:}{The function call.}
\item{formula:}{Same as input argument.}
\item{family:}{Same as input argument.}
\item{num_nuisance_perspp:}{The number of "nuisance" parameters per species. For example, if \code{family = binomial()} and species-specific spatial fields are not included, then there are no nuisance parameters. If \code{family = nbinom2()} and species-specific spatial fields are included, say, then there are three nuisance parameters per species.}
\item{trial_size:}{Same as input argument.}
\item{offset:}{Same as input argument.}
\item{add_spatial:}{Were species-specific spatial fields included?}
\item{mesh:}{Same as input argument.}
\item{num_archetypes:}{Same as input argument.}
\item{nlambda:}{Same as input argument.}
\item{lambda_min_ratio:}{Same as input argument.}
\item{lambda:}{The actual sequence of tuning parameter, lambda, values used.}
\item{parameters_df:}{The number of non-zero archetypal regression coefficients (if \code{selection_on = "betas"}) or number of estimated non-zero parameters (mixing proportions and archetypal regression coefficients, if \code{selection_on = "mixing_proportions"}) at each value of lambda. Note in both cases, the number of species-specific parameters is not counted since this does not change irrespective of the value of lambda.}
\item{parameters_path:}{The estimated archetypal regression coefficients (if \code{selection_on = "betas"}) or the estimated mixing proportions (if \code{selection_on = "mixing_proportions"}) at each value of lambda.
For the former, this takes the form of a three-dimensional array, where the third dimension corresponds to the lambda e.g., \code{betas_path[,,1]} corresponds to the estimated archetypal regression coefficients at the first value of lambda.
For the latter, this takes the form of a matrix where second dimension corresponds to the lambda e.g.,  \code{mixing_proportions_path[,1]} corresponds to the estimated mixing proportions at the first value of lambda}.
\item{logL:}{Estimated log-likelihood value of the asSAM i.e. the value of the approximated log-likelihood function, at each value of lambda.}
\item{AIC:}{The Akaike Information Criterion at each value of lambda. The AIC is defined as \eqn{-2\log(L) + 2df}, where \eqn{L} is the log-likelihood value and \eqn{df} is the given by \code{betas_df/mixing_proportions_df}.}
\item{BIC:}{The Bayesian Information Criterion at each value of lambda. The BIC is defined as \eqn{-2\log(L) + \log(M)df}, where \eqn{M} is the number of species.}
\item{BIC2:}{The Bayesian Information Criterion, version 2, at each value of lambda. The BIC is defined as \eqn{-2\log(L) + \log(MN)df}, where \eqn{M} is the number of species and \eqn{N} is the number of observational units. Compared to BIC, BIC2 uses a more severe model complexity penalty and thus is more likely to selection a sparse asSAM.}
\item{regularization_frame:}{A data frame combining some of the above information for easy "digestion".}
\item{betas_path/spp_effects_path/spp_nuisance_path/mixing_proportions_path/posterior_probability_path/logL_path}{Additional regularization path information for other parameters. This set of output can be safely ignored unless the user is interested in the path of these parameters, or wants to use some of these as warm starts for downstream applications of the \code{assam} function; see the example below.}
\item{control:}{Same as input argument.}
\item{beta_selection_control:}{Same as input argument.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Fits penalized approximate and scalable species archetype modeling (asSAMs) for model-based clustering of species based on their environmental response, into a small number of so-called archetypal responses. This is a modification of the main \code{\link[=assam]{assam()}} function, where a penalty is augmented to the approximate log-likelihood function to encourage sparsity in either the archetypal regression coefficients or the mixing proportions.
}
\details{
For the purposes of the package, the SAM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,M}, conditional on the species belong to archetype \eqn{k},

\deqn{g(\mu_{ij}) = \eta_{ij} = u_i^\top\alpha_j + x_i^\top\beta_k,}

where \eqn{g(.)} is a known link function, \eqn{u_i^\top\alpha_j} corresponds to a component that is to kept species-specific e.g., species-specific intercept, \eqn{x_i^\top\beta_k}  denotes the component corresponding to effect of archetypal response \eqn{k}. Additionally, species-specific spatial fields can be included in the linear predictor e.g., to account for residual spatial correlation above and beyond that explained by the archetypal responses. Conditional on the mean model above, the \eqn{y_{ij}} are assumed to be come from some response distribution using the additional dispersion and power parameters as appropriate. We refer the reader to the main \code{\link[=assam]{assam()}} function for more details and references for SAMs.

This function is specifically designed for performing variable selection, and building an associated regularization path, on either:
\enumerate{
\item the elements of \eqn{\beta_k} via the broken adaptive ridge \href{https://www.sciencedirect.com/science/article/pii/S0047259X17305067}{BAR} penalty. The BAR penalty can be interpreted as a kind of approximation to the \eqn{L_0} penalty, and encourages sparsity in the archetypal regression coefficients e.g., to uncover what covariates are informative for each of the archetypal responses;
\item the mixing proportions \eqn{\pi_k} via the log penalty of \href{https://www3.stat.sinica.edu.tw/statistica/J27N1/J27N17/J27N17.html}{equation (2.3) in}. The log penalty can be interpreted as like a \href{https://en.wikipedia.org/wiki/Lasso_(statistics)}{Lasso or \eqn{L_1} penalty} except more "steep", and encourages sparsity in the mixing proportions i.e., sends one or more of the \eqn{\pi_k}'s to zero and can be used to decide how many archetypes are needed in the SAM. Note when a mixing proportion is set to zero, the corresponding archetype is removed from the SAM.
}

After building a regularization path using this function, one can examine the path and select the tuning parameter minimizing one of the provided information criterion, say, to pass to \code{\link[=assam]{assam()}} function for fitting the final, sparse asSAM.
}
\section{Distributions and parallelization}{

Please see the \code{\link[=assam]{assam()}} function for more details on the distributions currently supported by the package, and how to properly parallelize the fitting process.
}

\examples{
\dontrun{
##----------------------
# Generate some multivariate abundance (count) data from a sparse SAM
##----------------------
library(tidyverse)
library(mvtnorm)
library(GGally)
library(doParallel)

set.seed(042025)

num_X <- 10
num_units <- 1000
num_spp <- 100
num_archetype <- 5
H <- outer(1:num_X, 1:num_X, "-")
H <- 0.5^abs(H)
covariate_dat <- rmvnorm(num_units, sigma = H) \%>\% 
    as.data.frame \%>\% 
    rename_with(., .fn = function(x) paste0("covariate", x))
rm(H)

true_betas <- runif(num_archetype * num_X, -1, 1) \%>\% matrix(nrow = num_archetype)
true_betas[which(abs(true_betas) < 0.4)] <- 0 # Making archetypal coefficients sparse
true_spp_effects <- matrix(runif(num_spp, -2, 0), ncol = 1)
true_dispparam <- 1/runif(num_spp, 1, 5) 
true_powerparam <- runif(num_spp, 1.4, 1.8)
true_mixprop <- c(0.2, 0.2, 0.3, 0.15, 0.15)
true_betas
 
simdat <- create_samlife(family = nbinom2(), 
formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) \%>\% as.formula, 
data = covariate_dat, 
betas = true_betas, 
spp_effects = true_spp_effects, 
spp_dispparam = true_dispparam, 
spp_powerparam = true_powerparam, 
mixing_proportion = true_mixprop,
seed = 042025)


##----------------------
# First construct regularization path for the mixing proportions in the asSAMs to decide 
# the number of archetypes. Then given this, construct regularization 
# path for the archetypal regression coefficients.
#' **Most users should start here**
##----------------------
## Construct the initial stacked species model fits that will be used throughout 
# the model selection process below
samfit_prefit <- assam(y = simdat$y,
formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) \%>\% as.formula,
data = covariate_dat,
family = nbinom2(),
num_archetypes = 2, #' This is arbitrary
num_cores = detectCores() - 2,
do_assam_fit = FALSE)


## Model selection for the number of archetypes
samfit_select <- passam(y = simdat$y,
formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) \%>\% as.formula,
data = covariate_dat,
family = nbinom2(),
num_archetypes = 10, #' Maximum number of archetypes is ten
selection_on = "mixing_proportions",
supply_quadapprox = samfit_prefit,
num_cores = detectCores() - 2,
beta_selection_control = list(min_df = 5))

samfit_select
samfit_select$regularization_frame
use_BIC <- which.min(samfit_select$BIC)
samfit_select$mixing_proportions_path[,use_BIC]
select_num_archetypes <- sum(samfit_select$mixing_proportions_path[,use_BIC] > 0)
#' The above suggests five archetypes should be included in the asSAM. 


## Now perform selection on the archetypal regression coefficients
# Minimum tuning parameter is such that there are at least five non-zero coefficients
# Note this can take a bit of time...apologies!
samfit_select <- passam(y = simdat$y,
formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) \%>\% as.formula,
data = covariate_dat,
family = nbinom2(),
num_archetypes = select_num_archetypes,
selection_on = "betas",
supply_quadapprox = samfit_prefit,
num_cores = detectCores() - 2,
beta_selection_control = list(min_df = 5))


samfit_select
samfit_select$regularization_frame
use_BIC <- which.min(samfit_select$BIC)
samfit_select$betas_path[,,use_BIC]


## Now fit the final asSAM given a chosen value of the tuning parameter
samfit <- assam(y = simdat$y,
formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) \%>\% as.formula,
data = covariate_dat,
family = nbinom2(),
beta_selection = TRUE,
num_archetypes = num_archetype,
uncertainty_quantification = TRUE,
supply_quadapprox = samfit_prefit,
beta_selection_control = list(lambda = samfit_select$lambda[use_BIC]),
bootstrap_control = list(method = "fast"),
num_cores = detectCores() - 2)

## Alternatively, the final assam can be fitted using "warm starts", which is often
# useful to consistency between the estimated beta and what is obtaind as passam 
samfit2 <- assam(y = simdat$y,
formula = paste("~ ", paste0(colnames(covariate_dat), collapse = "+")) \%>\% as.formula,
data = covariate_dat,
family = nbinom2(),
beta_selection = TRUE,
num_archetypes = num_archetype,
uncertainty_quantification = TRUE,
supply_quadapprox = samfit_prefit,
beta_selection_control = list(lambda = samfit_select$lambda[use_BIC],
warm_start = list(betas = samfit_select$betas_path[,,use_BIC-1],
spp_effects = samfit_select$spp_effects_path[,,use_BIC-1],
spp_nuisance = samfit_select$spp_nuisance_path[,,use_BIC-1],
mixing_proportions = samfit_select$mixing_proportions_path[,use_BIC-1],
posterior_probability = samfit_select$posterior_probability_path[,,use_BIC-1])),
bootstrap_control = list(method = "fast"),
num_cores = detectCores() - 2)


samfit
samfit$betas
true_betas

plot(true_spp_effects, samfit$spp_effects); abline(0,1)
plot(true_dispparam, samfit$spp_nuisance$dispersion, log = "xy"); abline(0,1)
# Note estimates for the archetypal responses and mixture proportions from (as)SAMs should be 
# close to the corresponding true values, *up to a reordering* of the mixture component
# s/archetypes (since the order is essentially arbitrary)
rbind(true_betas, samfit$betas) \%>\% 
t \%>\% 
as.data.frame \%>\%
GGally::ggpairs(.)
table(simdat$archetype_label, apply(samfit$posterior_probability, 1, which.max))


## Demonstrating basic use of functions for asSAM 
summary(samfit)

fitted(samfit)
 
simulate(samfit, data = covariate_dat)

residuals(samfit, type = "dunnsmyth")
 
# Basic residual analysis
plot(samfit, transform_fitted_values = TRUE, envelope = FALSE)
 
# Archetype-level predictions
predict(samfit, newdata = covariate_dat, type = "archetype", se_fit = TRUE) 

# Species-level predictions
predict(samfit, newdata = covariate_dat, type = "species_max", num_cores = 8, se_fit = FALSE) 
}
}
\author{
Francis K.C. Hui \href{mailto:fhui28@gmail.com}{fhui28@gmail.com}
}
