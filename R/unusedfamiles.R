#' @noMd
#' @noRd
betalogitfam <- function() {
    link <- "logit"
    linkfun <- function(mu) 
        return(log(mu/(1-mu)))
    linkinv <- function(eta) 
        return(exp(eta)/(1+exp(eta)))
    mu.eta <- function(eta) {
        ifelse(abs(eta)>30,.Machine$double.eps, exp(eta)/(1+exp(eta))^2) 
    }
    variance <- function(mu, phi) 
        return(mu*(1-mu)/(1+phi))
    
    structure(list(family = "beta", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
}


#' @noMd
#' @noRd
tweedielogfam <- function() {
        link <- "log"
        linkfun <- function(mu) 
            return(log(mu))
        linkinv <- function(eta) 
            return(pmax(exp(eta), .Machine$double.eps))
        mu.eta <- function(eta) 
            return(pmax(exp(eta), .Machine$double.eps))
        variance <- function(mu, power, phi) 
            return(pmax(phi * mu^power, .Machine$double.eps))
        
        structure(list(family = "tweedie", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
    }

#' @noMd
#' @noRd
nb2 <- function() {
        ##--------------------------
        ## Feed the family!
        ##--------------------------
        link <- "log"
        linkfun <- function(mu) 
            return(log(mu))
        linkinv <- function(eta) 
            return(pmax(exp(eta), .Machine$double.eps))
        mu.eta <- function(eta) 
            return(pmax(exp(eta), .Machine$double.eps))
        variance <- function(mu, phi) 
            return(pmax(mu + mu^2/phi, .Machine$double.eps))
        
        structure(list(family = "negative.binomial", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
    }
