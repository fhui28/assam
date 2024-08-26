#' @title Basic residual diagnostic plots for an approximate and scalable SAM object
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Five potential plots are currently available for some basic residual diagnostics for a fitted \code{assam} model: 1) a plot of residuals against (potentially transformed) fitted values; 2) a normal probability plot of residuals with simulated point-wise 95\% confidence interval envelope; 3) plot of residuals against observational unit index; 4) a plot of residuals again column index; 5) scale-location plot using (potentially transformed) fitted values.
#' 
#' @param x An object of class \code{assam}.
#' @param y The multivariate abundance data used as part of the \code{assam} object. This needs to be supplied since currently \code{assam} objects do not save the responses to save memory.
#' @param which_plot If a subset of the plots is desired, then a vector containing subset of the integers 1, 2, 3, 4, 5.
#' @param type The type of residuals to be used in constructing the plots. Currently the options available are: "PIT" and "dunnsmyth".
#' @param transform_fitted_values For plots 1 and 5, should the fitted values be transformed using the link function specified in the asSAM i.e.. what is available in \code{x$family$linkfun}. This can be useful visually so that observations are not too bunched up on the x-axis of plots 1 and 5.  
#' @param titles Titles to appear above each plot.
#' @param species_colors Either a scalar is supplied or a vector with length of number of species in the spatio-temporal multivariate abundance data. If the former than all species use this color in the plots. If the latter then the vector specified the colors to use for each species. Defaults to \code{NULL}, which results in each species having a unique color based on the [grDevices::rainbow()] palette.
#' @param smooth Should a smoother be added to each plot?
#' @param envelope Should approximate simulation envelopes be constructed for the normal probability plot? Default to \code{TRUE}. Note if \code{envelope = FALSE} then \code{envelope_col} and \code{envelope_K} are ignored.
#' @param envelope_col A vector of length 2, specifying the colors to use for the lines and shade respectively for the approximate simulation envelopes.
#' @param envelope_rep The number of simulations to use in order to build the approximate simulation envelope.
#' @param which_species A vector indexing the species to plot, if the residual plots should be constructed for only a subset of species. Defaults to \code{NULL}, in which case all species are plotted. This may be useful if the number of species is quite large.
#' @param seed This can be used set the seed when constructing the PIT and Dunn-Smyth residuals, which for discrete responses involve some degree of jittering.  
#' @param ... Additional graphical arguments.
#' 
#' @details 
#' This function is heavily adapted from [CBFM::plot.CBFM()] and [gllvm::plot.gllvm()]. As basic visual diagnostics, these plots should behave as follows: 
#' 1. the plot of residuals versus fitted values should not exhibit any noticeable pattern e.g., no (inverse) fan-shape or a trend; 
#' 2. the normal probability plot should have the residuals lying approximately on a straight line and almost all residuals lying within the approximate simulation envelopes; 
#' 3. a plot of residuals against observational unit index should not exhibit any noticeable pattern. The plot can also be used to look for potential outlying observational units; 
#' 4. a plot of residuals against species index should not exhibit any noticeable pattern. The plot can also be used to look for potential outlying species; 
#' 5. the scale-location plot should not exhibit any trend. 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>
#' 
#' @examples
#' \dontrun{
#' #' Please see the help file for assam for example.
#' }
#' 
#' @export
#' @importFrom graphics abline boxplot lines panel.smooth par plot points polygon
#' @importFrom grDevices rainbow
#' @importFrom mgcv gam predict.gam
#' @importFrom stats qnorm qqnorm qqline quantile rnorm
#' @importFrom tdigest tdigest tquantile
#'
#' @md

plot.assam <- function(x, 
                       y,
                       which_plot = 1:5, 
                       type = "dunnsmyth", 
                       transform_fitted_values = FALSE,
                       titles = c("Residuals vs. fitted values", "Normal probability plot", "Residuals vs. unit index", "Residuals vs. species index","Scale-Location plot"), 
                       species_colors = NULL, 
                       smooth = TRUE, 
                       envelope = TRUE, 
                       envelope_col = c("grey","lightgrey"), 
                       envelope_rep = 100, 
                       which_species = NULL, 
                       seed = NULL, 
                       ...) {
    
    num_units <- dim(x$linear_predictor)[1]
    num_spp <- dim(x$linear_predictor)[2]
    
    sppind <- 1:num_spp
    if(!is.null(which_species))
        sppind <- sort(which_species)
    
    if(length(sppind) > num_spp)
        stop("Length of which_species exceeded the number of species.")
    if(any(which_species > num_spp))
        stop("which_species should be a vector of integers ranging from 1 to the number of species.")
    
    if(any(which_plot > 5))
        stop("which_plot should be a vector of integers ranging from 1 to 5. There are only five possible plots that this function currently offers.")
    
    #------------------
    #' # Form plot titles
    #------------------
    mains <- rep("", 5)
    mains[which_plot] <- titles[which_plot]
    
    res <- residuals.assam(object = x, y = y, type = type, seed = seed)
    if(any(abs(res) > 1e3, na.rm = TRUE))
        warning("Some extremely large residuals (> 1000 in absolute value) will be left out of the plotting.")
    res[res < -1e3] <- -1e3
    res[res > 1e3] <- 1e3
    dsres <- res[, sppind]
    etamat <- fitted.assam(x, type = "mean")
    if(transform_fitted_values)
        etamat <- x$family$linkfun(etamat)
    xxx <- boxplot(c(etamat), outline = FALSE, plot = FALSE)$stats     
    yyy <- range(c(dsres[dsres > -1e3 & dsres < 1e3]), na.rm = TRUE)     
    
    
    #------------------
    #' # Form colors for species - done by prevalence
    #------------------
    csum <- order(colSums(as.matrix(y))[sppind])
    if(!is.null(species_colors)) {
        col <- rep(1, num_spp)
        col[1:num_spp] <- species_colors
    } 
    if(is.null(species_colors)) {
        if(num_spp < 8)
            col <- (1:num_spp)[csum]
        else
            col <- rainbow(num_spp + 1)[2:(num_spp + 1)][csum]
    }
    
    
    gr.pars <- list(...)
    par(...)
    
    
    #------------------
    #' # Residuals versus fitted values
    #------------------
    if(1 %in% which_plot) {
        if(is.null(gr.pars$xlim)) {
            plot(etamat, 
                 dsres, 
                 xlab = "Fitted values (potentially transformed)", 
                 ylab = "Residuals", 
                 type = "n", 
                 col = rep(col, each = num_units), 
                 main = mains[1], 
                 xlim = c(min(xxx), max(xxx)), 
                 ylim = yyy)
            abline(0, 0, col = "grey", lty = 3)
            } 
        else {
            plot(etamat, 
                 dsres, 
                 xlab = "Fitted values (potentially transformed)", 
                 ylab = "Residuals", 
                 type = "n", 
                 col = rep(col, each = num_units),
                 main = mains[1], 
                 ...)
            abline(0, 0, col = "grey", lty = 3)
            }
        
        if(smooth) 
            .gamEnvelope(etamat, dsres, col = rep(col, each = num_units), envelopes = TRUE, envelope.col = envelope_col, ...)
        
        }
    
    
    #------------------
    #' # Normal probability or quantile-quantile plot of residuals with an approximate point-wise 95\% confidence interval envelope          
    #------------------
    if(2 %in% which_plot) {
        qq.x <- qqnorm(c(dsres), 
                       main = mains[2], 
                       ylab = "Residuals", 
                       col = rep(col, each = num_units), 
                       cex = 0.5, 
                       xlab = "Theoretical quantiles", 
                       ylim = yyy, 
                       type = "n")
        
        num_obs <- num_units * num_spp
        if(envelope) {
            message("Constructing (approximate) simulation envelopes for normal probability plot...")
            
            yy <- quantile(dsres, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
            xx <- qnorm(c(0.25, 0.75))
            slope <- diff(yy) / diff(xx)
            int <- yy[1] - slope * xx[1]
            all_ris <- matrix(rnorm(sum(!is.na(qq.x$x)) * envelope_rep, mean = int, sd = slope), ncol = envelope_rep)
            Ym <- apply(all_ris, 2, sort)
            rm(all_ris)
            cis <- apply(Ym, 1, function(x) { 
                out <- try(tquantile(tdigest(x, 1000), probs = c(0.025, 0.975)), silent = TRUE)
                if(inherits(out, "try-error"))
                    out <- quantile(x, probs = c(0.025, 0.975))
                return(out)
            })
            rm(Ym)
            Xm <- sort(qq.x$x)
            
            polygon(Xm[c(1:length(Xm),length(Xm):1)], 
                    c(cis[1,],cis[2, length(Xm):1]), 
                    col = envelope_col[2], 
                    border = NA)
            }
        
        points(qq.x$x, 
               qq.x$y, 
               col = rep(col, each = num_units), 
               cex = 0.5)
        qqline(c(dsres), 
               col = envelope_col[1])
        }
    
    
    #------------------
    # Residuals against observational unit index          
    #------------------
    if(3 %in% which_plot) {
        plot(rep(1:num_units, num_spp), 
             dsres, 
             xlab = "Unit index", 
             ylab = "Residuals", 
             col = rep(col, each = num_units), 
             main = mains[3], 
             ylim = yyy,
             ...);
        abline(0, 0, col = "grey", lty = 3)
        if(smooth) 
            panel.smooth(rep(1:num_units, num_spp), 
                         dsres, 
                         col = rep(col, each = num_units), 
                         col.smooth = envelope_col[1], 
                         ...)
        }
    
    
    #------------------
    # Residuals against species index          
    #------------------
    if(4 %in% which_plot) {
        plot(rep(1:num_spp, each = num_units), 
             dsres, 
             xlab = "Species index", 
             ylab = " Residuals", 
             col = rep(col[csum], each = num_units), 
             main = mains[4], 
             ylim = yyy, 
             ...) 
        abline(0, 0, col = "grey", lty = 3)
        if(smooth) 
            panel.smooth(rep(1:num_spp, each = num_units), 
                         dsres, 
                         col = rep(col[csum], each = num_units), 
                         col.smooth = envelope_col[1], 
                         ...)
        }
    
    
    #------------------
    # Scale-location plot
    #------------------
    if(5 %in% which_plot) {
        sqres <- sqrt(abs(dsres))
        yyy <- range(sqres[dsres > -1e3 & dsres < 1e3], na.rm = TRUE)
        yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Residuals"))))
        
        if(is.null(gr.pars$xlim)) {
            plot(etamat, 
                 sqres, 
                 xlab = "Fitted values (potentially transformed)", 
                 ylab = yl, 
                 col = rep(col, each = num_units), 
                 main = mains[5], 
                 xlim = c(min(xxx), max(xxx)), 
                 ylim = yyy, 
                 ...)
            } 
        else {
            plot(etamat, 
                 sqres, 
                 xlab = "Fitted values", 
                 ylab = yl, 
                 col = rep(col, each = num_units), 
                 main = mains[5], 
                 ...)
            }
        
        if(smooth) 
            panel.smooth(etamat, 
                         sqres, 
                         col = rep(col, each = num_units), 
                         col.smooth = envelope_col[1], 
                         ...)
        }
    
}


## Modified from gllvm package. Thanks to Jenni Niku in the gllvm package for this function!
.gamEnvelope <- function(x, 
                         y, 
                         line.col = "black", 
                         envelope.col = c("grey","lightgrey"), 
                         col = 1, 
                         envelopes = TRUE, 
                         subsample = 5000, ...) {
    xSort <- sort(x, index.return = TRUE)
    gam.yx <- mgcv::gam(resp ~ cov, 
                        data = data.frame(resp = y[xSort$ix], cov = xSort$x))
    pr.y <- mgcv::predict.gam(gam.yx, 
                        se.fit = TRUE, 
                        newdata = data.frame(cov = xSort$x))
    
    prHi <- pr.y$fit + 1.96*pr.y$se.fit
    prLow <- pr.y$fit - 1.96*pr.y$se.fit
    n.obs <- length(prLow)     
    sel_x_index <- 1:n.obs
    
    if(envelopes) 
        polygon(xSort$x[c(sel_x_index,rev(sel_x_index))], 
                c(prHi,prLow[rev(sel_x_index)]), 
                col = envelope.col[2], 
                border = NA)
    
    lines(xSort$x, pr.y$fit, col = envelope.col[1])
    abline(h = 0, col = 1)
    points(x, y, col = col, ...)
}


