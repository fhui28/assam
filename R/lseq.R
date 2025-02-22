#' @title Produce a sequence of tuning parameters on the log (base 10) scale.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' @param from Minimum tuning parameter in sequence.
#' @param to Maximum tuning parameter in sequence.
#' @param length Length of tuning parameter sequence.
#' @param decreasing Should the sequence of tuning parameters be increasing or decreasing? Defaults to \code{FALSE} i.e., increasing.
#'
#' @return A sequence of tuning parameters.
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>
#'
#'
#' @export lseq
#' @md

lseq <- function (from, to, length, decreasing = FALSE) {
    stopifnot(from > 0)
    out <- 10^(seq(log10(from), log10(to), length.out = length))
    out <- out[order(out, decreasing = decreasing)]
    return(out)
}
