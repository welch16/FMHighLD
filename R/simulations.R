#' Simulated the allelic skew from a gaussian mixture with to components and
#' including a peak overlap event
#' @param n number of SNPs
#' @param p probability of being in each component
#' @param q probability of don't overlapping a peak overlap
#' @param mean mean of the gaussian mixtures (one is N(mean, sd) and the other
#'  is N(-mean, sd))
#' @param sd std. deviation of the allelic skew
#' @return the simulated allelic skew for n SNPs
#' @importFrom stats rnorm runif
#' @export
#' @examples
#' sim_skew(n = 20, p = .7, q = NULL, mean = 3, sd = 2)
#' sim_skew(n = 20, p = .7, mean = 3, sd = 2)
#' sim_skew(n = 20, p = .7, q = .8, mean = 3, sd = 2)
sim_skew <- function(n, p = .5, q = NULL, mean = 0, sd = 1) {

  stopifnot(
    is.numeric(n),
    is.numeric(p),
    is.numeric(q) | is.null(q),
    is.numeric(mean),
    is.numeric(sd))
  stopifnot(n > 0, p > 0, p < 1, q > 0, q < 1, sd > 0)

  n <- ceiling(n)

  u <- stats::runif(n)
  z1 <- stats::rnorm(n, mean = -mean, sd = sd)
  z2 <- stats::rnorm(n, mean = mean, sd = sd)

  # add allelic effect to skew
  idx <- u <= p
  z <- z1
  z[!idx] <- z2[!idx]

  # add peak overlap probability
  if (!is.null(q)) {
    z[stats::runif(n) < 1 - q] <- 0
  }
  z

}
