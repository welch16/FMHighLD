#' Compute the Kullback-Leibler (KL) divergence between two discrete
#' probability distributions
#' @param p1 the first probability vector
#' @param p2 the second probability vector
#' @return a scalar value with the KL-divergence between both
#'  distributions
#' @export
#' @examples
#'
#' u1 <- runif(100)
#' u1 <- u1 / sum(u1)
#' v1 <- runif(100)
#' v1 <- v1 / sum(v1)
#' kl(u1, v1)
kl <- function(p1, p2) {
  sum(p1 * (log(p1) - log(p2)))
}

#' Computes the Jensen-Shannon (JS) divergence between two discrete
#' probability distributions
#' @param p1 the first probability vector
#' @param p2 the second probability vector
#' @return a scalar value with the JS-divergence between both
#'  distributions
#' @export
#' @examples
#'
#' u1 <- runif(100)
#' u1 <- u1 / sum(u1)
#' v1 <- runif(100)
#' v1 <- v1 / sum(v1)
#' jsd(u1, v1)
jsd <- function(p1, p2) {
    m <- .5 * (p1 + p2)
    .5 * kl(p1, m) + .5 * kl(p2, m)
}
