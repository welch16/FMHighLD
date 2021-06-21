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

#' Simulates a simple dataset with one causal per LD block
#' @param nsnps Number of SNPs to simulate
#' @param nldblocks Number of LD blocks to simulate
#' @param ld_delta An additional parameter for the off diagonal correlations
#' @param ld_epsilon maximum entry wise random noise
#' @param ld_eidim dimension of the noise space, the smaller the value gets more
#'   noise
#' @param skew_p probability of being in each component
#' @param skew_q probability of don't overlapping a peak overlap
#' @param skew_mean mean of the gaussian mixtures (one is N(mean, sd) and the
#'  other is N(-mean, sd))
#' @param skew_sigma2 variance of the allelic skew
#' @param z_beta regression coefficients of the z-statistic as function of the
#'  annotation data
#' @param z_sigma2 variance of the z-statistic as function of the annotation
#'  data
#' @export
#' @importFrom stringr str_c
#' @importFrom purrr map map_chr
#' @importFrom dplyr mutate if_else
#' @importFrom tibble tibble
#' @examples
#' simulate_FMHighLD_simple(2000, 200, .2, .1, 2, .5, .9, 4, 1, c(0, 3), .5^2)
simulate_FMHighLD_simple <- function(nsnps, nldblocks, ld_delta, ld_epsilon,
  ld_eidim, skew_p, skew_q, skew_mean, skew_sigma2,
  z_beta, z_sigma2) {

  message("generating data with ", floor(nsnps / nldblocks),
    " variants per LD block")

  snp_names <- stringr::str_c("snp", seq_len(nsnps))

  ld_mat <- sim_block_cor(
    block_sizes = rep(floor(nsnps / nldblocks), nldblocks),
    corrs = rep(.8, nldblocks),
    delta = ld_delta,
    epsilon = ld_epsilon,
    eidim = ld_eidim)

  rownames(ld_mat) <- snp_names
  colnames(ld_mat) <- snp_names

  skew <- sim_skew(nsnps, p = skew_p, q = skew_q,
    mean = skew_mean, sd = sqrt(skew_sigma2))
  skew_mat <- as.matrix(skew, ncol = 1)
  rownames(skew_mat) <- snp_names
  colnames(skew_mat) <- "annot"
  skew_mat <- cbind(1, skew_mat)

  ld_clusters <- get_ld_clusters(ld_matrix = ld_mat, 0.7)

  causals <- split(ld_clusters, ld_clusters)
  causals <- purrr::map(causals, ~ sample(., 1))
  causals <- purrr::map_chr(causals, names)

  snp_data <- tibble::tibble(
    snp = snp_names,
    ld_cluster = ld_clusters)
  snp_data <- dplyr::mutate(snp_data,
      causal = dplyr::if_else(snp %in% causals, 1, 0),
      skew)

  dplyr::mutate(snp_data,
      z = rnorm(nsnps, dplyr::if_else(causal == 1,
        skew_mat %*% z_beta, 0.0), sd = sqrt(z_sigma2)))

}
