##########################################################
#### R code for adding noise to correlation matrices  ####
#### J. Hardin, S.R. Garcia, D. Golan                 ####
#### last updated 1/18/2013                           ####
##########################################################

#' adds noise to a user specific matrix
#' @param corr_matrix a pre-defined correlation matrix
#' @param epsilon maximum entry wise random noise
#' @param eidim dimension of the noise space, the smaller the value gets more
#'   noise
#' @return a simulated correlation matrix
#' @export
#' @examples
#' noise_cor(diag(15), epsilon = .5, eidim=2)
#' @note Based on Hardin et al, "A method for generating realistic
#' correlation matrices", Annals of Applied Statistics (2012)
#' @export
#' @importFrom Matrix crossprod
noise_cor <- function(corr_matrix, epsilon = .01, eidim = 2) {

  ndim <- nrow(corr_matrix)
  diag(corr_matrix) <- 1 - epsilon

  ### adding noise to the correlation matrix
  eivect <- replicate(ndim, {
    ei <- runif(eidim, -1, 1)
    sqrt(epsilon) * ei / sqrt(sum(ei^2))
  })

  error <- Matrix::crossprod(eivect)
  corr_matrix + error

}



#' Simulates a block correlation matrix with constant correlation per block
#' where the size and base correlation for each block is user specified
#' @param block_sizes A vector with the dimensions of each block
#' @param corrs A vector of same dimension than `block_sizes` with the base
#'  correlation for each block
#' @param delta An additional parameter for the off diagonal correlations
#' @param epsilon maximum entry wise random noise
#' @param eidim dimension of the noise space, the smaller the value gets more
#'   noise
#' @return a simulated correlation matrix
#' @export
#' @examples
#' rho_block <- c(.9, .2)
#' eps_block <- 0.99 - max(rho_block)
#' sim_block_cor(
#'  block_sizes = c(10, 5), corrs = rho_block, delta = 0.39,
#'  epsilon = eps_block, eidim = 2)
#' @note Based on Hardin et al, "A method for generating realistic
#' correlation matrices", Annals of Applied Statistics (2012)
#' @export
sim_block_cor <- function(block_sizes = c(10, 5, 8, 2, 15, 50),
  corrs = c(0.7, 0.7, 0.5, 0.9, 0.85, 0.4),
  delta = 0.39, epsilon = 0.99 - max(corrs), eidim = 2) {
  stopifnot(length(block_sizes) == length(corrs))
  k <- length(block_sizes)
  stopifnot(k > 0)
  ndim <- sum(block_sizes)
  big_cor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in seq_len(k)) {
    cor <- matrix(
      rep(corrs[i], block_sizes[i] * block_sizes[i]), ncol = block_sizes[i])
    if (i == 1) {
      big_cor[1:block_sizes[1], 1:block_sizes[1]] <- cor
    } else {
      big_cor[(sum(block_sizes[1:(i - 1)]) + 1):sum(block_sizes[1:i]),
        (sum(block_sizes[1:(i - 1)]) + 1):sum(block_sizes[1:i])] <- cor
    }
  }
  diag(big_cor) <- 1 - epsilon
  noise_cor(big_cor, epsilon, eidim)
}


#' Calculates an AR(1) Toeplitz matrix with a block diagonal structure. The size
#' and base correlation for each block is user specified
#' @param block_sizes A vector with the dimensions of each block
#' @param corrs A vector of same dimension than `block_sizes` with the base
#'  correlation for each block
#' @param epsilon maximum entry wise random noise
#' @param eidim dimension of the noise space, the smaller the value gets more
#'   noise
#' @return a simulated correlation matrix
#' @importFrom stats toeplitz
#' @examples
#' rho_top <- c(.7, .9)
#' eps_top = (1-max(rho_top))/(1+max(rho_top)) - .01
#' sim_toeplitz_cor(block_sizes=c(10, 5), corrs = rho_top,
#'  epsilon = eps_top, eidim = 2)
#' @note Based on Hardin et al, "A method for generating realistic
#' correlation matrices", Annals of Applied Statistics (2012)
#' @export
sim_toeplitz_cor <- function(block_sizes = c(10, 5, 8, 2, 15, 50),
  corrs = c(.7, .7, .5, .9, .85, .4), epsilon = .01, eidim = 2) {

  stopifnot(length(block_sizes) == length(corrs))
  k <- length(block_sizes)
  stopifnot(k > 0)
  ndim <- sum(block_sizes)
  big_cor <- matrix(rep(0, ndim * ndim), ncol = ndim)

  for (i in 1:k) {

    top <- corrs[i] ^ (seq_len(block_sizes[i]) - 1)
    cor <- stats::toeplitz(top)

    if (i == 1) {
      big_cor[1:block_sizes[1], 1:block_sizes[1]] <- cor
    } else {
      big_cor[(sum(block_sizes[1:(i - 1)]) + 1):sum(block_sizes[1:i]),
        (sum(block_sizes[1:(i - 1)]) + 1):sum(block_sizes[1:i])] <- cor
    }
  }

  diag(big_cor) <- 1 - epsilon
  noise_cor(big_cor, epsilon, eidim)

}

	
#' Auxiliar function to fill in the rest of the structure of the Hub
#' correlation matrix
#' @param r_max is the maximum user specified correlation
#' @param r_min is the minimum user specified correlation
#' @param power is the power at which the correlations descend
#' @param p is the size of the correlation block
rho_hub_aux <- function(r_max, r_min, power, p) {
  stopifnot(p > 2)
  rho_vec <- sapply(seq_len(p - 1) + 1,
    function(i) r_max - ((i - 2) / (p - 2)) ^ (power * (r_max - r_min)))
  c(1, rho_vec)
}



#' Calculates a Toeplitz matrix with values descending from a user specified
#' range of correlations per block. The matrix has a block diagonal structure
#' and the size and base correlation for each block is user specified.
#'
#' @param block_sizes A vector with the dimensions of each block
#' @param max_corrs A vector of same dimension than `block_sizes` with the max
#'  correlation for each block
#' @param min_corrs A vector of same dimension than `block_sizes` with the min
#'  correlation for each block
#' @param epsilon maximum entry wise random noise
#' @param eidim dimension of the noise space, the smaller the value gets more
#'   noise
#' @return a simulated correlation matrix
#' @importFrom stats toeplitz
#' @examples
#' rho_hub = c(.9, .7)
#' tau_hub = c( (.9 - .7) / (10-2), (.7-.6)/(5-2))
#' eps_hub = min(1-(rho_hub) - 0.75*(tau_hub) ) - .01
#' sim_hub_cor(block_sizes = c(10, 5),
#' max_corrs = rho_hub, min_corrs = c(.7,.6),
#'  power = 1, epsilon = eps_hub, eidim=2)
#' @note Based on Hardin et al, "A method for generating realistic
#' correlation matrices", Annals of Applied Statistics (2012)
#' @export
sim_hub_cor <- function(
  block_sizes = c(10, 5, 8, 7, 15, 50),
  max_corrs = c(.9, .7, .7, .5, .9, .3),
  min_corrs = c(.7, .7, .2, .3, .85, .2),
  power = 1, epsilon = .08, eidim = 2) {

  stopifnot(
    length(block_sizes) == length(max_corrs),
    length(max_corrs) == length(min_corrs))
  k <- length(block_sizes)
  stopifnot(k > 0, all(block_sizes > 2))
  stopifnot(all(min_corrs <= max_corrs))
  ndim <- sum(block_sizes)
  big_cor <- matrix(rep(0, ndim * ndim), ncol = ndim)

  for (i in seq_len(k)) {
    cor <- stats::toeplitz(
      rho_hub_aux(max_corrs[i], min_corrs[i], power, block_sizes[i]))

    if (i == 1) {
      big_cor[1:block_sizes[1], 1:block_sizes[1]] <- cor
    } else {
      big_cor[(sum(block_sizes[1:(i - 1)]) + 1):sum(block_sizes[1:i]),
        (sum(block_sizes[1:(i - 1)]) + 1):sum(block_sizes[1:i])] <- cor
    }
  }
  diag(big_cor) <- 1 - epsilon
  noise_cor(big_cor, epsilon, eidim)

}
