test_that("noise_cor is similar to matrix with low noise", {

  ident_mat <- diag(rep(1, 10))

  expect_true(norm(noise_cor(ident_mat, 1e-6, 2) - ident_mat) <= 1e-5)
})

test_that("noise_cor works", {

  set.seed(123)
  corr_matrix <- diag(rep(1:20))
  epsilon <- .3
  eidim <- 2
  ndim <- nrow(corr_matrix)

  diag(corr_matrix) <- 1 - epsilon

  ### adding noise to the correlation matrix
  eivect <- replicate(ndim, {
    ei <- runif(eidim, -1, 1)
    sqrt(epsilon) * ei / sqrt(sum(ei^2))
  })

  error <- Matrix::crossprod(eivect)
  set.seed(123)
  expect_equal(corr_matrix + error, noise_cor(corr_matrix, epsilon, eidim))

})

test_that("noise_cor error if the input is not matrix", {

  expect_error(noise_cor(1:10))

})

test_that("noise_cor error if epsilon < 0", {

  expect_error(noise_cor(diag(rep(1, 10)), epsilon = -1, 2))
  expect_error(noise_cor(diag(rep(1, 10)), epsilon = 0, 2))

})

test_that("noise_cor error if eidim is not numeric", {

  expect_error(noise_cor(diag(rep(1, 10)), epsilon = 3, "av"))

})

test_that("sim_block_cor error when block_sizes and corrs have diff. sizes", {
  block_sizes <- c(1, 2, 3)
  corrs <- c(.3, .5)
  expect_error(sim_block_cor(block_sizes, corrs))

})

test_that("sim_toeplitz_cor error when block_sizes and corrs have diff sizes", {
  block_sizes <- c(1, 2, 3)
  corrs <- c(.3, .5)
  expect_error(sim_toeplitz_cor(block_sizes, corrs))

})

test_that("sim_block_cor works", {

  block_sizes <- c(1, 2, 3)
  corrs <- c(.3, .5, .7)
  mat <- sim_block_cor(block_sizes, corrs)
  expect_true(ncol(mat) == nrow(mat))
  expect_true(nrow(mat) == sum(block_sizes))
  expect_true(is(mat, "matrix"))

})

test_that("sim_toeplitz_cor works", {

  block_sizes <- c(4, 2, 3)
  corrs <- c(.6, .5, .7)
  mat <- sim_toeplitz_cor(block_sizes, corrs)
  expect_true(ncol(mat) == nrow(mat))
  expect_true(nrow(mat) == sum(block_sizes))
  expect_true(is(mat, "matrix"))

})

test_that("sim_hub_cor works", {

  rho_hub <- c(.9, .7)
  tau_hub <- c((.9 - .7) / (10 - 2), (.7 - .6) / (5-2))
  eps_hub <- min(1 - (rho_hub) - 0.75 * (tau_hub)) - .01
  block_sizes <- c(10, 5)
  mat <- sim_hub_cor(block_sizes = block_sizes,
    max_corrs = rho_hub, min_corrs = c(.7, .6),
    power = 1, epsilon = eps_hub, eidim = 2)
  expect_true(ncol(mat) == nrow(mat))
  expect_true(nrow(mat) == sum(block_sizes))
  expect_true(is(mat, "matrix"))

})



test_that("get_ld_clusters works", {

  set.seed(12345)
  nsnps <- 20
  nldblocks <- 4

  nsnps / nldblocks

  snp_names <- stringr::str_c("snp", seq_len(nsnps))

  ld_mat <- sim_block_cor(
    block_sizes = rep(nsnps / nldblocks, nldblocks),
    corrs = rep(.8, nldblocks),
    delta = .2,
    epsilon = .1,
    eidim = 2)

  rownames(ld_mat) <- colnames(ld_mat) <- snp_names

  expect_equal(length(table(get_ld_clusters(ld_matrix = ld_mat, 0.7))),
    nldblocks)
  expect_true(all(
    table(get_ld_clusters(ld_matrix = ld_mat, 0.7)) == nsnps / nldblocks))

})

test_that("get_ld_clusters error when min_r2 is not correlation", {

  set.seed(12345)
  nsnps <- 20
  nldblocks <- 4

  nsnps / nldblocks

  snp_names <- stringr::str_c("snp", seq_len(nsnps))

  ld_mat <- sim_block_cor(
    block_sizes = rep(nsnps / nldblocks, nldblocks),
    corrs = rep(.8, nldblocks),
    delta = .2,
    epsilon = .1,
    eidim = 2)

  rownames(ld_mat) <- colnames(ld_mat) <- snp_names

  expect_error(get_ld_clusters(ld_mat, min_r2 = 3))
  expect_error(get_ld_clusters(ld_mat, min_r2 = -1))
})

test_that("get_ld_clusters tidy works", {

  set.seed(12345)
  nsnps <- 20
  nldblocks <- 4

  nsnps / nldblocks

  snp_names <- stringr::str_c("snp", seq_len(nsnps))

  ld_mat <- sim_block_cor(
    block_sizes = rep(nsnps / nldblocks, nldblocks),
    corrs = rep(.8, nldblocks),
    delta = .2,
    epsilon = .1,
    eidim = 2)

  rownames(ld_mat) <- colnames(ld_mat) <- snp_names
  out <- get_ld_clusters(ld_mat, .7, tidy = FALSE)
  expect_true(is(get_ld_clusters(ld_mat, .7, tidy = TRUE), "data.frame"))
  expect_equal(get_ld_clusters(ld_mat, .7, tidy = TRUE),
    tibble::tibble(snp = names(out), ld_cluster = out))

})

test_that("get_ld_clusters works without rownames or colnames", {

  set.seed(12345)
  nsnps <- 20
  nldblocks <- 4

  snp_names <- stringr::str_c("snp", seq_len(nsnps))

  ld_mat <- sim_block_cor(
    block_sizes = rep(nsnps / nldblocks, nldblocks),
    corrs = rep(.8, nldblocks),
    delta = .2,
    epsilon = .1,
    eidim = 2)

  rownames(ld_mat) <- snp_names
  expect_equal(length(table(get_ld_clusters(ld_matrix = ld_mat, 0.7))),
    nldblocks)
  expect_true(all(
    table(get_ld_clusters(ld_matrix = ld_mat, 0.7)) == nsnps / nldblocks))

  ld_mat <- sim_block_cor(
    block_sizes = rep(nsnps / nldblocks, nldblocks),
    corrs = rep(.8, nldblocks),
    delta = .2,
    epsilon = .1,
    eidim = 2)

  colnames(ld_mat) <- snp_names
  expect_equal(length(table(get_ld_clusters(ld_matrix = ld_mat, 0.7))),
    nldblocks)
  expect_true(all(
    table(get_ld_clusters(ld_matrix = ld_mat, 0.7)) == nsnps / nldblocks))



})

