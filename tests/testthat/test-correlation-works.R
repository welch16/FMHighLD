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