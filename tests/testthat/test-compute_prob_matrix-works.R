test_that("multiplication works", {

  nsnps <- 10
  probmat <- cbind(rep(1, 10), rep(0, 10), rep(0, 10))

  expect_equal(compute_mixture_prob(probmat), c(1, 0, 0))

  set.seed(1234)
  probmat <- matrix(runif(2 * nsnps), ncol = 2)
  probs <- colMeans(probmat)
  expect_equal(compute_mixture_prob(probmat), probs / sum(probs))

})
