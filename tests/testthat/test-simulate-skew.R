
test_that("expect error due to non-numeric input", {

  nsnps <- 10
  expect_error(sim_skew("a", p = 0.5, q = .9, mean = 4, sd = 1))
  expect_error(sim_skew(nsnps, p = "a", q = .9, mean = 4, sd = 1))
  expect_error(sim_skew(nsnps, p = 0.5, q = "a", mean = 4, sd = 1))
  expect_error(sim_skew(nsnps, p = 0.5, q = .9, mean = "a", sd = 1))
  expect_error(sim_skew(nsnps, p = 0.5, q = .9, mean = 4, sd = "a"))

})

test_that("expect error due to negative number of snps", {

  expect_error(sim_skew(-10, p = 0.5, q = .9, mean = 4, sd = 1))
})

test_that("expect error due to invalid probabilities", {

  expect_error(sim_skew(10, p = -1, q = .9, mean = 4, sd = 1))
  expect_error(sim_skew(10, p = 1.3, q = .9, mean = 4, sd = 1))
  expect_error(sim_skew(10, p = 1, q = .9, mean = 4, sd = 1))
  expect_error(sim_skew(10, p = 0.5, q = -1, mean = 4, sd = 1))
  expect_error(sim_skew(10, p = 0.5, q = 1.2, mean = 4, sd = 1))
  expect_error(sim_skew(10, p = 0.5, q = 1, mean = 4, sd = 1))

})

test_that("expect number of snps", {

  expect_length(sim_skew(10, p = .3, q = .9, mean = 3, sd = 1), 10)

})

test_that("expect same ouput", {

  set.seed(1234)
  n <- 20
  mean <- 3
  sd <- 0.3
  p <- 0.5
  q <- 0.8
  u <- runif(n)
  z1 <- stats::rnorm(n, mean = -mean, sd = sd)
  z2 <- stats::rnorm(n, mean = mean, sd = sd)

  idx <- u <= p
  z <- z1
  z[!idx] <- z2[!idx]
  z[stats::runif(n) < 1 - q] <- 0

  set.seed(1234)
  expect_equivalent(z, sim_skew(n, p, q, mean, sd))

})
