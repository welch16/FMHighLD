
test_that("causal rule works with default FMParam() rule", {

  set.seed(12345)
  nsnps <- 10
  nldblocks <- 2
  annot_beta <- c("(Intercept)" = 0, "skew" =  3)

  snp_data <- simulate_FMHighLD_simple(nsnps, nldblocks, 0.2, 0.1, 2,
    .5, .9, 4, 1^2, annot_beta, .5^2)

  res <- rep(1, nsnps)
  res[1] <- .1
  res[6] <- .1
  rule <- causal_rule(snp_data, residuals = res,
    FMParam(), group = "ld_cluster", annot_names = "skew")

  expect_equal(rule$which_snp, stringr::str_c("snp", c(1, 6)))

})

test_that("causal rule work  with 'all' strategy in FMParam()", {


  set.seed(12345)
  nsnps <- 30
  nldblocks <- 3
  annot_beta <- c("(Intercept)" = 0, "skew" =  3)

  snp_data <- simulate_FMHighLD_simple(nsnps, nldblocks, 0.2, 0.1, 2,
    .5, .9, 4, 1^2, annot_beta, .5^2)

  res <- rep(1, nsnps)
  res[1] <- .1
  res[11] <- .1
  res[21] <- .1
  rule <- causal_rule(snp_data, residuals = res,
    FMParam(strategy = "all", prob = .05),
      group = "ld_cluster", annot_names = "skew")

  expect_equal(rule$which_snp, stringr::str_c("snp", c(1, 11, 21)))

})
