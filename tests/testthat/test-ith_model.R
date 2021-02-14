set.seed(123)
nsnps <- 500
snp_data <- data.frame(x = rnorm(nsnps))
snp_data$y <- snp_data$x
snp_data$gg <- factor(if_else(seq_len(nsnps) <= nsnps / 2, "g1", "g2"))
gamma <- matrix(rep(1, 2 * nsnps), ncol = 2)

test_that("ith model works single trait", {


  model <- lm(y ~ 0 + x, data = snp_data)
  coef(model)

  expect_equivalent(
    coef(compute_ith_model(2, y ~ 0 + x, snp_data, gamma, TRUE)),
    coef(model))

})

test_that("ith model works multi trait", {

  model <- lme4::lmer(y ~ 0 + x + x | gg, data = snp_data)
  expect_equivalent(
    coef(compute_ith_model(2, y ~ 0 + x | gg, snp_data, gamma, FALSE)),
    coef(model))


})

test_that("ith model error, multitrait without random effect", {

  expect_error(
    compute_ith_model(2, y ~ 0 + x, snp_data, gamma, FALSE))

})
