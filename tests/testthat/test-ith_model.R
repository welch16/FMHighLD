set.seed(123)
nsnps <- 500
snp_data <- data.frame(snp = as.character(seq_len(nsnps)), x = rnorm(nsnps))
snp_data$y <- snp_data$x
snp_data$gg <- factor(if_else(seq_len(nsnps) <= nsnps / 2, "g1", "g2"))
gamma <- matrix(rep(1, 2 * nsnps), ncol = 2)

test_that("ith model works single trait", {

  model <- lm(y ~ 0 + x, data = snp_data)
  coef(model)

  expect_equivalent(
    coef(compute_ith_model(2, as.character(seq_len(nsnps)),
      y ~ 0 + x, snp_data, gamma, TRUE)),
    coef(model))

})

test_that("ith model works single trait with causal candidates", {

  causal_candidates <- sample(snp_data$snp, nsnps / 2)
  model <- lm(y ~ 0 + x,
    data = dplyr::filter(snp_data, snp %in% causal_candidates))
  coef(model)
  gamma <- matrix(rep(1, 2 * length(causal_candidates)), ncol = 2)

  expect_equivalent(
    coef(compute_ith_model(2, causal_candidates,
      y ~ 0 + x, snp_data, gamma, TRUE)),
    coef(model))


})

test_that("ith model works multi trait", {

  suppressWarnings({
    model <- lme4::lmer(y ~ 0 + x + x | gg, data = snp_data)
  })

  suppressWarnings({
    fmld_model <- compute_ith_model(2, as.character(seq_len(nsnps)),
      y ~ 0 + x | gg, snp_data, gamma, FALSE)
  })
  expect_equivalent(
    coef(fmld_model),
    coef(model))


})

test_that("ith model error, multitrait without random effect", {

  expect_error(
    compute_ith_model(2, as.character(seq_len(nsnps)), y ~ 0 + x, snp_data,
      gamma, FALSE))

})
