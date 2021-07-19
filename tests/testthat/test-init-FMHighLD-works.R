test_that("init_causal_candidates_coef works without intercept", {
  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  set.seed(123)
  init1 <- init_causal_candidates_coef(data, TRUE, response ~ 0 + annot,
    c("annot" =.5), 1,
    FMParam())
  set.seed(123)
  init2 <- init_causal_candidates_coef(data, TRUE, response ~ 0 + annot,
    c("(Intercept)" = 0, "annot" = .5), 1, FMParam())
  expect_equal(init1, init2)
})



test_that("init_causal_candidates_coef works with intercept", {
  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  set.seed(123)
  init1 <- init_causal_candidates_coef(data, TRUE, response ~ 1 + annot,
    c("(Intercept)" = 1, "annot" = .5), 1,
    FMParam())
  expect_equal(init1[[1]],
    c(ld1 = "snp5", ld2 = "snp8", ld3 = "snp6", ld4 = "snp2"))
})



test_that("init_causal_candidates_coef error when formula has no intcept", {

  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  set.seed(123)
  expect_error(init_causal_candidates_coef(data, TRUE, response ~ 0 + annot,
    c("(Intercept)" = 1, "annot" = .5), 1, FMParam()))

})

test_that("warning with unnamed coefficients", {

  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  set.seed(123)
  expect_warning(init_causal_candidates_coef(data, TRUE, response ~ 1 + annot,
    c(1, 3), 1, FMParam()))

})

test_that("init_causal_candidats_coef error when the dimensions differs", {

  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  set.seed(123)
  expect_error(init_causal_candidates_coef(data, TRUE, response ~ 1 + annot,
    c(1, .5, 4), 1, FMParam()))

})



test_that("init_causal_candidates_random works", {

  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  set.seed(123)
  init1 <- init_causal_candidates_random(data, TRUE, 1)[[1]]
  set.seed(123)
  init2 <- data %>%
    dplyr::group_by(ld_cluster) %>%
    dplyr::sample_n(1)
  init2 <- rlang::set_names(init2$snp, init2$ld_cluster)

  expect_equal(init1, init2)

})

test_that("init_iteration gets and FMIter object", {

  set.seed(1234)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)
  data <- build_fm_tibble(z, annot, ld_cluster, TRUE)

  expect_true(is(init_iteration(response ~ 1 + annot, NULL, data, TRUE, 1,
    FMParam()),
    "FMIter"))
  expect_true(is(init_iteration(response ~ 0 + annot,
    c("annot" = 3), data, TRUE, 1, FMParam()),
    "FMIter"))

})