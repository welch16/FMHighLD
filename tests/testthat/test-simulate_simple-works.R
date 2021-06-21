test_that("simple simulation works", {

  set.seed(12345)
  nsnps <- 2000
  nldblocks <- 200
  delta_ld <- 0.2
  epsilon_ld <- 0.1
  eidim_ld <- 2

  snp_names <- stringr::str_c("snp", seq_len(nsnps))

  ld_mat <- sim_block_cor(
    block_sizes = rep(nsnps / nldblocks, nldblocks),
    corrs = rep(.8, nldblocks),
    delta = delta_ld,
    epsilon = epsilon_ld,
    eidim = eidim_ld)

  rownames(ld_mat) <- colnames(ld_mat) <- snp_names

  skew <- sim_skew(nsnps, p = 0.5, q = .9, mean = 4, sd = 1)
  skew_mat <- as.matrix(skew, ncol = 1)
  rownames(skew_mat) <- snp_names
  colnames(skew_mat) <- "annot"

  ld_clusters <- get_ld_clusters(ld_matrix = ld_mat, 0.7)

  causals <- split(ld_clusters, ld_clusters) %>%
    purrr::map(~ sample(., 1)) %>%
    purrr::map_chr(names)

  snp_data <- tibble::tibble(
    snp = snp_names,
    ld_cluster = ld_clusters) %>%
    dplyr::mutate(
      causal = if_else(snp %in% causals, 1, 0), skew) %>%
    dplyr::mutate(
      z = rnorm(nsnps, if_else(causal == 1, 0 +  3 * skew, 0.0), sd = .5))

  set.seed(12345)
  snp_data_fmld <-
    simulate_FMHighLD_simple(nsnps, nldblocks, delta_ld, epsilon_ld,
      eidim_ld, .5, .9, 4, 1, c(0, 3), .5^2)

  expect_equal(snp_data, snp_data_fmld)
})
