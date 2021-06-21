test_that("build_design_matrix works", {

  set.seed(12345)
  nsnps <- 2000
  nldblocks <- 20
  annot_beta <- c("(Intercept)" = 0, "skew" =  3)

  snp_data <- simulate_FMHighLD_simple(nsnps, nldblocks, 0.2, 0.1, 2,
    .5, .9, 4, 1^2, annot_beta, .5^2)

  snp_names <- snp_data$snp

  skew_mat <- build_design_matrix(snp_data, annot_beta)
  rownames(skew_mat) <- snp_names

  mat <- Matrix::Matrix(cbind(1, snp_data$skew))
  rownames(mat) <- snp_names
  colnames(mat) <- colnames(skew_mat)

  expect_equal(skew_mat, mat)

})
