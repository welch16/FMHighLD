mat <- Matrix::bdiag(matrix(rep(1, 9), nrow = 3),
  matrix(rep(.8, 16), nrow = 4))
mat <- as.matrix(mat)
out1 <- c(stringr::str_c("ld", rep(1, 3)),
  stringr::str_c("ld", 2:5))
names(out1) <- stringr::str_c("snp", seq_along(out1))
out2 <- c(stringr::str_c("ld", rep(1, 3)),
  stringr::str_c("ld", rep(2, 4)))
names(out2) <- stringr::str_c("snp", seq_along(out1))

test_that("baseline ld cluster works", {

  expect_equal(get_ld_clusters(mat, 1), out1)
  expect_equal(get_ld_clusters(mat, .8), out2)

})

test_that("baseline ld cluster works tidy", {

  tib1 <- tibble::tibble(snp = names(out1), ld_cluster = out1)
  tib2 <- tibble::tibble(snp = names(out2), ld_cluster = out2)

  expect_equal(get_ld_clusters(mat, 1, tidy = TRUE), tib1)
  expect_equal(get_ld_clusters(mat, .8, tidy = TRUE), tib2)

})

test_that("correlation error", {

  expect_error(get_ld_clusters(mat, -.3))
  expect_error(get_ld_clusters(mat, 3))

})

test_that("correlation error: different names", {

  mat2 <- mat
  rownames(mat2) <- letters[seq_len(nrow(mat2))]
  colnames(mat2) <- LETTERS[seq_len(ncol(mat2))]

  expect_error(get_ld_clusters(mat2, 1))

})

mat3 <- mat
rownames(mat3) <- letters[seq_len(nrow(mat3))]
colnames(mat3) <- letters[seq_len(nrow(mat3))]

snps <- c(seq(25, by = 25, length.out = 3),
  seq(150, by = 50, length.out = 2),
  seq(300, by = 50, length.out = 2))
 
names(snps) <- rownames(mat3)
out3 <- c(
    stringr::str_c("ld", rep(1, 3)),
    stringr::str_c("ld", rep(2, 2)),
    stringr::str_c("ld", rep(3, 2)))
names(out3) <- names(snps)

test_that("ld clusters with distance works vector", {

  expect_equal(
    get_ld_clusters(mat3, .8, snps, 50),
    out3)

  snps2 <- seq(25, by = 25, length.out = 7)
  names(snps2) <- rownames(mat3)

  out4 <- stringr::str_c("ld", seq_along(snps2))
  names(out4) <- names(snps2)

  expect_equal(
    get_ld_clusters(mat3, .8, snps2, 24),
    out4)

})

test_that("ld clusters with distance works IRanges", {

  snps <- IRanges::IRanges(start = snps, width = 1)
  names(snps) <- rownames(mat3)

  expect_equal(
    get_ld_clusters(mat3, .8, snps, 50),
    out3)

})

test_that("ld clusters with distance works GRanges", {

  snps <- IRanges::IRanges(start = snps, width = 1)
  names(snps) <- rownames(mat3)

  snps <- GenomicRanges::GRanges(seqnames = "chrom", ranges = snps)

  expect_equal(
    get_ld_clusters(mat3, .8, snps, 50),
    out3)

})
