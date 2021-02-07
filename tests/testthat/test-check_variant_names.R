 nsnps <- 10
 snps <- stringr::str_c("snp", seq_len(nsnps))
 z <- rnorm(nsnps)
 names(z) <- snps
 z_list <- list(
   t1 = z[1:4],
   t2 = z[5:6],
   t3 = z[7:10],
   t4 = z[3:8])

annot <- matrix(runif(nsnps * 5), nrow = nsnps)
rownames(annot) <- snps
colnames(annot) <- stringr::str_c("annot", seq_len(5))

ldmat <- noise_cor(diag(nsnps), epsilon = .5, eidim = 2)
rownames(ldmat) <- snps
colnames(ldmat) <- snps
ld <- get_ld_clusters(ldmat, min_r2 = .1)

test_that("expect error if any argument is null", {

  expect_error(check_variant_names(response = z, annot_matrix = annot))
  expect_error(check_variant_names(response = z_list, annot_matrix = annot))
  expect_error(check_variant_names(annot_matrix = annot, ld_clusters = ld))
  expect_error(check_variant_names(response = z, ld_clusters = ld))

})

test_that("expect `FALSE` if response is unnamed", {

  expect_equal(check_variant_names(response = set_names(z, NULL),
    annot_matrix = annot, ld_clusters = ld), FALSE)
  aux_z_list <- purrr::map(z_list, set_names, NULL)
  expect_equal(check_variant_names(response = aux_z_list,
    annot, ld), FALSE)

})

test_that("expect `FALSE` if annot rows are unnamed", {

  aux_annot <- annot
  rownames(aux_annot) <- NULL
  expect_equal(check_variant_names(response = z,
    annot_matrix = aux_annot, ld_clusters = ld), FALSE)

})

test_that("expect `FALSE` if ld is unnamed", {

  expect_equal(check_variant_names(response = z,
    annot_matrix = annot, ld_clusters = set_names(ld, NULL)), FALSE)

})

test_that("expect `FALSE` if a name is missing in multitrait", {

  aux_z_list <- z_list
  names(aux_z_list[[2]]) <- NULL

  expect_equal(check_variant_names(response = aux_z_list,
    annot, ld), FALSE)

})

test_that("expect `TRUE` if all names are correct", {

  expect_equal(check_variant_names(z, annot, ld), TRUE)
  expect_equal(check_variant_names(z_list, annot, ld), TRUE)

})
