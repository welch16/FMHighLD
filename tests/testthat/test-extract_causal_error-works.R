test_that("error if input is not a list", {

  expect_error(extract_causal_vector("a"))
})

test_that("error if not all the arguments are `tibble::tibble`", {

  causal_list <- list(tibble::tibble(), "a")
  expect_error(extract_causal_vector(causal_list))

})

test_that(
  "error if one `tibble` doesn't contains `which_snp` and `ld_cluster` columns", {

  causal_list <- list(
    tibble::tibble(which_snp = "snp1", ld_cluster = "ld1"),
    tibble::tibble(which_snp = "snp1", ld_cluster0 = "ld1"))
  expect_error(extract_causal_vector(causal_list))

})

test_that("extract_causal_vector works", {

  causal_list <- list(tibble::tibble(
    ld_cluster = stringr::str_c("ld", 1:2),
    which_snp = stringr::str_c("snp", 1:2)))
  expect_identical(extract_causal_vector(causal_list),
    list(c(ld1 = "snp1", ld2 = "snp2")))

})