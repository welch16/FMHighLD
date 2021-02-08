test_that("get_snp_names works", {

  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  snp_names <- sort(stringr::str_c("snp", seq_len(nsnps)))
  expect_identical(snp_names, get_snp_names(z, TRUE))
  z_list <- list()
  z_list[["a"]] <- z[1:3]
  z_list[["b"]] <- z[3:5]
  z_list[["c"]] <- z[6:10]
  expect_identical(snp_names, get_snp_names(z_list,  FALSE))

})

test_that("build_fm_tibble works", {

  library(magrittr)
  library(dplyr)
  set.seed(123)
  nsnps <- 10
  z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
  annot <- matrix(rnorm(nsnps), ncol = 1)
  rownames(annot) <- names(z)
  colnames(annot) <- "annot"
  ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
  names(ld_cluster) <- names(z)

  fm_tib <- tibble::tibble(snp = names(z), response = z,
    annot = annot[, 1], ld_cluster = ld_cluster)

  expect_equivalent(
    dplyr::arrange(fm_tib, snp), build_fm_tibble(z, annot, ld_cluster, TRUE))

  z_list <- list()
  z_list[["a"]] <- z[1:3]
  z_list[["b"]] <- z[3:5]
  z_list[["c"]] <- z[6:10]

  fm_tib2 <- purrr::map(z_list,
    ~ tibble::tibble(snp = names(.), response = .)) %>%
    purrr::map2(names(.), ~ dplyr::mutate(.x, trait = .y)) %>%
    dplyr::bind_rows()

  fm_tib2 <- dplyr::inner_join(fm_tib2, dplyr::select(fm_tib, -response),
    by = "snp")

  expect_equivalent(
    dplyr::arrange(fm_tib2, snp),
    build_fm_tibble(z_list, annot, ld_cluster, FALSE))

})