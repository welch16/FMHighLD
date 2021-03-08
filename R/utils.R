#' Compute the Kullback-Leibler (KL) divergence between two discrete
#' probability distributions
#' @param p1 the first probability vector
#' @param p2 the second probability vector
#' @return a scalar value with the KL-divergence between both
#'  distributions
#' @export
#' @examples
#'
#' u1 <- runif(100)
#' u1 <- u1 / sum(u1)
#' v1 <- runif(100)
#' v1 <- v1 / sum(v1)
#' kl(u1, v1)
kl <- function(p1, p2) {
  p1 <- p1 / sum(p1)
  p2 <- p2 / sum(p2)
  sum(p1 * (log(p1) - log(p2)))
}

#' Computes the Jensen-Shannon (JS) divergence between two discrete
#' probability distributions
#' @param p1 the first probability vector
#' @param p2 the second probability vector
#' @return a scalar value with the JS-divergence between both
#'  distributions
#' @export
#' @examples
#'
#' u1 <- runif(100)
#' u1 <- u1 / sum(u1)
#' v1 <- runif(100)
#' v1 <- v1 / sum(v1)
#' jsd(u1, v1)
jsd <- function(p1, p2) {
  m <- 0.5 * (p1 + p2)
  0.5 * kl(p1, m) + 0.5 * kl(p2, m)
}


#' Checks that the variant names exists and are the same across the response,
#'  the annotation matrix and the ld_clusters vector
#' @param response a vector for the `singletrait` case or a named list for
#'  the multitrait case with an element for each variant. In the multitrait
#'  case, each element of the list is a named vector where each element
#'  is named after the trait for which the association was tested.
#' @param annot_matrix a matrix with same number of rows as variants and number
#'  of columns as annotations
#' @param ld_clusters a character vector with the ld cluster to which each
#'  variant belongs
#' @return a logical indicator determining if the response, annotation matrix
#'  and ld_clusters variables are named, and whether the names are the same
#' @export
check_variant_names <- function(response, annot_matrix, ld_clusters) {

  if (is.list(response)) {
    response_names <- purrr::map(response, names)
    response_names <- unique(unlist(response_names))
    have_name <- all(!purrr::map_lgl(response, ~ is.null(names(.))))
  } else {
    response_names <- names(response)
    have_name <- TRUE
  }

  annot_names <- rownames(annot_matrix)
  ld_names <- names(ld_clusters)

  variant_names_exists <- all(
    ! is.null(response_names),
    ! is.null(annot_names),
    ! is.null(ld_names))

  resp_annot_names <- identical(sort(response_names),
    sort(annot_names))

  resp_ldclust_names <- identical(
    sort(response_names),
    sort(ld_names))

  have_name & variant_names_exists & resp_annot_names & resp_ldclust_names

}

#' Gets the name of the response based on the formula
#' @param formula a formula object with the underlying linear / linear mixed
#' model
#' @export
#' @return the name of the response
#' @importFrom stats formula as.formula
#' @examples
#' get_response_name(y ~ x)
#' get_response_name("y ~ x")
get_response_name <- function(formula) {

if (is.character(formula)) formula <- as.formula(formula)

out <- attr(stats::terms(formula), "factors")
if (all(rowSums(out) > 0)) {
  stop("There is not response in formula")
}

rownames(out)[1]

}

#' Checks whether the formula contains the intercept
#'
#' @param formula a lm / lmm `formula` for the underlying linear models
#' @return a logical value indicating whether the fixed effects have an
#'  intercept term
#' @export
#' @importFrom stats terms
has_intercept <- function(formula) {

  stopifnot(class(formula) == "formula")

  attr(stats::terms(formula), "intercept") == 1

}

#' extract causal candidate vectors
#'
#' @param causal_list a list of `tibble::tibble` with columns ld_cluster and
#'  which_snp
#' @return a list of the same length as `causal_list` with two names character
#'  vectors with the causal candidate snps
#' @examples
#' causal_list <- list(tibble::tibble(
#'  ld_cluster = stringr::str_c("ld", 1:2),
#'  which_snp = stringr::str_c("snp", 1:2)))
#' extract_causal_vector(causal_list)
#' @importFrom rlang set_names
#' @importFrom purrr map map2 map_lgl
#' @importFrom tibble is_tibble
extract_causal_vector <- function(causal_list) {

  stopifnot(is.list(causal_list))
  stopifnot(all(purrr::map_lgl(causal_list, tibble::is_tibble)))
  all_names <- purrr::map(causal_list, names)
  stopifnot(all(purrr::map_lgl(all_names, ~ "which_snp" %in% .)))
  stopifnot(all(purrr::map_lgl(all_names, ~ "ld_cluster" %in% .)))

  candidates <- purrr::map(causal_list, "which_snp")
  ld_clusters <- purrr::map(causal_list, "ld_cluster")
  purrr::map2(candidates, ld_clusters, rlang::set_names)

}

#' Wrap around the `norm` operator to work with `dplyr::rowwise`
#'
#' @param ... names of the variables to use when computing the norm2
#' @return the norm2 of the columns picks
#' @examples
#' library(magrittr)
#' tibble::tribble( ~ "vec1", ~ "vec2", 1, 1, 0, 1, 1, 0, 0, 0) %>%
#'  dplyr::rowwise() %>%
#'  dplyr::mutate( norm2_wrap(vec1, vec2))
norm2_wrap <- function(...) {

  vec <- list(...)
  vec <- unlist(vec)
  norm(matrix(vec), "F")

}
