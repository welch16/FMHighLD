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
#' @importFrom stats formula
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