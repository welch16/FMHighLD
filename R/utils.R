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

  variant_names_exists <- all(
    ! is.null(names(response)),
    ! is.null(rownames(annot_matrix)),
    ! is.null(names(ld_clusters)))

  resp_annot_names <- identical(sort(names(response)),
    sort(rownames(annot_matrix)))
  resp_ldclust_names <- identical(
    sort(names(response)),
    sort(names(ld_clusters)))

  variant_names_exists & resp_annot_names & resp_ldclust_names

}
