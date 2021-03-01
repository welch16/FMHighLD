#' @rdname FMIter-methods
#' @aliases models
#' @docType methods
setMethod("models",
  signature = signature(object = "FMIter"),
  definition = function(object) object@models)

#' @rdname FMIter-methods
#' @aliases probmatrix
#' @docType methods
setMethod("probmatrix",
  signature = signature(object = "FMIter"),
  definition = function(object) object@gamma)

#' @rdname FMIter-methods
#' @aliases sigma0
#' @docType methods
setMethod("sigma0",
  signature = signature(object = "FMIter"),
  definition = function(object, response) {

    gamma <- probmatrix(object)
    sumgamma <- colSums(as.matrix(gamma))

    causal_candidates <- causal_candidates(object)
    response_list <- purrr::map(causal_candidates, ~ response[.])
    out_list <- purrr::map_dbl(response_list,
      ~ crossprod(.^2, gamma[, 1]) / sumgamma[1])
    sqrt(mean(out_list))

})

#' @rdname FMIter-methods
#' @aliases causal_candidates
#' @docType methods
setMethod("causal_candidates",
  signature = signature(object = "FMIter"),
  definition = function(object)object@causal_candidates)
