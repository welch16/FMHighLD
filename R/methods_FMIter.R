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

#' @rdname FMIter-methods
#' @aliases prob_metric
#' @docType methods
setMethod("prob_metric",
  signature = signature(object = "FMIter", prev = "FMIter"),
  definition = function(object, prev) {

    nmixtures <- ncol(probmatrix(object))

    jsd_vec <- purrr::map_dbl(
      seq_len(nmixtures),
      ~ jsd(probmatrix(object)[, .], probmatrix(prev)[, .]))

    jsd_vec

})

#' @rdname FMIter-methods
#' @aliases mccl
#' @docType methods
setMethod("mccl",
  signature = signature(object = "FMIter", prev = "FMIter"),
  definition = function(object, prev) {

    purrr::map2_dbl(
      causal_candidates(object), causal_candidates(prev),
      ~ mean(.x != .y))

})

#' @rdname FMIter-methods
#' @aliases coeff_diff
#' @docType methods
setMethod("coeff_diff",
  signature = signature(object = "FMIter", prev = "FMIter"),
  definition = function(object, prev) {

    if (object@singletrait) {
      current_coefs <- purrr::map(models(object), coef)
      prev_coefs <- purrr::map(models(prev), coef)

      current_errs <- purrr::map(models(object), lm_std_err)
      prev_errs <- purrr::map(models(object), lm_std_err)

      current_background_errs <- object@sigma[1, 1]
      if (nrow(prev@sigma) > 0) {
        prev_background_errs <- prev@sigma[1, 1]
      } else {
        prev_background_errs <- 0
      }

      out1 <- purrr::map2_dbl(current_coefs, prev_coefs, ~ (.x - .y)^2)
      out2 <- purrr::map2_dbl(current_errs, prev_errs, ~ (.x - .y)^2)
      out3 <- purrr::map2_dbl(current_background_errs, prev_background_errs,
        ~ (.x - .y)^2)
      out <- c(out1, out2, out3)
      # out <- sqrt(sum(out) / (1 + max(out)))


    } else {
      browser()
      object
    }
    out
})
