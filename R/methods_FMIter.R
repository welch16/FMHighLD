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


substitute_weights <- function(i, weight, causal_candidate, gamma) {

  gamma_vec <- rlang::set_names(as.vector(gamma[, i]), causal_candidate)
  weight[causal_candidate] <- gamma_vec
  weight

}

#' @rdname FMIter-methods
#' @aliases loglikelihood
#' @docType methods
setMethod("loglikelihood",
  signature = signature(object = "FMIter", fmld_data = "data.frame"),
  definition = function(object, fmld_data) {

    if (object@singletrait) {
      fmld_data <- as.data.frame(fmld_data)
      fmld_data <- tibble::column_to_rownames(fmld_data, "snp")
      preds <- purrr::map(models(object),
        stats::predict, fmld_data)
      sigma <- object@sigma[1, ]
      f0 <- stats::dnorm(x = fmld_data[["response"]],
        mean = 0, sd = sqrt(sigma[1]))
      fz <- purrr::map2(preds, sigma[-1],
        ~ stats::dnorm(x = fmld_data[["response"]],
        mean = .x, sd = sqrt(.y)))
      probs <- compute_mixture_prob(object@gamma)
      f0 <- probs[1] * f0
      fz <- purrr::map2(fz, probs[-1], `*`)
      fz[[length(fz) + 1]] <- f0
      fz <- do.call(cbind, fz)
      loglike <- sum(log(rowSums(fz)))

    } else {

      browser()
      x

    }

    loglike
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
#' @aliases coef_diff
#' @docType methods
setMethod("coef_diff",
  signature = signature(object = "FMIter", prev = "FMIter",
    model_error = "logical", background_error = "logical"),
  definition = function(object, prev, model_error, background_error) {

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

      out1 <- purrr::map2_dbl(current_coefs, prev_coefs, ~ sum((.x - .y)^2))
      out <- out1
      if (model_error) {
        out2 <- purrr::map2_dbl(current_errs, prev_errs, ~ sum((.x - .y)^2))
        out <- c(out, out2)
      }
      if (background_error) {
        out3 <- purrr::map2_dbl(current_background_errs, prev_background_errs,
          ~ (.x - .y)^2)
        out <- c(out, out3)
      }

    } else {
      browser()
      object
    }
    out
})

#' @rdname FMIter-methods
#' @aliases coef_diff
#' @docType methods
#' @importFrom flexmix parameters
setMethod("coef_diff_em",
  signature = signature(object = "FMIter", prev = "FMIter",
    model_error = "logical", background_error = "logical"),
  definition = function(object, prev, model_error, background_error) {

    if (object@singletrait) {

      current_params <- purrr::map(models(object), flexmix::parameters)
      prev_params <- purrr::map(models(prev), flexmix::parameters)
      out1 <- purrr::map2_dbl(current_params, prev_params,
        ~ sum((.x[[1]][1:2, 1] - .y[[1]][1:2, 1])^2))
      out <- out1
      if (model_error) {
        out2 <- purrr::map2_dbl(current_params, prev_params,
          ~ sum((.x[[1]][3, 1] - .y[[1]][3, 1])^2))
        out <- c(out, out2)
      }
      if (background_error) {
        out3 <- purrr::map2_dbl(current_params, prev_params,
          ~ (.x[[2]][2, 2] - .y[[2]][2, 2])^2)
        out <- c(out, out3)
      }

    } else {
      browser()
      object
    }
    out
})
