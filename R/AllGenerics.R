
#' strategy method
#'
#' `strategy` returns the randomization strategy
#' @param object A `FMParam` object.
#' @return The randomization strategy to be used by `FMHighLD`
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' strategy(FMParam())
setGeneric("strategy",
  function(object) standardGeneric("strategy"))

#' params_all method
#'
#' `params_all` returns a list with the randomization probability to be used
#' by `FMHighLD` and the  absolute residual rank to select another causal
#' candidate
#' @param object A `FMParam` object.
#' @return The randomization probability to be used by `FMHighLD` and the
#' absolute residual rank to select another causal candidate
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' params_all(FMParam(strategy = "all", prob = .9, k = 2))
setGeneric("params_all",
  function(object) standardGeneric("params_all"))

#' params_pickm method
#'
#' `params_pickm` returns a list with the randomization probability to
#' be used by `FMHighLD`, the  absolute residual rank to select another
#' causal candidate, the number of ld groups to apply the randomization
#' strategy and the method used to pick those ld groups
#' @param object A `FMParam` object.
#' @return The randomization probability to be used by `FMHighLD` and the
#' absolute residual rank to select another causal candidate
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' params_pickm(FMParam(strategy = "pick_m", prob = .9, k = 2))
setGeneric("params_pickm",
  function(object) standardGeneric("params_pickm"))

#' error_bound method
#'
#' `error_bound` returns a constant indicating the min. values used to avoid
#' zero denominators
#' @param object a `FMParam` object
#' @return the error bound
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' error_bound(FMParam())
setGeneric("error_bound",
  function(object) standardGeneric("error_bound"))

#' max_iter method
#'
#' `max_iter` returns the maximum number of iteration used by the FMHighLD
#' algorithm
#' @param object a `FMParam` object
#' @return the max. number of iterations
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' max_iter(FMParam())
setGeneric("max_iter",
  function(object) standardGeneric("max_iter"))

#' min_tol method
#'
#' `min_tol` returns the minimum tolerance used by the FMHighLD
#' algorithm
#' @param object a `FMParam` object
#' @return the max. number of iterations
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' min_tol(FMParam())
setGeneric("min_tol",
  function(object) standardGeneric("min_tol"))

#' models method
#'
#' `models` return the list of models estimated in the `FMHighLD` iteration
#' @param object a `FMIter` object
#' @return a list of linear models for single trait or linear mixed models for
#'   the multitrait case
#' @docType methods
#' @rdname FMIter-methods
setGeneric("models",
  function(object) standardGeneric("models"))
  
#' probmatrix method
#'
#' `probmatrix` returns the posterior probability matrix computed in the EM
#' algorithm latest iteration
#' @param object a `FMIter` object
#' @return a `Matrix::Matrix` object with the posterior probabilities computed
#' in the latest EM iteration
#' @docType methods
#' @rdname FMIter-methods
setGeneric("probmatrix",
  function(object) standardGeneric("probmatrix"))