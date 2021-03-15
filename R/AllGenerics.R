
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
#' @return the min. change tolerance used by the FMHighLD algorithm
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' min_tol(FMParam())
setGeneric("min_tol",
  function(object) standardGeneric("min_tol"))

#' annot_tol method
#'
#' `annot_tol` returns the minimum tolerance to accept variants with annotation
#' equal to zero
#' @param object a `FMParam` object
#' @references the min. tolerance used to accept an annotation equal to zero
#' @docType methods
#' @rdname FMParam-methods
#' @export
#' @examples
#' annot_tol(FMParam())
setGeneric("annot_tol",
  function(object) standardGeneric("annot_tol"))

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
  
#' sigma0 method
#'
#' `sigma0` obtains the std. deviations of the non-causal components
#' @param object a `FMIter` object
#' @param response the vector of all the z-values in the FMHighLD model
#' @return the std. deviation of the non-causal mixture
#' @docType methods
#' @rdname FMIter-methods
#' @export
setGeneric("sigma0",
  function(object, response) standardGeneric("sigma0"))

#' causal_candidates method
#'
#' Gets lists of causal candidates for every mixture
#' @param object a `FMIter` object
#' @return a list of causal candidates for every mixture componet
#' @docType methods
#' @rdname FMIter-methods
#' @export
setGeneric("causal_candidates",
  function(object) standardGeneric("causal_candidates"))

#' Compute the column-wise JSD metric between the gamma matrix of the current
#' and previous FMIter objects
#' @param object The current `FMIter` object
#' @param prev The previous `FMIter` iteration of the object
#' @return a vector with the jsd metric for every mixture components
#' @docType methods
#' @rdname FMIter-methods
setGeneric("prob_metric",
  function(object, prev) standardGeneric("prob_metric"))
  
#' Compute the mean one-zero loss between the causal candidates selected
#' @param object The current `FMIter` object
#' @param prev The previous `FMIter` iteration of the object
#' @return the mean one-zero loss of causal candidates variants between the
#'  current and previous iteration
#' @docType methods
#' @rdname FMIter-methods
setGeneric("mccl",
  function(object, prev) standardGeneric("mccl"))

#' Compute the difference between the current and previous beta coefficients of
#'  the underlying linear (mixed) models
#' @param object The current `FMIter` object
#' @param prev The previous `FMIter` iteration of the object
#' @param model_error a logical value indicating whehter including the model
#'  errors' difference to the output vector
#' @param background_error a logical value indicating whether including the
#'  background errors' difference to the output vector
#' @return the distance between the fixed effect coefficients
#' @docType methods
#' @rdname FMIter-methods
setGeneric("coef_diff",
  function(object, prev, model_error, background_error)
    standardGeneric("coef_diff"))
