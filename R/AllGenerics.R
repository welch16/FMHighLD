
#' strategy method
#'
#' `strategy` returns the randomization strategy
#' @param object A `FMRandParam` object.
#' @return The randomization strategy to be used by `FMHighLD`
#' @docType methods
#' @rdname FMRandParam-methods
setGeneric("strategy",
  function(object) standardGeneric("strategy"))

#' params_all method
#'
#' `params_all` returns a list with the randomization probability to be used
#' by `FMHighLD` and the  absolute residual rank to select another causal
#' candidate
#' @param object A `FMRandParam` object.
#' @return The randomization probability to be used by `FMHighLD` and the
#' absolute residual rank to select another causal candidate
#' @docType methods
#' @rdname FMRandParam-methods
setGeneric("params_all",
  function(object) standardGeneric("params_all"))

#' params_pickm method
#'
#' `params_pickm` returns a list with the randomization probability to
#' be used by `FMHighLD`, the  absolute residual rank to select another
#' causal candidate, the number of ld groups to apply the randomization
#' strategy and the method used to pick those ld groups
#' @param object A `FMRandParam` object.
#' @return The randomization probability to be used by `FMHighLD` and the
#' absolute residual rank to select another causal candidate
#' @docType methods
#' @rdname FMRandParam-methods
setGeneric("params_pickm",
  function(object) standardGeneric("params_pickm"))
