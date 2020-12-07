
#' @rdname FMRandParam-methods
#' @aliases strategy
#' @docType methods
setMethod("strategy",
  signature = signature(object = "FMRandParam"),
  definition = function(object) object@strategy)

#' @rdname FMRandParam-methods
#' @aliases  params_all
#' @docType methods
setMethod("params_all",
  signature = signature(object = "FMRandParam"),
  definition = function(object)
    list(prob = object@prob, k = object@k))

#' @rdname FMRandParam-methods
#' @aliases params_pickm
#' @docType methods
setMethod("params_pickm",
  signature = signature(object = "FMRandParam"),
  definition = function(object) {
    params <- params_all(object)
    params[["m"]] <- object@m
    params[["msel"]] <- object@m_sel
    params
})