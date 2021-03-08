
#' @rdname FMParam-methods
#' @aliases strategy
#' @docType methods
setMethod("strategy",
  signature = signature(object = "FMParam"),
  definition = function(object) object@strategy)

#' @rdname FMParam-methods
#' @aliases  params_all
#' @docType methods
setMethod("params_all",
  signature = signature(object = "FMParam"),
  definition = function(object)
    list(prob = object@prob, k = object@k))

#' @rdname FMParam-methods
#' @aliases params_pickm
#' @docType methods
setMethod("params_pickm",
  signature = signature(object = "FMParam"),
  definition = function(object) {
    params <- params_all(object)
    params[["m"]] <- object@m
    params[["msel"]] <- object@m_sel
    params
})

#' @rdname FMParam-methods
#' @aliases error_bound
#' @docType methods
setMethod("error_bound",
  signature = signature(object = "FMParam"),
  definition = function(object) object@error_bound)

#' @rdname FMParam-methods
#' @aliases max_iter
#' @docType methods
setMethod("max_iter",
  signature = signature(object = "FMParam"),
  definition = function(object) object@max_iter)
  
#' @rdname FMParam-methods
#' @aliases min_tol
#' @docType methods
setMethod("min_tol",
  signature = signature(object = "FMParam"),
  definition = function(object) object@min_tol)

#' @rdname FMParam-methods
#' @aliases annot_tol
#' @docType methods
setMethod("annot_tol",
  signature = signature(object = "FMParam"),
  definition = function(object) object@annot_tol)
