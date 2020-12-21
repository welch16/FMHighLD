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
