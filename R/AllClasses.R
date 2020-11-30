#' @import  methods
#' @import utils
NULL


#' @rdname RandParam
setClass("RandParam",
  representation = representation(strategy = "character"),
  prototype = prototype(strategy = "none"))



#' @rdname FmIter
#' @export
#' @importClassesFrom Matrix Matrix
#' @importFrom purrr map_chr
setClass("FmIter",
  representation = representation(
    singletrait = "logical",
    models = "list",
    gamma = "Matrix",
    mu = "Matrix",
    sigma = "Matrix"),
  prototype = prototype(
    singletrait = TRUE,
    models = list(),
    gamma = Matrix::Matrix(nrow = 0, ncol = 0),
    mu = Matrix::Matrix(nrow = 0, ncol = 0),
    sigma = Matrix::Matrix(nrow = 0, ncol = 0)))

setValidity("FmIter",
  function(object) {

    if (length(object@models) == 0) {
      out <- TRUE
    } else {
      out <- all(purrr::map_chr(object@models, class) == "lm") |
        all(purrr::map_chr(object@models, class) == "lmerMod")
    }
    out

  }
)

#' @name FmIter object
#'
#' \code{FmIter} contains the elements of an iteration for the FMHighLD
#' algorithm.
#'
#' @aliases FmIter FmIter-class
#'
#' @docType class
#' @rdname FmIter
#'
NULL