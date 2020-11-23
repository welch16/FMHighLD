#' @import  methods
#' @import utils
NULL


#' @rdname RandParam
#' @export
setClass("RandParam",
  representation = representation(strategy = "character"),
  prototype = prototype(strategy = "none"))

setValidity("RandParam",
  function(object) {
  })
