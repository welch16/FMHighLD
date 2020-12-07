#' @import  methods
#' @import utils
NULL


#' @rdname FMRandParam
setClass("FMRandParam",
  representation = representation(
    strategy = "character",
    prob = "numeric",
    k = "numeric",
    m = "numeric",
    m_sel = "character"),
  prototype = prototype(
    strategy = NA_character_,
    prob = NA_real_,
    k = NA_integer_,
    m = NA_integer_,
    m_sel = NA_character_))
    
setValidity("FMRandParam",
  function(object) {

    if (is.na(object@strategy)) {
      out <- TRUE
    } else {
      valid_strat <- object@strategy %in% c("none", "all", "pick_m")
      if (object@strategy == "none") {
        valid_params <- TRUE
      } else if (object@strategy == "all") {
        valid_params <- object@prob > 0 & object@prob < 1 &
          object@k > 0
      } else if (object@strategy == "pick_m") {
        valid_params <- object@prob > 0 & object@prob < 1 &
          object@k > 0 & object@m > 0 &
          object@m_sel %in% c("any", "large_ld")
      }
      out <- valid_strat & valid_params
    }
    out
  }
)


#' Constructor for the `FMRandParam` class
#'
#' `FMRandParam` contains the configuration for randomized selection of causal
#' candidate variants
#'
#' @param strategy name for the randomization strategy to be used. It could be
#' one of `none`, `all` or `pick_m`. By default, if all parameters of the
#' function are `NULL` the strategy is going to be `none`.
#' For the `all` strategy, it will pick the the `k`-th variant (in terms of
#' their absolute(residual) with probability `prob`. Finally, for the `pick_m`
#' strategy, the `FMHighLD` algorithm will apply a similar strategy but for `m`
#' ld groups. The `m` ld groups will be selected randomly if `m_strat` is equal
#' to any, or like the top `m` ld groups with highest average LD score
#' @param prob probability used to randomize the causal candidate, by default
#' it is 0.95
#' @param k absolute residual index of the causal candidate to select with
#' probability `prob`.  The default this strategy will selected the variant with
#' the second rank
#' @param m number of ld groups for m_sel the randomization strategy is going to
#' be applied. By default, it will apply the strategy to only one ld group.
#' @param m_sel Selection strategy to pick the `m` ld groups to be randomized.
#' By default, it will pick the `m` groups randomly. Alternatively, it could
#' also select the `m` groups by their highest average ld score. Either one of
#' `random` or `large_ld`.
#' @return A `FMRandParam` object
#' @export
#' @importFrom methods new
#' @importFrom stringr str_to_lower
FMRandParam <- function(strategy = NULL, prob = NULL, k = NULL,
  m = NULL, m_sel = NULL) {

  if (is.null(strategy)) {
    out <- methods::new("FMRandParam", strategy = "none")
  } else {
    strategy <- stringr::str_to_lower(strategy)
    if (strategy == "none") {
      out <- methods::new("FMRandParam", strategy = "none")
    } else if (strategy == "all") {
      if (is.null(prob)) prob <- 0.95
      if (is.null(k)) k <- 2
      out <- methods::new("FMRandParam", strategy = "all", prob = prob, k = k)
    } else if (strategy == "pick_m") {
      if (is.null(prob)) prob <- 0.95
      if (is.null(k)) k <- 2
      if (is.null(m)) m <- 1
      if (is.null(m_sel)) m_sel <- "random"
      m_sel <- stringr::str_to_lower(m_sel)
      out <- methods::new("FMRandParam", strategy = "pick_m", prob = prob,
        k = k, m = m, m_sel = m_sel)
    }
  }
  out
}


#' @rdname FMIter
#' @importClassesFrom Matrix Matrix
#' @importFrom purrr map_chr
setClass("FMIter",
  representation = representation(
    nassoc = "numeric",
    singletrait = "logical",
    models = "list",
    gamma = "Matrix",
    mu = "Matrix",
    sigma = "Matrix"),
  prototype = prototype(
    nassoc = 0,
    singletrait = TRUE,
    models = list(),
    gamma = Matrix::Matrix(nrow = 0, ncol = 0),
    mu = Matrix::Matrix(nrow = 0, ncol = 0),
    sigma = Matrix::Matrix(nrow = 0, ncol = 0)))

setValidity("FMIter",
  function(object) {

    if (length(object@models) == 0) {
      out <- TRUE
    } else {
      out <- all(purrr::map_chr(object@models, class) == "lm") |
        all(purrr::map_chr(object@models, class) == "lmerMod")
      out <- out & object@nassoc == nrow(gamma) & object@nassoc == nrow(mu)
    }
    out

  }
)

#' `FMIter` object
#'
#' `FMIter` contains the elements of an iteration for the FMHighLD algorithm
#' @param nassoc a numeric variable with the number of association tested
#'  in the experiment
#' @param singletrait a logical indicator determininig if we are considering a
#'  single trait or multiple traits
#' @param models a list with the models in the mixture
#' @param gamma a matrix with the posterior probabilities of a snps belonging to
#'  any element in the mixture
#' @param mu a matrix with the predicted values of every snp with each column
#'  being an element in the mixture
#' @param sigma a matrix with the variances of every normal distribution in the
#'  mixture.
#'
#' @aliases FMIter FMIter-class
#'
#' @docType class
#' @rdname FMIter
#' @importFrom methods new
FMIter <- function(nassoc = NULL, singletrait = NULL, models = NULL,
  gamma = NULL, mu = NULL, sigma = NULL) {

  if (is.null(nassoc)) nassoc <- 0
  if (is.null(singletrait)) singletrait <- TRUE
  if (is.null(models)) models <- list()
  if (is.null(gamma)) gamma <- Matrix::Matrix(nrow = 0, ncol = 0)
  if (is.null(mu)) mu <- Matrix::Matrix(nrow = 0, ncol = 0)
  if (is.null(sigma)) sigma <- Matrix::Matrix(nrow = 0, ncol = 0)

  methods::new("FMIter",
    nassoc = nassoc, singletrait = singletrait, models = models,
      gamma = gamma, mu = mu, sigma = sigma)

}
