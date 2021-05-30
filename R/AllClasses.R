#' @import  methods
#' @import utils
NULL


#' @rdname FMParam
setClass("FMParam",
  representation = representation(
    strategy = "character",
    prob = "numeric",
    k = "numeric",
    m = "numeric",
    m_sel = "character",
    error_bound = "numeric",
    max_iter = "numeric",
    min_tol = "numeric",
    annot_tol = "numeric"),
  prototype = prototype(
    strategy = NA_character_,
    prob = NA_real_,
    k = NA_integer_,
    m = NA_integer_,
    m_sel = NA_character_,
    error_bound = 1e-4,
    max_iter = 100,
    min_tol = 1e-6,
    annot_tol = 1e-6))
    
setValidity("FMParam",
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
          object@m_sel %in% c("random", "large_ld")
      }
      out <- valid_strat & valid_params
    }
    out <- out & object@error_bound > 0 & object@error_bound < 1 &
      object@max_iter > 1 & object@min_tol > 0 & object@min_tol < 1
    out
  }
)


#' Constructor for the `FMParam` class
#'
#' `FMParam` contains the configuration for randomized selection of causal
#' candidate variants
#'
#' @param strategy name for the randomization strategy to be used. It could be
#' one of `none`, `all` or `pick_m`. By default, if all parameters of the
#' function are `NULL` the strategy is going to be `none`.
#' For the `all` strategy, it will pick the the `k`-th variant (in terms of
#' their absolute(residual) with probability `prob`. Finally, for the `pick_m`
#' strategy, the `FMHighLD` algorithm will apply a similar strategy but for `m`
#' LD groups. The `m` LD groups will be selected randomly if `m_strat` is equal
#' to any, or like the top `m` LD groups with highest average LD score
#' @param prob probability used to randomize the causal candidate, by default
#' it is 0.95
#' @param k absolute residual index of the causal candidate to select with
#' probability `prob`.  The default this strategy will selected the variant with
#' the second rank
#' @param m number of LD groups for m_sel the randomization strategy is going to
#' be applied. By default, it will apply the strategy to only one LD group.
#' @param m_sel Selection strategy to pick the `m` LD groups to be randomized.
#' By default, it will pick the `m` groups randomly. Alternatively, it could
#' also select the `m` groups by their highest average LD score. Either one of
#' `random` or `large_ld`.
#' @param error_bound a double constant indicating the min value used to avoid
#'  zero denominators
#' @param max_iter the max. number of iteration by the algorithm used to stop
#'  the algorithm
#' @param min_tol the min. tolerance used to stop the algorithm
#' @param annot_tol the min. tolerance to accept prediction with small
#'  annotation values
#' @return A `FMParam` object with the `FMHighLD` algorithm configuration
#' @export
#' @importFrom methods new
#' @importFrom stringr str_to_lower
FMParam <- function(strategy = NULL, prob = NULL, k = NULL,
  m = NULL, m_sel = NULL, error_bound = 1e-4, max_iter = 100,
  min_tol = 1e-6, annot_tol = 1e-6) {

  if (is.null(strategy)) {
    out <- methods::new("FMParam", strategy = "none",
        error_bound = error_bound, max_iter = max_iter, min_tol = min_tol)
  } else {
    strategy <- stringr::str_to_lower(strategy)
    if (strategy == "none") {
      out <- methods::new("FMParam", strategy = "none",
        error_bound = error_bound, max_iter = max_iter, min_tol = min_tol)
    } else if (strategy == "all") {
      if (is.null(prob)) prob <- 0.95
      if (is.null(k)) k <- 2
      out <- methods::new("FMParam", strategy = "all", prob = prob, k = k,
        error_bound = error_bound, max_iter = max_iter, min_tol = min_tol)
    } else if (strategy == "pick_m") {
      if (is.null(prob)) prob <- 0.95
      if (is.null(k)) k <- 2
      if (is.null(m)) m <- 1
      if (is.null(m_sel)) m_sel <- "random"
      m_sel <- stringr::str_to_lower(m_sel)
      out <- methods::new("FMParam", strategy = "pick_m", prob = prob,
        k = k, m = m, m_sel = m_sel, error_bound = error_bound,
        max_iter = max_iter, min_tol = min_tol, annot_tol = annot_tol)
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
    causal_candidates = "list",
    gamma = "Matrix",
    mu = "Matrix",
    sigma = "Matrix"),
  prototype = prototype(
    nassoc = 0,
    singletrait = TRUE,
    models = list(),
    causal_candidates = list(),
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
      out <- out & object@nassoc == nrow(gamma)
    }
    out

  }
)

#' `FMIter` object
#'
#' `FMIter` contains the elements of an iteration for the FMHighLD algorithm
#' @param nassoc a numeric variable with the number of association tested
#'  in the experiment
#' @param singletrait a logical indicator determining if we are considering a
#'  single trait or multiple traits
#' @param models a list with the models in the mixture
#' @param causal_candidates a list of the same length as `models` made with
#'  vectors of causal candidate SNPs per each valid trait - LD cluster
#'  combination
#' @param gamma a matrix with the posterior probabilities of a SNPs belonging to
#'  any element in the mixture
#' @param mu a matrix with the predicted values of every SNP with each column
#'  being an element in the mixture
#' @param sigma a matrix with the variances of every normal distributions in the
#'  mixture.
#'
#' @aliases FMIter FMIter-class
#'
#' @docType class
#' @rdname FMIter
#' @importFrom methods new
FMIter <- function(nassoc = NULL, singletrait = NULL, models = NULL,
  causal_candidates = NULL, gamma = NULL, mu = NULL, sigma = NULL) {

  if (is.null(nassoc)) nassoc <- 0
  if (is.null(singletrait)) singletrait <- TRUE
  if (is.null(models)) models <- list()
  if (is.null(causal_candidates)) causal_candidates <- list()
  if (is.null(gamma)) gamma <- Matrix::Matrix(nrow = 0, ncol = 0)
  if (is.null(mu)) mu <- Matrix::Matrix(nrow = 0, ncol = 0)
  if (is.null(sigma)) sigma <- Matrix::Matrix(nrow = 0, ncol = 0)

  methods::new("FMIter",
    nassoc = nassoc, singletrait = singletrait, models = models,
      causal_candidates = causal_candidates, gamma = gamma, mu = mu,
      sigma = sigma)

}
