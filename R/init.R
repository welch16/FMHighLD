
#' Computes the i-th mixture model for fitting the FM-HighLD model
#' @param i a numeric index, usually i = 2 only that correspond to the causal
#'  variants, but the model is performed for a list in case more models are
#'  being fitted for the locus
#' @param formula formula used for fitting the underlying linear model used by
#'  FM-HighLD
#' @param data `data.frame` used for the FM-HighLD model
#' @param gamma a `Matrix::Matrix` with the conditional probabilities of the
#'  variants being in each model
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @importFrom stats lm
#' @importFrom lme4 lmer
compute_ith_model <- function(i, formula, data, gamma, singletrait = TRUE) {

  if (singletrait) {
    data$weights <- gamma[, i]
    model <- stats::lm(formula, data = data, weights = weights)
  } else {
    weights <- gamma[, i]
    model <- lme4::lmer(formula, data = data, weights = weights)
  }
  model
}


#' Performs the initial iteration of the FM-HighLD model
#' @param formula formula used for fitting the underlying linear model used by
#'  FM-HighLD
#' @param data `data.frame` used for the FM-HighLD model
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @return a list with the models, and probability matrix for the EM-algorithm
#' @importFrom Matrix Matrix
init_iteration <- function(formula, data, singletrait = TRUE) {

  n <- nrow(data)
  gamma_mat <- do.call(rbind, lapply(seq_len(n), function(x)c(1, 1)))
  gamma_mat <- Matrix::Matrix(gamma_mat)

  p <- ncol(gamma_mat)
  idxs <- seq_len(p)

  models <- lapply(idxs[-1],
    compute_ith_model,
    formula, data, gamma, multi_trait)

  FMIter(nassoc = nrow(data), singletrait = singletrait, models = models,
    gamma = gamma_mat)

}

#' Utility to build the formula of the model of the underlying FM-HighLD model
#' @param response name of the response variable
#' @param fixed_effects names of the fixed effects of the model
#' @param fixed_intercept logical indicating if there is a fixed effect
#'  intercept, by default this values is `TRUE`
#' @param random_effect_group a character variable or NULL indicating if there
#'  are random effects for the multi_trait case
#' @param random_effects names of the random effects of the model
#' @param random_intercept logical indicating if there is a random effect
#'  intercept, by default this values is `TRUE`
#' @export
#' @examples
#' build_formula("z", c("skew", "var1", "var2"), FALSE)
#' build_formula("z", c("skew", "var1", "var2"), TRUE)
#' build_formula("z", c("skew", "var1", "var2"), TRUE, "gene", c("skew"), TRUE)
#' build_formula("z", c("skew", "var1", "var2"), TRUE, "gene", c("skew"), FALSE)
#' @importFrom dplyr if_else
#' @importFrom glue glue
#' @importFrom stringr str_c
build_formula <- function(response, fixed_effects, fixed_intercept,
  random_effect_group = NULL, random_effects, random_intercept) {

  collapse_sign <- "+"

  if (is.null(random_effect_group)) {
    out_formula <- glue::glue(
      "{response} ~ {fixint} + {fixed}",
      response = response,
      fixint = dplyr::if_else(fixed_intercept, "1", "0"),
      fixed = stringr::str_c(fixed_effects, collapse = collapse_sign))
  } else {
    out_formula <- glue::glue(
      "{response} ~ {fixint} + {fixed} + ({ranint} + {random} || {group})",
      response = response,
      fixint = dplyr::if_else(fixed_intercept, "1", "0"),
      fixed = stringr::str_c(fixed_effects, collapse = collapse_sign),
      ranint = dplyr::if_else(random_intercept, "1", "0"),
      random = stringr::str_c(random_effects, collapse = collapse_sign),
      group = random_effect_group)
  }
  as.character(out_formula)

}
