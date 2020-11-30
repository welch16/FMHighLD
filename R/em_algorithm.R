#' Performs an Expectation-Maximization iterations for the multi trait model
#' @param formula a formula object with the underlying linear model
#' @param data a `data.frame` object with the annotation matrix and z-values
#' @param pi a vector of probabilities in the FMHighLD mixture part
#' @param model_list a list with the underlying linear models of FMHighLD
#' @param sigma0 the error variance estimate of the background model
#' @param verbose a logical indicator of adding process messages
#' @return a list with the updated parameters in one EM iteration
#' @importFrom purrr map map_int
#' @importFrom lme4 lFormula
#' @importFrom rlang rep_along
#' @importFrom dplyr filter mutate
#' @importFrom Matrix crossprod Matrix
#' @importFrom stats model.matrix
em_iteration_multi <-
  function(formula, data, pi, model_list, sigma0, verbose = TRUE) {

  if (verbose) message("--starting EM iter")

  error_bound <- 1e-4
  # N <- nrow(data)

  # include_intercept <- attr(stats::terms(formula), "intercept")
  response <- get_response_name(formula)

  # convert pi to probabilities
  pi <- pi / sum(pi)

  # generate mean matrix
  if (verbose) message("--calculating means")
  mu <- purrr::map(model_list, stats::predict, newdata = data)
  names(mu) <- NULL
  mu <- do.call(cbind, mu)
  # add background mean
  mu <- cbind(0, mu)

  # generate variance matrix
  if (verbose) message("--calculating variances")

  Zt <- lme4::lFormula(formula, data)[["reTrms"]][["Zt"]]
  re_error_list <- purrr::map(model_list, broom::tidy)

  # fix zero or almost zero estimated variances
  residual_error_list <- purrr::map(
    re_error_list, dplyr::filter, group == "Residual")
  residual_error_list <- purrr::map(residual_error_list, "estimate")
  residual_error_list <- purrr::map(residual_error_list,
    ~ pmax(., rlang::rep_along(., error_bound)))

  re_error_list <- purrr::map(re_error_list, dplyr::filter,
    !group %in% c("fixed", "Residual"))
  re_error_list <- purrr::map(re_error_list, "estimate")
  re_error_list <- purrr::map(re_error_list,
    ~ pmax(., rep_along(., error_bound)))

  nrandom_effects <- unique(purrr::map_int(re_error_list, length))
  ntraits <- nrow(Zt) / nrandom_effects

  re_variances <- purrr::map(re_error_list, parse_re_variance, ntraits)
  sigma <- purrr::map(re_variances, ~ Matrix::crossprod(Zt, .))

  sigma <- purrr::map(sigma, ~ . %*% Zt)
  sigma <- purrr::map(sigma, diag)

  sigma <- purrr::map2(sigma, residual_error_list, ~ sqrt(.x + .y ^ 2))
  sigma <- do.call(cbind, sigma)
  sigma <- cbind(sigma0, sigma)

  if (verbose) message("--calculating weights")

  gamma_mat <- estep(
    purrr::pluck(data, response), as.vector(pi), as.matrix(mu), sigma)
  gamma_mat <- pmax(gamma_mat, 1e-12)
  idxs <- seq_len(ncol(gamma_mat))

  X_matrices <- purrr::map(model_list, stats::model.matrix)

  # p = X_matrices %>% map_int(ncol) %>% unique()

  if (verbose) message("--fitting new models")
  models <- purrr::map2( idxs[-1], X_matrices, underlying_lme_model,
    gamma_mat, error_bound)

  group <- NULL
  list(
    "models" = models,
    "gamma" = gamma_mat,
    "mu" = mu,
    "sigma" = sigma)
}

#' Gets the name of the response based on the formula
#' @param formula a formula object with the underlying linear / linear mixed
#' model
#' @return the name of the response
#' @importFrom stats formula
get_response_name <- function(formula) {

out <- attr(stats::terms(formula), "factors")
rownames(out)[1]

}

#' Parses the random error variance matrix
#' @param re_error estimates of the RE errors
#' @param ntraits integer with the number of trais
#' @return a sparse `Matrix` object
#' @importFrom Matrix Matrix
parse_re_variance <- function(re_error, ntraits) {

  re_var <- re_error ^ 2
  re_vars <- rep(re_var, ntraits)
  Matrix::Matrix(diag(re_vars), sparse = TRUE)
}

#' Performs an Expectation-Maximization iterations for the single
#' trait model
#' @param formula a formula object with the underlying linear model
#' @param data a `data.frame` object with the annotation matrix and z-values
#' @param pi a vector of probabilities in the FMHighLD mixture part
#' @param model_list a list with the underlying linear models of FMHighLD
#' @param sigma0 the error variance estimate of the background model
#' @param verbose a logical indicator of adding process messages
#' @return a list with the updated parameters in one EM iteration
#' @importFrom purrr map map2 map_dbl map_int
#' @importFrom stats predict update model.matrix
#' @importFrom broom glance
em_iteration_single <- function(
  formula, data, pi, model_list, sigma0, verbose) {

  if (verbose) message("--starting EM iter")

  error_bound <- 1e-4
  N <- nrow(data)

  # build a mean matrix to build E-step probabilities
  # include_intercept <- attr(stats::terms(formula), "intercept")
  response <- get_response_name(formula)

  # convert to probabilities if they are not
  pi <- pi / sum(pi)

  # generate mean matrix

  if (verbose) message("--calculating means")
  mu <- purrr::map(model_list, stats::predict, newdata = data)
  names(mu) <- NULL
  mu <- do.call(cbind, mu)
  # add background mean
  mu <-  cbind(0, mu)

  # generate variance matrix
  if (verbose) message("--calculating variances")
  sigma <- purrr::map(model_list, broom::glance)
  sigma <- purrr::map_dbl("sigma")

  # need to convert into a matrix for use with Rcpp estep function
  sigma <- purrr::map(sigma, rep, N)
  sigma <- cbind(sigma0, sigma)

  if (verbose) message("--calculating weights")

  gamma_mat <- estep(
    purrr::pluck(data, response), as.vector(pi), as.matrix(mu), sigma)
  gamma_mat <- pmax(gamma_mat, 1e-12)
  idxs <- seq_len(ncol(gamma_mat))

  X_matrices <- purrr::map(model_list,
    ~ stats::update(., data = mutate(data, w = 1)))
  X_matrices <- purrr::map(X_matrices, stats::model.matrix)

  # p <- unique(purrr::map_int(X_matrices, ncol))

  if (verbose) message("--fitting new models")
  models <- purrr::map2( idxs[-1], X_matrices, underlying_linear_model,
    gamma_mat, error_bound)

  list(
    "models" = models,
    "gamma" = gamma_mat,
    "mu" = mu,
    "sigma" = sigma)
}

#' performs the underlying linear model in FMHighLD
#' @param i index of the model
#' @param X covariate matrix
#' @param gamma estep probabilities
#' @param error_bound numerical constant to avoid zero eigenvalues
#' @return an `lm` model
#' @importFrom dplyr mutate
#' @importFrom stats lm
underlying_linear_model <- function(i, X, gamma, error_bound) {

  w <- NULL
  weights <- gamma[, i]
  eigvals <- eigen(crossprod(X, diag(weights)) %*% X)

  if (any(eigvals$values <= 0)){
    # fix very little eig values, which may cause numerical difficulties
    weights <- pmax(weights, error_bound)
  }

  stats::lm(formula, data = dplyr::mutate(data, w = weights),
    weights = w)
}


#' performs the underlying linear model in FMHighLD
#' @param i index of the model
#' @param X covariate matrix
#' @param gamma estep probabilities
#' @param error_bound numerical constant to avoid zero eigenvalues
#' @return an `lmer` model
#' @importFrom dplyr mutate
#' @importFrom lme4 lmer
#' @importFrom Matrix crossprod
underlying_lme_model <- function(i, X, gamma, error_bound){

  weights <- gamma[, i]
  eigvals <- eigen(Matrix::crossprod(X, diag(weights)) %*% X)
  if (any(eigvals$values <= 0)) {
    weights <- pmax(weights, error_bound)
  }
  lme4::lmer(formula, data, weights = weights)

}
