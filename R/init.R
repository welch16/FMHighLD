
#' Computes the i-th mixture model for fitting the FM-HighLD model
#' @param i a numeric index, usually i = 2 only that correspond to the causal
#'  variants, but the model is performed for a list in case more models are
#'  being fitted for the locus
#' @param causal_candidates a vector of the causal candidates selected for each
#'  LD cluster
#' @param formula formula used for fitting the underlying linear model used by
#'  FM-HighLD
#' @param data `data.frame` used for the FM-HighLD model
#' @param gamma a `Matrix::Matrix` with the conditional probabilities of the
#'  variants being in each model
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @importFrom stats lm
#' @importFrom lme4 lmer
compute_ith_model <- function(i, causal_candidates, formula, data, gamma,
  singletrait = TRUE) {

  data <- dplyr::filter(data, snp %in% causal_candidates)

  if (singletrait) {
    data$weights <- gamma[, i]
    model <- stats::lm(formula, data = data, weights = weights)
  } else {
    weights <- gamma[, i]
    model <- lme4::lmer(formula, data = data, weights = weights)
  }
  snp <- NULL
  model
}

#' Performs the initial iteration of the FM-HighLD model
#' @param formula formula used for fitting the underlying linear model used by
#'  FM-HighLD
#' @param init_coef_mean a coefficients vector to init the algorithm with causal
#'  candidates selected from models with underlying coefficients sampled from a
#'  multivariate normal with mean `init_coef_mean`
#' @param data `data.frame` used for the FM-HighLD model
#' @param singletrait  a logical indicator determining if the model is for
#' multi_trait or single-trait fine-mapping
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @param fm_param a `FMParam` object with the parameters used to run `FMHighLD`
#' @return a list with the models, and probability matrix for the EM-algorithm
#' @importFrom Matrix Matrix
#' @importFrom stats predict
#' @importFrom purrr map_int map2
init_iteration <- function(formula, init_coef_mean, data, singletrait,
  ncausal_mixt, fm_param) {

  if (is.null(init_coef_mean)) {
    causal_candidates <- init_causal_candidates_random(data, singletrait,
      ncausal_mixt)
  } else {
    causal_candidates <- init_causal_candidates_coef(data, singletrait,
      formula, init_coef_mean, ncausal_mixt, fm_param)
  }

  n <- unique(purrr::map_int(causal_candidates, length))
  gamma_mat <- matrix(rep(1, (ncausal_mixt + 1) * n), nrow = n)
  gamma_mat <- Matrix::Matrix(gamma_mat)
  rownames(gamma_mat) <- names(causal_candidates)
  p <- ncol(gamma_mat)
  idxs <- seq_len(p)

  models <- purrr::map2(idxs[-1], causal_candidates,
    compute_ith_model,
    formula, data, gamma_mat, singletrait)

  FMIter(nassoc = nrow(data), singletrait = singletrait, models = models,
    causal_candidates = causal_candidates, gamma = gamma_mat, mu = NULL,
    sigma = NULL)

}

#' Initializes the causal candidates before fitting the model by picking the
#' SNPs based on a model with sampled random coefficients from a normal
#' distribution centered around the `init_coef_mean`
#' @param data `data.frame` used for the FM-HighLD model
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @param formula the formula of the underlying linear model
#' @param init_coef_mean the mean vector to sample the init. coefficients from
#'  a normal distribution
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @param fm_param a `FMParam` object with the parameters used to run `FMHighLD`
#' @return a list of length `ncausal_mixt` with a vector of one causal
#'  candidate per LD group (if `singletrait = TRUE`) or per combination
#'  of ld_group
#'  and trait if (`singletrait = FALSE`)
#' @importFrom rlang syms set_names
#' @importFrom dplyr group_by sample_n
#' @importFrom purrr map map2
#' @importFrom stringr str_c
#' @importFrom MASS mvrnorm
#' @export
#' @examples
#' nsnps <- 10
#' z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
#' annot <- matrix(rnorm(nsnps), ncol = 1)
#' rownames(annot) <- names(z)
#' colnames(annot) <- "annot"
#' ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
#' names(ld_cluster) <- names(z)
#' data <- build_fm_tibble(z, annot, ld_cluster, TRUE)
#' init_causal_candidates_coef(data, TRUE, response ~ 0 + annot, .5, 1,
#'  FMParam())
#' init_causal_candidates_coef(data, TRUE, response ~ 0 + annot, c(0,.5), 1,
#'  FMParam())
#' init_causal_candidates_coef(data, TRUE, response ~ 1 + annot, c(.3,.5), 1,
#'  FMParam())
init_causal_candidates_coef <- function(data, singletrait, formula,
  init_coef_mean, ncausal_mixt, fm_param) {

  if (singletrait) {
    init_model <- stats::lm(formula, data = data)
  } else {
    init_model <- lme4::lmer(formula, data = data)
  }

  if (! has_intercept(formula) & length(init_coef_mean) > 1 &
    init_coef_mean[1] != 0) {
    stop("the propossed intecept is different from zero")
  }

  init_coef <- stats::coef(init_model)
  if (!has_intercept(formula) & length(init_coef_mean) > 1) {
    init_coef_mean <- init_coef_mean[-1]
  }
  if (length(init_coef) != length(init_coef_mean)) {
    stop(stringr::str_c("The dimension of `init_coef_mean` is not the",
      "same than the parameters suggested by the formula", sep = " "))
  }
  if (is.null(names(init_coef_mean))) {
    warning("`init_coef_mean` is not named, will use names from formula")
    names(init_coef_mean) <- names(init_coef)
  }

  init_models <- replicate(ncausal_mixt, {
    if (has_intercept(formula)) {
      init_model$coefficients <- MASS::mvrnorm(1, mu = init_coef_mean,
        Sigma = diag(rep(1e-2, length(init_coef_mean))))
    } else {
      init_model$coefficients <- rnorm(1, init_coef_mean, sd = 1e-2)
      names(init_model$coefficients) <- names(init_coef_mean)
    }
    init_model
  }, simplify = FALSE)

  causal_candidate_list <- purrr::map(init_models,
    ~ select_causals_single(data, ., fm_param,
      "ld_cluster", names(init_coef)))
  causal_candidate_snps <- purrr::map(causal_candidate_list, "which_snp")
  causal_candidate_ld <- purrr::map(causal_candidate_list, "ld_cluster")
  causal_candidate_snps <- purrr::map2(
    causal_candidate_snps, causal_candidate_ld, rlang::set_names)
  causal_candidate_snps
}


#' Initializes the causal candidates before fitting the model by sampling a
#' random variant from each LD cluster
#' @param data `data.frame` used for the FM-HighLD model
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @return a list of length `ncausal_mixt` with a vector of one causal
#'  candidate per LD group (if `singletrait = TRUE`) or per combination of
#'  ld_group and trait if (`singletrait = FALSE`)
#' @importFrom rlang syms set_names
#' @importFrom dplyr group_by sample_n
#' @importFrom purrr map map2
#' @importFrom stringr str_c
#' @export
#' @examples
#' nsnps <- 10
#' z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
#' annot <- matrix(rnorm(nsnps), ncol = 1)
#' rownames(annot) <- names(z)
#' colnames(annot) <- "annot"
#' ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
#' names(ld_cluster) <- names(z)
#' data <- build_fm_tibble(z, annot, ld_cluster, TRUE)
#' init_causal_candidates_random(data, TRUE, 1)
#' init_causal_candidates_random(data, TRUE, 2)
#' z_list <- list()
#' z_list[["a"]] <- z[1:3]
#' z_list[["b"]] <- z[3:5]
#' z_list[["c"]] <- z[6:10]
#' data = build_fm_tibble(z_list, annot, ld_cluster, FALSE)
#' init_causal_candidates_random(data, FALSE, 3)
init_causal_candidates_random <- function(data, singletrait, ncausal_mixt) {

  if (singletrait) {
    data <- dplyr::mutate(data, sample_var = ld_cluster)
    ld_vars <- "ld_cluster"
  } else {
    data <- dplyr::mutate(data,
      sample_var = stringr::str_c(trait, ld_cluster, sep = ":"))
    ld_vars <- c("ld_cluster", "trait")
  }

  stopifnot(all(ld_vars %in% names(data)))

  group_data <- dplyr::group_by(data, !!!rlang::syms(ld_vars))
  group_data <- replicate(ncausal_mixt, {
    dplyr::sample_n(group_data, 1)
  }, simplify = FALSE)

  causal_cand <- purrr::map(group_data, "snp")
  samp_names <- purrr::map(group_data, "sample_var")
  purrr::map2(causal_cand, samp_names, rlang::set_names)

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

  collapse_sign <- " + "

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
