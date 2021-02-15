
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
#' multi_trait or single-trait fine-mapping
#' @return a list with the models, and probability matrix for the EM-algorithm
#' @importFrom Matrix Matrix
#' @importFrom stats predict
init_iteration <- function(formula, data, singletrait = TRUE) {

  n <- nrow(data)
  gamma_mat <- do.call(rbind, lapply(seq_len(n), function(x)c(1, 1)))
  gamma_mat <- Matrix::Matrix(gamma_mat)

  p <- ncol(gamma_mat)
  idxs <- seq_len(p)

  models <- lapply(idxs[-1],
    compute_ith_model,
    formula, data, gamma_mat, singletrait)

  FMIter(nassoc = nrow(data), singletrait = singletrait, models = models,
    gamma = gamma_mat, mu = NULL, sigma = NULL)

}


#' Initializes the causal candidates before fitting the model
#' @param data `data.frame` used for the FM-HighLD model
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @return a list of length `ncausal_mixt` with a vector of one causal candida-
#'  te per ld group (if `singletrait = TRUE`) or per combination of ld_group
#'  and trait if (`singletrait = FALSE`)
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
#' init_causal_candidates(data, TRUE, 1)
#' init_causal_candidates(data, TRUE, 2)
#' z_list <- list()
#' z_list[["a"]] <- z[1:3]
#' z_list[["b"]] <- z[3:5]
#' z_list[["c"]] <- z[6:10]
#' data = build_fm_tibble(z_list, annot, ld_cluster, FALSE)
#' init_causal_candidates(data, FALSE, 3)
init_causal_candidates <- function(data, singletrait, ncausal_mixt) {

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
