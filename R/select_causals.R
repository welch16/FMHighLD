
#' Checks whether the intercept is include in the underlying linear model
#' @param features names of the estimated coefficient
#' @return a logical indicating whether include the intercept in the model
#' @importFrom stringr str_detect
include_intercept <- function(features) {
  any(stringr::str_detect("Intercept", features))
}


#' Builds the design matrix
#' @param data a `data.frame` used to fit the FM-HighLD model
#' @param beta the estimated linear model coefficient. In the multi-trait case,
#'  this is the estimated coefficient for the fixed effects
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom Matrix Matrix
#' @importFrom stringr str_detect
#' @export
#' @examples
#' set.seed(1234)
#' nsnps <- 20
#' data <- data.frame(var1 = rnorm(nsnps), var2 = rnorm(nsnps))
#' beta <- c("(Intercept)" = 2, "var1" = 1, "var2" = 2)
#' build_design_matrix(data, beta)
build_design_matrix <- function(data, beta) {

  features <- names(beta)
  include_intercept <- include_intercept(features)
  if (include_intercept) features <- features[-1]

  design_matrix <- dplyr::select(data, tidyselect::all_of(features))
  design_matrix <- Matrix::Matrix(as.matrix(design_matrix))

  if (include_intercept) {
    design_matrix <- cbind(1, design_matrix)
  }
  design_matrix
}

#' Selects the causal variants in the locus for a single-trait model
#' @param data a `data.frame` used to fit the FM-HighLD model
#' @param model a fitted model with the previously selected causal candidates
#' @param group name of the grouping variable used to pick the causal
#'  candidates
#' @param rand_param configuration parameter used to pick the causal candidate
#'  per group
#' @importFrom stats formula terms coef model.matrix
select_causals_single <- function(data, model, group, rand_param) {

  model_formula <- stats::formula(model)
  response <- as.character(model_formula)[2]
  beta <- stats::coef(model)
  design_matrix <- build_design_matrix(data, beta)
  fitted <- as.numeric(design_matrix %*% beta)
  residuals <- data[[response]] - fitted

  causal_rule(data, residuals, rand_param, group)
}

#' Selects the causal variants in the locus for a multi-trait model
#' @param data a `data.frame` used to fit the FM-HighLD model
#' @param model a fitted model with the previously selected causal candidates
#' @param group name of the grouping variable used to pick the causal
#'  candidates
#' @param rand_param configuration parameter used to pick the causal candidate
#'  per group
#' @param cond_res a logical indicating whether to use conditional residuals
#'  of the linear mixed model
#' @importFrom stats formula
#' @importFrom lme4 fixef
select_causals_multi <- function(data, model, group,
  rand_param, cond_res = FALSE) {

  model_formula <- stats::formula(model)
  response <- as.character(model_formula)[2]

  beta <- lme4::fixef(model)
  design_matrix <- build_design_matrix(data, beta)
  fitted <- as.numeric(design_matrix %*% beta)
  residuals <- data[[response]] - fitted

  if (cond_res) {
    stop("need to implement conditional residuals")

  }

  causal_rule(data, residuals, rand_param, group)
}


#' Applies the causal selection rule to the residual vector
#' @param data a `data.frame` used to fit the FM-HighLD model
#' @param residuals a vector of residuals for all the SNPs
#' @param group name of the grouping variable used to pick the causal
#'  candidates
#' @param rand_param configuration parameter used to pick the causal candidate
#'  per group
#' @importFrom dplyr select group_by mutate summarize ungroup bind_rows
#' @importFrom dplyr top_n sample_n filter inner_join anti_join
#' @importFrom tidyselect one_of
#' @importFrom nnet which.is.max
#' 
causal_rule <- function(data, residuals, rand_param, to_group) {

  new_data <- dplyr::select(SNP, tidyselect::one_of(to_group))
  new_data <- dplyr::mutate(new_data, res_seq = residuals ^ 2)
  new_data <- dplyr::group_by(new_data, to_group)

  if (is.null(rand_param)) {
    out <- dplyr::summarize(new_data,
      which_snp = SNP[nnet::which.is.max(-res_seq)], .groups = "drop")
  } else if (rand_param$strat == "all") {

    out <- dplyr::summarize(new_data,
      which_snp =
        SNP[select_kth_random(res_seq, rand_param$prob, rand_param$select)],
        .groups = "drop")

  } else if (rand_param$strat == "pick_M"){

    out <- dplyr::summarize(new_data,
      n = n(),
      which_snp = SNP[nnet::which.is.max(-res_seq)])

    if (rand_param$which == "any"){

      picks <- dplyr::ungroup(out)
      picks <- dplyr::filter(picks, n > 1)
      picks <- dplyr::select(picks, -n, -which_snp)
      picks <- dplyr::sample_n(picks, rand_param$M)

    } else if (rand_param$which == "largeLD"){

      picks <- dplyr::ungroup(out)
      picks <- dplyr::top_n(picks, rand_param$M, wt = n)
      picks <- dplyr::select(picks, -n, -which_snp)

    }

    pick_data <- dplyr::inner_join(picks, new_data, by = to_group)
    pick_data <- dplyr::group_by(pick_data, !! rlang::syms(to_group))
    pick_data <- dplyr::summarize(pick_data,
      which_snp =
        SNP[select_kth_random(res_seq, rand_param$prob, rand_param$select)],
      .groups = "drop")

    out <- dplyr::ungroup(out)
    out <- dplyr::select(out, -n)
    out <- dplyr::anti_join(out, picks, by = to_group)
    out <- dplyr::bind_rows(out, pick_data)

  }

  n <- which_snp <- SNP <- res_seq <- NULL

  out
}

#' Selects the k-th smallest squared residual randomly with probability
#'  prob
#' @param res a vector of residuals
#' @param prob a probability to randomize the selection of the k-th
#'  smallest residual
#' @param k an integer value
#' @return the id of the selected k-th smallest residual
#' @importFrom dplyr if_else
#' @importFrom stats runif
select_kth_random <- function(res, prob = 0.5, k = 2) {

  nsnps <- length(res)
  if (nsnps == 1) {
    1
  } else {
    nsnps <- length(res)
    dplyr::if_else(stats::runif(1) <= prob,
      which(order(res) == min(k, nsnps)),
      nnet::which.is.max(-res))
  }
}
