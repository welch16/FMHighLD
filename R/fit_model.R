#' Builds a `tibble::tibble` to be used by `FMHighLD`
#' @param response a vector for the `singletrait` case or a named list for
#'  the multitrait case with an element for each variant. In the multitrait
#'  case, each element of the list is a named vector where each element
#'  is named after the trait for which the association was tested.
#' @param annot_matrix a matrix with same number of rows as variants and number
#'  of columns as annotations
#' @param ld_clusters a character vector with the LD cluster to which each
#'  variant belongs
#' @param singletrait  a logical indicator determining if the model is for
#' multi_trait or single-trait fine-mapping
#' @return a `tibble::tibble` with at least columns SNP, response, ld_cluster
#'  and the column names of `annot_matrix`. If `singletrait == FALSE`, it will
#'  also have a `trait column`.
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map2 partial reduce
#' @importFrom dplyr mutate bind_rows inner_join
#' @examples
#' nsnps <- 10
#' z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
#' annot <- matrix(rnorm(nsnps), ncol = 1)
#' rownames(annot) <- names(z)
#' colnames(annot) <- "annot"
#' ld_cluster <- sample(stringr::str_c("ld", seq_len(5)), nsnps, replace = TRUE)
#' names(ld_cluster) <- names(z)
#' build_fm_tibble(z, annot, ld_cluster, TRUE)
#' z_list <- list()
#' z_list[["a"]] <- z[1:3]
#' z_list[["b"]] <- z[3:5]
#' z_list[["c"]] <- z[6:10]
#' build_fm_tibble(z_list, annot, ld_cluster, FALSE)
build_fm_tibble <- function(response, annot_matrix, ld_clusters, singletrait) {

  fm_tibble <- tibble::tibble(snp = get_snp_names(response, singletrait))
  if (singletrait) {
    response_tibble <- tibble::tibble(snp = names(response), response)
  } else {
    response_tibble <- purrr::map(response,
      ~ tibble::tibble(snp = names(.), response = .))
    response_tibble <- purrr::map2(response_tibble, names(response_tibble),
      ~ dplyr::mutate(.x, trait = .y))
    response_tibble <- dplyr::bind_rows(response_tibble)
  }

  annot_tibble <- as.data.frame(annot_matrix)
  annot_tibble <- tibble::as_tibble(annot_tibble, rownames = "snp")
  ld_clust_tibble <- tibble::tibble(snp = names(ld_clusters),
    ld_cluster = ld_clusters)

  purrr::reduce(list(fm_tibble, response_tibble, annot_tibble, ld_clust_tibble),
    purrr::partial(dplyr::inner_join, by = "snp"))
}

#' Gets a sorted vector with the SNP names
#' @param response Either a numeric vector when `singletrait == TRUE` or a
#'   list with named vectors when `singletrait == FALSE`
#' @param singletrait  a logical indicator determining if the model is for
#' multi_trait or single-trait fine-mapping
#' @return a sorted character vector with the names of all the SNPs in the
#'   response
#' @importFrom purrr map
#' @export
#' @examples
#' nsnps <- 10
#' z <- rlang::set_names(rnorm(nsnps), stringr::str_c("snp", seq_len(nsnps)))
#' get_snp_names(z, TRUE)
#' z_list <- list()
#' z_list[["a"]] <- z[1:3]
#' z_list[["b"]] <- z[3:5]
#' z_list[["c"]] <- z[6:10]
#' get_snp_names(z_list,  FALSE)
get_snp_names <- function(response, singletrait) {

  if (singletrait) {
    snp_names <- names(response)
  } else {
    snp_names <- purrr::map(response, names)
    snp_names <- unique(do.call(c, snp_names))
  }
  sort(snp_names)
}

#' Compute the mixture probabilities from a probability matrix
#'
#' @param prob_matrix a matrix of dimensions nassoc x nmixtures
#' @return a vector of length nmixtures
#' @export
#' @examples
#' compute_mixture_prob(matrix(rep(1, 100), ncol = 2))
compute_mixture_prob <- function(prob_matrix) {

  prob_matrix <- as.matrix(prob_matrix)
  causal_prob <- colMeans(prob_matrix)
  causal_prob / sum(causal_prob)
}

#' Fit the FMHighLD model
#' @param response a vector for the `singletrait` case or a named list for
#'  the multitrait case with an element for each variant. In the multitrait
#'  case, each element of the list is a named vector where each element
#'  is named after the trait for which the association was tested.
#' @param annot_matrix a matrix with same number of rows as variants and number
#'  of columns as annotations
#' @param ld_clusters a character vector with the LD cluster to which each
#'  variant belongs
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @param formula a `formula` for the underlying linear models, by default is
#'  `NULL`, in which case `fmhighld_fit` will create the formula.
#' @param init_coef_mean a coefficients vector to init the algorithm with causal
#'  candidates selected from models with underlying coefficients sampled from a
#'  multivariate normal with mean `init_coef_mean`
#' @param fm_param a `FMParam` object with the parameters used to run `FMHighLD`
#' @param save_iter a logical indicator determining whether the iteration data
#'  is going to be returned
#' @param verbose a logical indicator determining whether messages are going to
#'  be used
#' @return results
#' @export
fmhighld_fit <- function(response, annot_matrix, ld_clusters,
  singletrait = TRUE, ncausal_mixt = 1, formula = NULL, init_coef_mean = NULL,
  fm_param = FMParam(), save_iter = FALSE,
  verbose = FALSE) {

  stopifnot(is.matrix(annot_matrix) | is.vector(annot_matrix))
  if (is.vector(annot_matrix)) {
    warning("will convert `annot_matrix` argument into a matrix")
    annot_matrix <- as.matrix(annot_matrix, ncol = 1)
  }

  same_names <- check_variant_names(response, annot_matrix, ld_clusters)

  if (! same_names) {
    warning("There variants are not named or the names don't match, will match
      'response', 'annot_matrix' and 'ld_clusters' by position")
    if (is.list(response)) {
      rlang::abort("`response` is a list, please check the variant names")
    }
    names(response) <- stringr::str_c("snp", seq_along(response))
    rownames(annot_matrix) <- names(response)
    names(ld_clusters) <- names(response)
  }

  snp_names <- get_snp_names(response, singletrait)
  annot_names <- colnames(annot_matrix)

  stopifnot(
    length(snp_names) == nrow(annot_matrix),
    length(snp_names) == length(ld_clusters))

  ## init data structure
  fmld_data <- build_fm_tibble(response, annot_matrix, ld_clusters, singletrait)

  # init algorithm parameteres
  iter <- 0
  error_bound <- error_bound(fm_param)
  max_iter <- max_iter(fm_param)
  min_tol <- min_tol(fm_param)
  annot_tol <- annot_tol(fm_param)

  if (is.null(formula)) {
    formula <- build_formula("response", colnames(annot_matrix), TRUE)
    warning("using formula ", as.character(formula))
  }

  init <- init_iteration(formula, init_coef_mean, fmld_data,
    singletrait, ncausal_mixt, fm_param)
  # causal_prob <- compute_mixture_prob(probmatrix(init))

  continue <- TRUE
  current_iter <- init
  if (save_iter) {
    like_vec <- matrix(rep(NA, max_iter * ncausal_mixt), ncol = ncausal_mixt)
    full_like_vec <- matrix(rep(NA, max_iter * ncausal_mixt), ncol = 
      ncausal_mixt)
    all_models <- vector(mode = "list", length = max_iter)
  }

  aux_df <- as.data.frame(fmld_data)

  while (continue) {

    iter <- iter + 1
    prev_iter <- current_iter
    model_list <- models(prev_iter)
    sigma0 <- sigma0(prev_iter, fmld_data[[get_response_name(formula)]])
    if (verbose) {
      message("starting iter ", iter)
    }

    if (singletrait) {
      causal_list <- purrr::map(model_list,
        ~ select_causals_single(fmld_data, .x, fm_param, "ld_cluster",
          annot_names))

    } else {
      browser()
      debugonce(select_causals_multi)
      causal_list <- purrr::map(model_list,
        ~ select_causals_multi(train_data, .x, response, to_group,
          fm_param, cond_res))
    }



    if (singletrait) {
      current_iter <- em_iteration_single(formula, fmld_data, causal_list,
        prev_iter, sigma0, fm_param, verbose)
    } else {
      current_iter <- em_iteration_multi(formula, fmld_data, causal_list,
        prev_iter, sigma0, fm_param, verbose)
    }

    ## compute metrics
    entropy_vec <- prob_metric(current_iter, prev_iter)
    mccl_vec <- mccl(current_iter, prev_iter)
    coef_diff <- coef_diff(current_iter, prev_iter, FALSE, FALSE)
    pl <- philips(c(entropy_vec[2], mccl_vec, coef_diff))

    curr_data <- get_causal_data(current_iter, fmld_data)
    curr_like <- purrr::map_dbl(curr_data, ~ loglikelihood(current_iter, .))

    prev_data <- get_causal_data(prev_iter, fmld_data)
    if (iter == 1) {
      prev_like <- rep(-Inf, ncausal_mixt)
    } else {
      prev_like <- purrr::map_dbl(prev_data, ~ loglikelihood(prev_iter, .))
    }
    if (all(prev_like > curr_like)) {
      current_iter <- prev_iter
    }
    if (save_iter) {
      like_vec[iter, ] <- curr_like
      full_like_vec[iter, ] <- prev_like
      all_models[iter] <- current_iter
    }

    print(iter)
    print(stringr::str_c("pl: ", pl))
    print(stringr::str_c("loglikeli: ", max(curr_like)))

    continue <- iter < max_iter & pl >= min_tol

  }

  # remember to add final causal candidates to object
  if (save_iter) {
    out <- list(all_models = all_models, loglike = like_vec,
      full_loglike = full_like_vec, final = current_iter)
  } else {
    out <- current_iter
  }
  out
}

#' Fit the FMHighLD model using flexmix
#' @param response a vector for the `singletrait` case or a named list for
#'  the multitrait case with an element for each variant. In the multitrait
#'  case, each element of the list is a named vector where each element
#'  is named after the trait for which the association was tested.
#' @param annot_matrix a matrix with same number of rows as variants and number
#'  of columns as annotations
#' @param ld_clusters a character vector with the LD cluster to which each
#'  variant belongs
#' @param singletrait  a logical indicator determining if the model is for
#'  multi_trait or single-trait fine-mapping
#' @param formula a `formula` for the underlying linear models, by default is
#'  `NULL`, in which case `fmhighld_fit` will create the formula.
#' @param fm_param a `FMParam` object with the parameters used to run `FMHighLD`
#' @param save_iter a logical indicator determining whether the iteration data
#'  is going to be returned
#' @param verbose a logical indicator determining whether messages are going to
#'  be used
#' @return results
#' @export
fmhighld_fit_em <- function(response, annot_matrix, ld_clusters,
  singletrait = TRUE, formula = NULL,
  fm_param = FMParam(), save_iter = FALSE,
  verbose = FALSE) {

  ncausal_mixt <- 1

  stopifnot(is.matrix(annot_matrix) | is.vector(annot_matrix))
  if (is.vector(annot_matrix)) {
    warning("will convert `annot_matrix` argument into a matrix")
    annot_matrix <- as.matrix(annot_matrix, ncol = 1)
  } else if (is.matrix(annot_matrix) & ncol(annot_matrix) > 1) {
    # check if the matrix contains the intercept vector
    int_idx <- apply(annot_matrix, 2, function(x)all(x == 1))
    cnms <- colnames(annot_matrix)
    annot_matrix <- as.matrix(annot_matrix[, 2], ncol = 1)
    colnames(annot_matrix) <- cnms[!int_idx]
  }

  same_names <- check_variant_names(response, annot_matrix, ld_clusters)

  if (! same_names) {
    warning("There variants are not named or the names don't match, will match
      'response', 'annot_matrix' and 'ld_clusters' by position")
    if (is.list(response)) {
      rlang::abort("`response` is a list, please check the variant names")
    }
    names(response) <- stringr::str_c("snp", seq_along(response))
    rownames(annot_matrix) <- names(response)
    names(ld_clusters) <- names(response)
  }

  snp_names <- get_snp_names(response, singletrait)
  annot_names <- colnames(annot_matrix)

  stopifnot(
    length(snp_names) == nrow(annot_matrix),
    length(snp_names) == length(ld_clusters))

  ## init data structure
  fmld_data <- build_fm_tibble(response, annot_matrix, ld_clusters, singletrait)

  # init algorithm parameters
  iter <- 0
  error_bound <- error_bound(fm_param)
  max_iter <- max_iter(fm_param)
  min_tol <- min_tol(fm_param)
  annot_tol <- annot_tol(fm_param)

  if (is.null(formula)) {
    formula <- build_formula("response", colnames(annot_matrix), TRUE)
    warning("using formula ", as.character(formula))
  }

  init <- init_iteration_em(formula, fmld_data, singletrait, ncausal_mixt,
    fm_param)

  continue <- TRUE
  current_iter <- init
  if (save_iter) {
    like_vec <- rep(NA, max_iter)
    full_like_vec <- rep(NA, max_iter)
    all_models <- vector(mode = "list", length = max_iter)
  }

  aux_df <- as.data.frame(fmld_data)


  while (continue) {
    iter <- iter + 1
    prev_iter <- current_iter
    model_list <- models(prev_iter)
    sigma0 <- sigma0(prev_iter, fmld_data[[get_response_name(formula)]])
    if (verbose) {
      message("starting iter ", iter)
    }

    if (singletrait) {
      causal_list <- purrr::map(model_list,
        ~ select_causals_single(fmld_data, .x, fm_param, "ld_cluster",
          annot_names))
    } else {
      browser()
      debugonce(select_causals_multi)
      causal_list <- purrr::map(model_list,
        ~ select_causals_multi(train_data, .x, response, to_group,
          fm_param, cond_res))
    }


    if (singletrait) {
      current_iter <- fm_iteration_single(formula, fmld_data, causal_list,
        prev_iter, sigma0, fm_param, ncausal_mixt, verbose)
    } else {
      current_iter <- em_iteration_multi(formula, fmld_data, causal_list,
        prev_iter, sigma0, fm_param, verbose)
    }

    ## compute metrics
    entropy_vec <- prob_metric(current_iter, prev_iter)[-1]
    mccl_vec <- mccl(current_iter, prev_iter)
    coef_diff <- coef_diff_em(current_iter, prev_iter, FALSE, FALSE)
    pl <- philips(c(entropy_vec, mccl_vec, coef_diff))

    curr_like <- purrr::map_dbl(models(current_iter), flexmix::logLik)
    prev_like <- purrr::map_dbl(models(prev_iter), flexmix::logLik)
    if (prev_like > curr_like) {
      current_iter <- prev_iter
    }
    if (save_iter) {
      like_vec[iter] <- curr_like
      full_like_vec[iter] <- purrr::map_dbl(models(current_iter),
        flexmix::logLik, newdata = aux_df)
      all_models[iter] <- models(current_iter)[[1]]
    }

    # tibble::tibble(like = like_vec, id = seq_along(like_vec)) %>%
    #   ggplot(aes(id, like_vec)) +
    #     geom_line()
    # all_models[[which.max(like_vec)]] %>%
    #   flexmix::parameters()

    if (verbose) {
      print(iter)
      print(purrr::map(models(current_iter), flexmix::parameters)[[1]])
      print(stringr::str_c("pl: ", pl))
      print(stringr::str_c("loglikeli: ", curr_like))
    }

    # browser()

    continue <- iter < max_iter & pl >= min_tol

    # we have causal candidates already, now we need to apply the rest of the
    # EM-algorithm

  }

  if (save_iter) {
    out <- list(all_models = all_models, loglike = like_vec,
      full_loglike = full_like_vec, final = current_iter)
  } else {
    out <- current_iter
  }
  out
}


#' Performs the initial iteration of the FM-HighLD model
#' @param formula formula used for fitting the underlying linear model used by
#'  FM-HighLD
#' @param data `data.frame` used for the FM-HighLD model
#' @param singletrait  a logical indicator determining if the model is for
#' multi_trait or single-trait fine-mapping
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @param fm_param a `FMParam` object with the parameters used to run `FMHighLD`
#' @return a list with the models, and probability matrix for the EM-algorithm
#' @importFrom Matrix Matrix
#' @importFrom flexmix predict
#' @importFrom stats predict
#' @importFrom purrr map_int map2
init_iteration_em <- function(formula, data, singletrait, ncausal_mixt,
  fm_param) {

  causal_candidates <- init_causal_candidates_random(data, singletrait,
      ncausal_mixt)

  n <- unique(purrr::map_int(causal_candidates, length))
  gamma_mat <- matrix(rep(1, (ncausal_mixt + 1) * n), nrow = n)
  gamma_mat <- Matrix::Matrix(gamma_mat)
  rownames(gamma_mat) <- names(causal_candidates)

  flexmix_model <- flexmix_fit(formula, data, causal_candidates,
    ncausal_mixt + 1)
  flexmix_model <- flexmix::relabel(flexmix_model)
  models <- list(flexmix_model)

  FMIter(nassoc = nrow(data), singletrait = singletrait, models = models,
    causal_candidates = causal_candidates, gamma = gamma_mat, mu = NULL,
    sigma = NULL)

}


#' Fits the EM algorithm part of the FMHighLD model using flexmix
#' @param formula Formula of the non-zero mean part of the FMHighLD model
#' @param data `data.frame` used for the FM-HighLD model
#' @param causal_candidates a list with vectors of causal candidates
#' @param ncausal_mixt number of components in the mixture model
#' @param ... other parameters using to fit `flexmix::flexmix`
#' @return a `flexmix::flexmix` object
#' @importFrom dplyr filter
#' @importFrom stringr str_c
#' @importFrom flexmix flexmix FLXMRglm
flexmix_fit <- function(formula, data, causal_candidates,
  ncausal_mixt, ...) {

  fit_data <- dplyr::filter(data, snp %in% causal_candidates[[1]])
  char_formula <- as.character(formula)
  mixt_formula <- as.formula(stringr::str_c("~", char_formula[3]))

  fmodel1 <- flexmix::FLXMRglm(mixt_formula, family = "gaussian")
  fmodel2 <- flexmix::FLXMRglm(family = "gaussian")

  snp <- NULL
  flexmix::flexmix(response ~ 1, k = ncausal_mixt, data = fit_data,
    model = list(fmodel1, fmodel2), ...)
}

#' Performs an Expectation-Maximization iterations for the single
#' trait model
#' @param formula a formula object with the underlying linear model
#' @param data a `data.frame` object with the annotation matrix and z-values
#' @param causal_current A list of `tibble::tibble` with columns `ld_cluster`
#'  and `which_snp` obtained in the current iteration of the algorithm
#' @param prev_iter a `FMIter` object obtained in the previous iteration of the
#'  FMHighLD algorithm
#' @param sigma0 the error variance estimate of the background model
#' @param fm_param configuration parameter used to pick the causal candidate
#'  per group
#' @param ncausal_mixt the number of mixtures used in the causal models
#' @param verbose a logical indicator of adding process messages
#' @return a list with the updated parameters in one EM iteration
#' @importFrom methods new
#' @importFrom purrr map map2 map_dbl map_int
#' @importFrom stats predict update model.matrix
#' @importFrom broom glance
fm_iteration_single <- function(
  formula, data, causal_current, prev_iter, sigma0, fm_param, ncausal_mixt,
   verbose) {

  if (verbose) message("--starting EM iter")

  error_bound <- error_bound(fm_param)

  # convert to probabilities if they are not
  # pi <- compute_mixture_prob(probmatrix(prev_iter))
  model_list <- models(prev_iter)

  # generate mean matrix
  if (verbose) message("--calculating means")
  new_data <- purrr::map(causal_current, dplyr::select, which_snp)
  new_data <- purrr::map(new_data, dplyr::left_join, data,
    by = c(which_snp = "snp"))
  model_list <- purrr::map(model_list, flexmix::relabel)

  predictions <- purrr::map2(model_list, new_data, flexmix::predict)
  predictions <- predictions[[1]]
  params <- flexmix::parameters(model_list[[1]])
  # print(params)
  mu <- cbind(predictions[[2]][, 2], predictions[[1]][, 1])

  # generate variance matrix
  if (verbose) message("--calculating variances")
  sigma <- cbind(rep(params[[2]]["sigma", 2], nrow(mu)),
    rep(params[[1]]["sigma", 1], nrow(mu)))

  if (verbose) message("--calculating weights")
  gamma_mat <- flexmix::posterior(model_list[[1]],
    newdata = as.data.frame(new_data))
  pi <- colSums(gamma_mat)
  pi <- pi / sum(pi)


  if (verbose) message("--fitting new models")
  candidates <- purrr::map(causal_current, "which_snp")
  new_data <- purrr::map(new_data, dplyr::rename, snp = which_snp)
  flexmix_model <- flexmix_fit(formula, as.data.frame(new_data[[1]]),
    candidates, ncausal_mixt + 1)
  flexmix_model <- flexmix::relabel(flexmix_model)
  models <- list(flexmix_model)

  which_snp <- NULL

  FMIter(nassoc = nrow(new_data), singletrait = TRUE, models = models,
    causal_candidates = extract_causal_vector(causal_current),
    gamma = Matrix::Matrix(gamma_mat), mu = Matrix::Matrix(mu),
    sigma = Matrix::Matrix(sigma))
}
