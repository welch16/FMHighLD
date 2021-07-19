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
