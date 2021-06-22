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

  # if (save_iter) {
  #   states <- list()
  # } else {
  #   states <- NULL
  # }

  if (is.null(formula)) {
    formula <- build_formula("response", colnames(annot_matrix), TRUE)
    warning("using formula ", as.character(formula))
  }

  init <- init_iteration(formula, init_coef_mean, fmld_data,
    singletrait, ncausal_mixt, fm_param)
  # causal_prob <- compute_mixture_prob(probmatrix(init))

  continue <- TRUE
  current_iter <- init
  

  while (continue) {
    # browser()
    iter <- iter + 1
    # if (iter > 10) browser()
    prev_iter <- current_iter
    model_list <- models(prev_iter)
    sigma0 <- sigma0(prev_iter, fmld_data[[get_response_name(formula)]])
    if (verbose) {
      message("starting iter ", iter)
    }

# em_iteration_single <- function(
#   formula, data, pi, model_list, sigma0, fm_param, verbose) {
# em_iteration_multi <-
#   function(formula, data, pi, model_list, sigma0, fm_param, verbose = TRUE) {


    if (singletrait) {
      causal_list <- purrr::map(model_list,
        ~ select_causals_single(fmld_data, .x, fm_param, "ld_cluster",
          annot_names))
      causal_wide <- causal_list %>%
        purrr::map2(
          stringr::str_c("which_snp", seq_len(ncausal_mixt), sep = "_"),
            ~ rlang::set_names(.x, c("ld_cluster", .y))) %>%
        purrr::reduce(purrr::partial(dplyr::inner_join, by = "ld_cluster"))
      # print(dplyr::pull(causal_wide, which_snp_1))

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

    curr_like <- likelihood(current_iter, fmld_data, 1 - 1e-3)
    print(iter)
    print(purrr::map(models(current_iter), coef))
    print(stringr::str_c("pl: ", pl))
    print(stringr::str_c("loglikeli: ", log(curr_like)))

    continue <- iter < max_iter & pl >= min_tol

        # mce = map2(causal_list,prev_causals,
        #            inner_join,by = to_group) %>%
        #     map(mutate, diff = which_snp != SNP) %>%
        #     map(pull,diff) %>%
        #     map_dbl(mean) %>%
        #     max()
       
        # weight_list = iter_list %>% map( ~ .$gamma[,2]) 
        
        # entropy = prev_causals %>%
        #     map(pull,prob) %>%
        #     map2_dbl(weight_list, ~ kl(.y,.x))
        # jsdist = prev_causals %>%
        #     map(pull,prob) %>%
        #     map2_dbl(weight_list,jsd)

        # metrics = c("js_dist" = jsdist)

        # continue = any(metrics >= tol) & iter < max_iter
        # iter = iter + 1
    # we have causal candidates already, now we need to apply the rest of the
    # EM-algorithm

  }

  # remember to add final causal candidates to object
  print(round(compute_mixture_prob(probmatrix(current_iter)), 2))
  current_iter
}

iterative_select_causal_init <- function(train_data,
                                      model_formula,
                                      causal_candidates,
                                      features = c("ATAC","TF_alle"),
                                      response = "eQTL_tStat",                            
                                      to_group = c("eQTL_gene","cluster"),
                                      cond_res = FALSE,
                                      rand_param = NULL,
                                      max_iter = 100,
                                      tol = 1e-6,save_state = FALSE,
                                      univariate = FALSE,
                                      verbose = FALSE)
{
    ## init iteration parameters
    error_bound = 1e-4
    continue = TRUE
    iter = 0
  
    states = NULL
    if(save_state){
        states = list()
    }

    init = init_iteration(model_formula,                   
                         causal_candidates,
                         response,univariate)

    model_list = init$models
    pi = init$gamma %>% colMeans()

    background_error = rep_along(model_list,1)
    
    prev_causals = list(
        causal_candidates %>%
        mutate(prob = 1/2)) 


    while(continue){
    
        if(verbose) message("Iter ", iter)
        
        if(univariate){
            causal_list = map(model_list,
                              ~ select_causals_single(train_data,.x,
                                                   response,to_group,
                                                   rand_param))
        }else{
        	
            causal_list = map(model_list,
                              ~ select_causals_multi(train_data,.x,
                                                response,to_group,rand_param,
                                                cond_res))
        }

        filtered_list = map(causal_list,
                            ~ rebuild_data(train_data,.,to_group))

        if(univariate){
            iter_list = map( filtered_list,
                            ~ em_iteration_single(formula = model_formula,
                                           data = .,
                                           response = response,
                                           pi = pi,
                                           model_list = model_list,
                                           background_error = background_error
                                           ))
            
            param_list  = map2(iter_list, filtered_list,
                               parse_estimated_parameter_single,response) %>%
                unlist(recursive = FALSE)

            param_metric = map2(
                param_list[!str_detect(names(param_list),"random_eff")],
                list(
                    pi = pi,
                    param_list = map(model_list,coef),
                    residual_error_list = map(model_list,broom::glance) %>%
                        map("sigma"),
                    background_error = background_error),
                ~ ( unlist(.x) - unlist(.y))^2) %>%
                map_dbl(sum) %>% { sum(.) / (1 + max(.))} %>% sqrt()

            pi = param_list$pi
            model_list = iter_list %>% map("models") %>%
                map(pluck,1)
            background_error = param_list$background_error 
            
        }else{
            iter_list = map( filtered_list,
                            ~ em_iteration_multi(formula = model_formula,
                                           data = .,
                                           response = response,
                                           pi = pi,
                                           model_list = model_list,
                                           background_error = background_error,
                                           REML = TRUE,
                                           verbose = verbose))
            param_list  = map2(iter_list, filtered_list,
                               parse_estimated_parameter_multi,response) %>%
                unlist(recursive = FALSE)

            error_list = map(model_list,
                                  broom::tidy) %>%
                map(filter,group != "fixed")
            
            param_metric = map2(
                param_list[!str_detect(names(param_list),"random_eff")],
                list(
                    pi = pi,
                    param_list = map(model_list,fixef),
                    gene_error_list = map(error_list,
                                          filter,group != "Residual") %>%
                        map(pull,estimate),
                    residual_error_list = map(error_list,filter,group == "Residual") %>%
                        map(pull,estimate),
                    background_error = background_error),
                ~ ( unlist(.x) - unlist(.y))^2) %>%
                map_dbl(sum) %>% { sum(.) / (1 + max(.))} %>% sqrt()

            pi = param_list$pi
            model_list = iter_list %>%
                map(pluck,"models")
            background_error = param_list$background_error 
            
            
            model_list %<>% unlist()
        }            
      
        mce = map2(causal_list,prev_causals,
                   inner_join,by = to_group) %>%
            map(mutate, diff = which_snp != SNP) %>%
            map(pull,diff) %>%
            map_dbl(mean) %>%
            max()
       
        weight_list = iter_list %>% map( ~ .$gamma[,2]) 
        
        entropy = prev_causals %>%
            map(pull,prob) %>%
            map2_dbl(weight_list, ~ kl(.y,.x))
        jsdist = prev_causals %>%
            map(pull,prob) %>%
            map2_dbl(weight_list,jsd)

        metrics = c("js_dist" = jsdist)

        continue = any(metrics >= tol) & iter < max_iter
        iter = iter + 1
                  
        prev_causals = causal_list %>%
            map(dplyr::rename,SNP = which_snp) %>%
            map2(iter_list,
                 ~ mutate(.x,prob = .y$gamma[,2]))

        metrics["miss_class_err"] = mce
        metrics["philips"] = param_metric        

        if(univariate){

            model_metrics = model_list %>%
                map(broom::glance)
                     
        }else{

            model_metrics = model_list %>%
                map(broom::glance)
            
        }

        metrics["filtered_err"] = model_metrics %>% map_dbl("sigma")
        metrics["filtered_logLik"] = model_metrics %>% map_dbl("logLik")
        metrics["filtered_AIC"] = model_metrics %>% map_dbl("AIC")
        metrics["filtered_BIC"] = model_metrics %>% map_dbl("BIC")
        metrics["kl"] = entropy

        if(save_state){
            ## states[[iter]] = param_list
            states[[iter]] = list(
                "param" = param_list,
                "causal" = prev_causals,
                "metrics" = metrics)            
            
        }

        

    }


    metrics["iter_exit"] = iter - 1

    if(univariate){
        final_causals = map(model_list,
                             ~ select_causals_single(train_data,.x,
                                                  response,to_group,
                                                  rand_param))
        ranef_list = NULL
    }else{
        final_causals = map( model_list,
                          ~ select_causals_multi(train_data,.x,
                                                 response,to_group,rand_param,
                                                 cond_res))
        ranef_list = map(iter_list,"models") %>%
            unlist(recursive =FALSE) %>%
            map(ranef) %>%
            map(as_tibble) %>%
            map(spread,term,condval)
    }
    
     
    filtered_data = map(final_causals,
                        ~ rebuild_data(train_data,.,to_group))


    
    list(
        estimated_parameters = map(param_list,unlist,recursive = FALSE),
        metrics = metrics,
        final_causal = final_causals[[1]],
        gamma = iter_list[[1]]$gamma,
        states = states,
        ranef = ranef_list)
    
  
}
