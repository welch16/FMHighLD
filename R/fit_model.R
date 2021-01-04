

# set.seed(12345)


#' Fit the FMHighLD model
#' @param formula formula used for fitting the underlying linear model used by
#'  FM-HighLD
#' @param response a vector for the `singletrait` case or a named list for
#'  the multitrait case with an element for each variant. In the multitrait
#'  case, each element of the list is a named vector where each element
#'  is named after the trait for which the association was tested.
#' @param annot_matrix a matrix with same number of rows as variants and number
#'  of columns as annotations
#' @param ld_clusters a character vector with the ld cluster to which each
#'  variant belongs
#' @param singletrait  a logical indicator determining if the model is for
#' multi_trait or single-trait fine-mapping
#' @param fm_param a `FMParam` object with the parameters used to run `FMHighLD`
#' @param saveIter a logical indicator determining whether the iteration data is
#'  going to be returned
#' @param verbose a logical indicator determining whether messages are going to
#'  be used
#' @return results
fmhighld_fit <- function(formula, response, annot_matrix, ld_clusters,
  singletrait = TRUE, fm_param = FMParam(), saveIter = FALSE,
  verbose = FALSE) {

  stopifnot(
    length(response) == nrow(annot_matrix),
    length(response) == length(ld_clusters))

  same_names <- check_variant_names(response, annot_matrix, ld_clusters)
  if (! same_names) {
    warning("There variants are not named or the names don't match, will match
      'response', 'annot_matrix' and 'ld_clusters' by position")
  }


  # init algorithm parameteres
  iter <- 0
  error_bound <- error_bound(fm_param)
  max_iter <- max_iter(fm_param)
  min_tol <- min_tol(fm_param)

  if (saveIter) {
    states <- list()
  } else {
    states <- NULL
  }

  init <- init_iteration(formula, data, singletrait)
  model_list <- models(init)
  pi <- colMeans(probmatrix(init))
  background_error <- rlang::rep_along(model_list, 1)
  prev_causals <- dplyr::mutate(init_candidates, prob = 0.5)



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
