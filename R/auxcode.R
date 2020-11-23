
rebuild_data <- function(new_data,causals,to_group)
{

    causals %>%
        dplyr::rename(SNP = which_snp) %>%
        inner_join(new_data, by = c("SNP",to_group))
    

}



parse_estimated_parameter_multi <- function(em_iteration,data,response)
{
    gamma_matrix = em_iteration %>%
        pluck("gamma")

    gamma_sums = gamma_matrix %>%
        colSums()
    
    pi = gamma_sums / sum(gamma_sums)
        
    models = em_iteration %>%
        pluck("models")
    
    param_list = models %>%
        map(fixef)

    background_error = pluck(data,response) %>%
        {crossprod(.^2,gamma_matrix[,1])} %>%
        {  . / gamma_sums[1]} %>%
        as.numeric()  %>% sqrt()

    X_matrices = models %>%
        map(model.matrix)

    N = X_matrices %>%
        map_dbl(nrow)

    include_intercept = models %>%
        map(formula) %>%
        map(terms) %>%
        map_int(attr,"intercept")
    
    P = X_matrices %>%
        map2_dbl(include_intercept, ~ ncol(.x) - ( 1 - .y))

    factors = rep_along(gamma_sums,N - P) / gamma_sums


    error_list = models %>%
        map(VarCorr) %>%
        map(as.data.frame) %>%
        map(as_tibble) %>%
        map2( factors[-1],
             ~ mutate(.x, mix_vcov = .y * vcov))

    residual_error_list = error_list %>%
        map(filter,grp == "Residual") %>%
        map(pull,mix_vcov)

    gene_error_list = error_list %>%
        map(filter, grp != "Residual") %>%
        map(pull,mix_vcov)

    ## random_effect_list = map(models,ranef) %>%
    ##     map(as.data.frame) %>%
    ##     map(as_tibble)
    
    list(
        pi = pi,
        param_list = param_list,
        gene_error_list = gene_error_list,
        residual_error_list = residual_error_list,
        background_error = background_error)
}

parse_estimated_parameter_single <- function(em_iteration,data,response)
{
    gamma_matrix = em_iteration %>%
        pluck("gamma")

    gamma_sums = gamma_matrix %>%
        colSums()
    
    pi = gamma_sums / sum(gamma_sums)
        
    models = em_iteration %>%
        pluck("models")
    
    param_list = models %>%
        map(coef)

    background_error = pluck(data,response) %>%
        {crossprod(.^2,gamma_matrix[,1])} %>%
        {  . / gamma_sums[1]} %>%
        as.numeric()  %>% sqrt()

    X_matrices = models %>%
        map(model.matrix)

    N = X_matrices %>%
        map_dbl(nrow)

    include_intercept = models %>%
        map( ~ formula(.)) %>%
        map(terms) %>% map_int(attr,"intercept")

    P = X_matrices %>%
        map2_dbl(include_intercept, ~ ncol(.x) - (1 - .y))
    
    factors = rep_along(gamma_sums,N - P) / gamma_sums

    residual_error_list = models %>%
        map(broom::glance) %>%
        map_dbl("sigma") %>% 
        map2( factors[-1],
             ~ .x * .y)
    
    list(
        pi = pi,
        param_list = param_list,
        residual_error_list = residual_error_list,
        background_error = background_error)

}

multiply_factor_corrected <- function(annot_mat, model_coeff) {

  colnames(annot_mat)[1] <- names(model_coeff)[1]
  annot_mat[, names(model_coeff)] %*% model_coeff

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


lmm_mixture <- function(train_data,
                        model_formula,
                        features = c("ATAC","TF_alle"),
                        response = "eQTL_tStat",
                        max_iter = 100,
                        tol = 1e-6,save_state = FALSE,
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
                          train_data,
                          response,FALSE)

    model_list = init$models
    pi = init$gamma %>% colMeans()
    background_error = rep_along(model_list,1)
    gamma = init$gamma[,2]
     
    while(continue){

        if(verbose) message(iter)
        em_iter = em_iteration_multi(
            formula = model_formula,
            data = train_data,
            response = response,
            pi = pi,
            model_list  = model_list,
            background_error = background_error,
            REML = TRUE,
            verbose = verbose)
        if(verbose) message("EM step done")
        
        params = parse_estimated_parameter_multi(
            em_iter,train_data,response) 
        
        error_list = map(model_list,broom::tidy) %>%
            map(filter,group != "fixed")
        
        param_diff = map2(
            params,
            list(
                pi = pi,
                param_list = map(model_list,fixef),
                gene_error_list = map(error_list, filter,group != "Residual") %>%
                    map(pull,estimate),
                residual_error_list = map(error_list, filter,group == "Residual") %>%
                    map(pull,estimate),
                background_error = background_error),
            ~ ( unlist(.x) - unlist(.y))^2) %>%
            map_dbl(sum) 
        
        param_metric = sum(param_diff) / ( 1 + max(param_diff)) %>% sqrt()
        
        pi = params$pi
        model_list = pluck(em_iter,"models")
        background_error = params$background_error 
       
        weights = pluck(em_iter,"gamma") %>% {.[,2]}

        entropy = kl(weights,gamma)
        jdist = jsd(weights,gamma)

        gamma = weights

        metrics = c("js_dist" = jdist)

        continue = any(abs(metrics) >= tol) & iter <= max_iter
        iter = iter + 1
                  
        metrics["philips"] = param_metric        

        
        metrics["err"] = model_list %>%
            pluck(1) %>% broom::glance() %>%
            pluck("sigma")

        metrics["kl"] = entropy
        if(save_state){
            ## states[[iter]] = param_list
            states[[iter]] = list(
                "param" = params,
                "metrics" = metrics)                        
        }

    }

    if(verbose) message("Algorithm done")
    
    metrics["iter_exit"] = iter - 1

    ranef_list = model_list %>%
        map(ranef) %>%
        map(as.data.frame) %>%
        map(as_tibble) %>%
        map(spread,term,condval)

    
    list(
        estimated_parameters = map(params,unlist,recursive = FALSE),
        metrics = metrics,
        gamma = em_iter$gamma[,2],
        states = states,
        ranef = ranef_list)
    

}

generate_init <- function(seed,all_data)
{
    set.seed(seed)

    all_data %>%
        group_by(eQTL_gene,cluster) %>%
        sample_n(1) %>%
        ungroup()

}


perform_test_init <- function(train_data,rand_param,init,
							model_formula,
							features = "skew",
							verbose = FALSE)
{

    iterative_select_causal_init(
        train_data,
        model_formula,
        init,
        features,
        response = "stat",
        to_group = c("eQTL_gene","cluster"),
        cond_res = FALSE,
        rand_param = rand_param,
        max_iter = max_iter,
        tol = tol,
        save_state = FALSE,
        univariate = FALSE,
        verbose = verbose)
}
