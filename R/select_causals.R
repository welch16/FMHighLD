
#' Checks whether the intercept is include in the underlying linear model
#' @param features names of the estimated coefficient
#' @return a logical indicating whether include the intercept in the model
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
#' @importFrom stats formula terms
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
select_causals_multi <- function(data, model, group, rand_param = NULL,
  cond_res = FALSE) {

  include_intercept = beta %>%
        names() %>%
        str_detect("Intercept") %>%
        any()
  model_formula <- stats::formula(model)
  response <- as.character(model_formula)[2]

    include_intercept = formula(model) %>%
        terms() %>% attr("intercept") 
        
    beta = fixef(model)
    
    des_mat = build_design_matrix(new_data,beta,include_intercept)
    fitted = as.numeric(des_mat %*% beta)

    res = pluck(new_data,response) - fitted
    


    causal_rule(new_data,res,rand_param,to_group)
}


causal_rule <- function(new_data,res,rand_param,to_group)
{
    aux_data = new_data %>%
        select("SNP",to_group) %>%
        mutate(
            res.sq = res^2) %>%
        group_by(.dots = to_group)
    
    if(is.null(rand_param)){
        out = aux_data %>% 
            summarize(
                which_snp = SNP[nnet::which.is.max(-res.sq)]
            )  %>%
          ungroup()    
    }else if(rand_param$strat == "all"){
        
        out = aux_data %>% 
            summarize(
                which_snp = SNP[select_kth_random(res.sq,rand_param$prob,rand_param$select)]
            ) %>%
            ungroup()
        
    }else if(rand_param$strat == "pick_M"){
    
        out = aux_data %>%
            summarize(
                n = n(),
                which_snp = SNP[nnet::which.is.max(-res.sq)])
        
        if(rand_param$which == "any"){
            picks = out %>%
                ungroup() %>% 
                filter( n > 1) %>%
                select(-n,-which_snp) %>% 
            sample_n(rand_param$M)
            
        }else if(rand_param$which == "largeLD"){
            picks = out %>%
                ungroup() %>% 
                arrange(desc(n)) %>%
                head(rand_param$M) %>%
                select(-n,-which_snp)
            
        }
        
        pick_data = picks %>%
            inner_join(
                aux_data,by = to_group) %>%
          group_by(.dots = to_group) %>% 
            summarize(
                which_snp = SNP[select_kth_random(res.sq,rand_param$prob,rand_param$select)]
            ) %>%
            ungroup()
        
        out = out %>%
            ungroup() %>%
            select(-n) %>%
          anti_join(picks,by = to_group) %>%
            bind_rows(pick_data)
        
    }
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