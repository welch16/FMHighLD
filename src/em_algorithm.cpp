# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' E-step
//' 
//' Calculates the E-step probabilities in the implementation
//' of the EM-algorithm for mixture models used by FM-HighLD
//' 
//' @param y a vector with the association statistics with length n
//' @param pi a vector with the prior prob. of allocation into the
//'   components with length k
//' @param mu a nxk matrix with the conditional means of Y | Z = k
//' @param sigma a nxk matrix with the conditional std. dev of Y | Z = k
//' @return a nxk matrix with entries
//'   P(Z_i = k | Y_i = y_i , mu = mu_ik ,sigma = sigma_ik)
//' @export
// [[Rcpp::export]]
arma::mat estep(
		arma::vec y,
		arma::vec pi,
		arma::mat mu,
		arma::mat sigma) {
  
  const int N = mu.n_rows;
  const int K = mu.n_cols;
  
  arma::vec sum_vector(N);
  arma::mat gamma_mat(N,K);
    
  for(int i = 0; i < N ; i++){
    for(int j = 0; j < K; j++){
      gamma_mat(i,j) = pi(j) * R::dnorm4(y(i),mu(i,j),sigma(i,j),0);
    }
    sum_vector(i) = sum(gamma_mat.row(i));
    for(int j = 0; j < K; j++ ){
      gamma_mat(i,j) = gamma_mat(i,j) / sum_vector(i);
    }
  }
 
  return gamma_mat;
}


//' Computes the variance of a linear mixed model
//'
//' @param var_factors a vector with the variances for each component of
//'   the linear mixed model random effects
//' @param residual_error a vector with the residual errors of the
//'   linear mixed model used by FM-HighLD
//' @param re_error a vector with the random effect erros of the linear
//'   mixed model used by FM-HighLD
//' @return The variance matrix of the model
// [[Rcpp::export]]
arma::mat lmm_comp_var(
		       arma::vec var_factors,
		       arma::vec residual_error,
		       arma::vec re_error) {

  const int L = residual_error.n_elem;
  const int N = var_factors.n_elem;

  arma::mat var_mat(N,L);

  for(int i = 0; i < N; i++){
    for(int l = 0; l < L; l++){
      var_mat(i,l) = std::pow(residual_error(l),2) + var_factors(i) * std::pow(re_error(l),2);
    }    
  }
  
  return var_mat;

}

// arma::vec mstep_pi(arma::mat gamma)
// {
//   // const int K = gamma.n_cols;
//   const int N = gamma.n_rows;
//   
//   arma::vec pi = trans(sum(gamma,0))/ N;
//   
//   // for(int j = 0; j < K; j++){
//   //   pi(j) = sum(gamma.col(j)) / N;
//   // }
//   
//   
//   
//   return pi;
// }

// /*** R
// y = rnorm(5)
// sigma2 = rexp(3)
// mu = matrix(rnorm(15),nrow = 5)
// pi = c(3,1,1)
// ll = estep(y,pi/sum(pi),sigma2, mu)
// # pp = mstep_pi(ll)
// */
