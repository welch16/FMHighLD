# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

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

// [[Rcpp::export]]
arma::mat lmm_comp_var(
		       arma::vec var_factors,
		       arma::vec residual_error,
		       arma::vec gene_error) {

  const int L = residual_error.n_elem;
  const int N = var_factors.n_elem;

  arma::mat var_mat(N,L);

  for(int i = 0; i < N; i++){
    for(int l = 0; l < L; l++){
      var_mat(i,l) = std::pow(residual_error(l),2) + var_factors(i) * std::pow(gene_error(l),2);
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
