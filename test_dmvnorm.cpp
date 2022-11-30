#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]

arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}
// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat x,  arma::rowvec mean,  arma::mat sigma, bool log = false) { 
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = std::log(2.0 * M_PI);
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if (log){ 
    return(logretval);
  }else { 
    return(exp(logretval));
  }
}

// [[Rcpp::export]]
NumericVector wrap(NumericMatrix x, NumericVector mean, NumericMatrix sigma){
  NumericVector ans(x.ncol());
  
 ans = Rcpp::wrap((dmvnorm(as<arma::mat>(x), mean, as<arma::mat>(sigma), true)));
  
  
  return ans;
}




