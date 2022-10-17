#include <Rcpp.h>
using namespace Rcpp;

//log q function

//[[Rcpp::export]]
#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
double log_q_cpp(NumericMatrix theta, NumericVector lambda, 
               int n_sources, int n_tracers){
  
  
  
  NumericVector thetavec = theta(0,_);
  
  
  NumericMatrix thetaminusmean(1, n_sources);
  for(int i = 0; i <n_sources; i++){
    thetaminusmean(0,i) = thetavec(i) - lambda(i);
  }
  
  NumericMatrix choldecomp(n_sources, n_sources);
  int count = 0;
  for(int j = 0; j< n_sources; j++){ 
    for(int i = 0; i<n_sources; i++){
      if (i <= j){
        count +=1;
        choldecomp((i),(j)) = lambda(n_sources -1 +count);
        
        
      }
      else{
        choldecomp(i,j) = 0;
      }
      
    }
  }
  NumericMatrix tcholdecomp = transpose(choldecomp);
  
  NumericMatrix p1(1, n_sources);
  
  for (int i = 0; i < thetaminusmean.nrow(); i++) 
  {
    for (int j = 0; j < tcholdecomp.ncol(); j++) 
      for (int k = 0; k < tcholdecomp.nrow(); k++){
        {
          p1(i,j) =p1(i,j) + thetaminusmean(i,k) * tcholdecomp(k,j) ;
        }
      }
  }
  
  
  NumericMatrix pt = transpose(p1);
  
  double p1tp1 = 0;
  for (int i = 0; i < p1.nrow(); i++) 
  {
    for (int j = 0; j < pt.ncol(); j++) 
      for (int k = 0; k < pt.nrow(); k++){
        {
          p1tp1 += p1(i,k) * pt(k,j) ;
        }
      }
  }
  
  
  
  
  double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) * 
      log(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i)) 
    - log(tgamma(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)))) 
    +(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) - 1) * log(theta(0,(i+n_sources))) - 
    lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i) * theta(0,(i+n_sources));
  }
  
  double sumlogdiag = 0;
  for(int i = 0; i<n_sources; i++){
    for(int j = 0; j<n_sources; j++){
      if(i == j){
        sumlogdiag += log(choldecomp(i,j));
      }
    }}
  
  
  double x = -0.5 * n_sources * log(2 * M_PI) - 0.5 * sumlogdiag 
    - 0.5 * p1tp1 + gamman;
  
  return x
    ;
  
}















// [[Rcpp::export]]
List delta_lqltcpp(NumericVector lambda, NumericMatrix theta, double eps, int n_sources, int n_tracers) {
  eps = 0.001;
  double k = lambda.length();
  NumericVector ans(k);
  NumericVector d(k);
  NumericVector lambdaplusd(k);
  NumericVector lambdaminusd(k);
  

  
  for(int i = 0; i<k; i++){
    
    for (int j = 0; j<k; j++){
      d(j) = 0;
    }
    d(i) = eps;

    
    for (int j = 0; j<k; j++){
      lambdaplusd(j) = lambda(j) + d(j);
      lambdaminusd(j) = lambda(j) - d(j);
    }
    ans(i) = (log_q_cpp(theta, lambdaplusd, n_sources, n_tracers) -  
      log_q_cpp(theta, lambdaminusd, n_sources, n_tracers))/(2 * eps);
  }
  return  Rcpp::List::create(Rcpp::Named("ans") = ans,
                             Rcpp::Named("lambdaplusd") = lambdaplusd,
                             Rcpp::Named("d") = d,
                             Rcpp::Named("lambdaminusd") = lambdaminusd);
}


// delta_lqlt <- function(lambda, theta, eps = 0.001) {
//   k <- length(lambda)
//   ans <- rep(NA, k)
//   for (i in 1:k) {
//     d <- rep(0, k)
//     d[i] <- eps
//     ans[i] <- (log_q(lambda + d, theta) - log_q(lambda - d, theta)) / (2 * max(d))
//   }
//   return(ans)}