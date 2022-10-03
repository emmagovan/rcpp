#include <Rcpp.h>
using namespace Rcpp;



 // x <- t(rMVNormC(S, mu = mean, U = chol_prec))
// i think insert r code to do first part here
// then call it into the function and output
// r code and maths all bound together
//should be able to do with a loop

//[[Rcpp::export]]
NumericMatrix sim_thetacpp(int S, NumericVector lambda, int n_sources, 
                           int n_tracers, NumericMatrix rMVNorm){
NumericMatrix theta(S, (n_sources + n_tracers));
  
  for (int i = 0; i<n_sources; i++){
    
      theta(_,i) = rMVNorm(_,i);
    
}
  NumericMatrix gammam(S, n_tracers);
  
  for(int i = 0; i<n_tracers; i++){
 gammam(_,i) = Rcpp::rgamma(S,  lambda((n_sources + (n_sources * (n_sources + 1)) / 2) + i),
                              1/lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i));
  }
  
  
    for(int i = 0; i < (n_tracers); i++){
      theta(_,i+n_sources) = gammam(_,i);
      }

  return theta;
  
  }