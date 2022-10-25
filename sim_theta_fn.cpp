#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]

arma::mat rMVNormC(const double n,
                   const arma::vec mu,
                   const arma::mat U) {
  
  // Dimension of MVN
  int p = mu.size();
  
  // Simulate iid standard normals
  arma::mat Z(p, n);
  Z.imbue(norm_rand);
  
  // Now backsolve and add back on the means
  arma::mat X = solve(U, Z);
  for ( int i = 0; i < n; ++i ) {
    X.col(i) += mu;
  }
  
  return X.t();
}



 // x <- t(rMVNormC(S, mu = mean, U = chol_prec))
// i think insert r code to do first part here
// then call it into the function and output
// r code and maths all bound together
//should be able to do with a loop

//[[Rcpp::export]]
NumericMatrix sim_thetacpp(int S, NumericVector lambda, int n_sources, 
                           int n_tracers){
NumericMatrix theta(S, (n_sources + n_tracers));
  
NumericVector mean(n_sources);

for(int i = 0; i<n_sources; i++){
  mean(i) = lambda(i);
  }

  
  NumericMatrix chol_prec(n_sources, n_sources);
  int count = 0;
  for(int j = 0; j< n_sources; j++){ 
    for(int i = 0; i<n_sources; i++){
      if (i <= j){
        count +=1;
        chol_prec((i),(j)) = lambda(n_sources -1 +count);
        
        
      }
      else{
        chol_prec(i,j) = 0;
      }
      
    }
  }
  NumericMatrix tchol_prec = transpose(chol_prec);
  
  // 
  // NumericMatrix rMVNorm(theta.nrow(), n_sources);
  // 
  // // now need to make rMVNormC(S, mu = mean, U = chol_preg)
double p = (mean.length());
  // NumericMatrix Z(p, S);
  // for(int i = 0; i<(p*S); i++){
  //   Z(_,i) = Rcpp::rnorm(p);
  //   }
  // 
  // x = U^-1 Z
  
  
  
  for (int i = 0; i<n_sources; i++){
    
      theta(_,i) = rMVNormC(S*p, mean, tcholprec);
    
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



