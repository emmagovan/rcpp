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

/*** R
OLDrMVNormC <- function(n, mu, U){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  # U <- chol(Omega) # By default R's chol fxn returns upper cholesky factor
  X <- backsolve(U, Z) # more efficient and stable than actually inverting
  X <- sweep(X, 1, mu, FUN=`+`)
  return(t(X))
}
MASSrMVN <- function(n, mu, U){
  S <- solve(crossprod(U))
  return(MASS::mvrnorm(n, mu, S))
}
mvnfastrMVN <- function(n, mu, U){
  S <- solve(crossprod(U))
  return(mvnfast::rmvn(n, mu, S))
}

U <- matrix(c(1, 0, 0, 4, 5, 0, 7, 8, 9), 3, 3)
mu <- 1:3
ansOLD <- OLDrMVNormC(1000, mu, U)
ansNEW <- rMVNormC(1000, mu, U)
microbenchmark::microbenchmark(
  OLDrMVNormC(1000, mu, U),
  rMVNormC(1000, mu, U),
  MASSrMVN(1000, mu, U),
  mvnfastrMVN(1000, mu, U)
)
*/
