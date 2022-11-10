#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double test(double eps_0, double tau, double t){
NumericVector alpha_min(2);
double alpha_t;

for(int i=0; i<2; i++){
  alpha_min(0) = eps_0;
  alpha_min(1) = eps_0 * tau/(t+1);
}

alpha_t = Rcpp::min(alpha_min);

return alpha_t;
}