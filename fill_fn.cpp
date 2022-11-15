#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector fill(NumericVector x) {
  NumericVector vec(100);
  for (int i = 0; i<vec.length(); i++){
    
  vec(i) = x(i % x.length());
 
  }
  return vec;
}

