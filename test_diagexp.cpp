#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix test(int n_sources, NumericVector lambda){
NumericMatrix diagexp(n_sources, n_sources);

for(int i =0; i<n_sources; i++){
  for(int j = 0; j<n_sources; j++){
    if(i ==j){
      diagexp(i,i) = exp(lambda(n_sources +i));
    }
    else{diagexp(i,j) = 0;
    }
  }
}
return diagexp;
}