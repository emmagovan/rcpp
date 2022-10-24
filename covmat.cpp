#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix cov_mat_cpp(NumericMatrix x, NumericMatrix y) {
NumericVector meanx(x.ncol());
NumericVector meany(y.ncol()); 
NumericMatrix covmat(x.ncol(), y.ncol());

for(int i = 0; i<x.ncol(); i++){
  meanx(i) = mean(x(_,i));
  }
for(int i = 0; i<y.ncol(); i++){
  meany(i) = mean(y(_,i));
  }

NumericMatrix xminusmean(x.nrow(), x.ncol());
NumericMatrix yminusmean(y.nrow(), y.ncol());

for(int j = 0; j<x.ncol(); j++){
  for(int i=0; i<x.nrow(); i++){
  xminusmean(i,j) = x(i,j) - meanx(j);
}
}

for(int j = 0; j<y.ncol(); j++){
  for(int i =0; i<y.nrow(); i++){
  yminusmean(i,j) = y(i,j) - meany(j);
  }
}

NumericMatrix sumxy(x.ncol(), y.ncol());

NumericVector xcol(x.ncol());
NumericVector ycol(y.ncol());

for(int i = 0; i<x.ncol(); i++){
  for(int j=0; j<y.ncol(); j++){
    for(int n =0; n<x.nrow(); n++){
      
      sumxy(i,j) += xminusmean(n,i) * yminusmean(n,j);
      }

    

    
    }}


for(int i=0; i<x.ncol(); i++){
  for(int j = 0; j<y.ncol(); j++){
    covmat(i,j) = sumxy(i,j)/(x.nrow()-1);
  }
}

return covmat;
}

