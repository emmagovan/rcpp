#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
List log_q_cpp(NumericMatrix theta, NumericVector lambda, 
                        int n_sources, int n_tracers){
  
  NumericVector mean(n_sources);
  
  NumericVector thetavec = theta(1,_);
  
  for(int i = 0; i <n_sources; i++){
    mean(i) = thetavec(i);
    }
  
  NumericMatrix thetaminusmean(1, n_sources);
  for(int i = 0; i <n_sources; i++){
    thetaminusmean(1,i) = thetavec(i) - mean(i);
    }
  
  NumericMatrix choldecomp(n_sources, n_sources);
  int count = 0;
  for(int i = 0; i< n_sources; i++){ 
    for(int j = 0; j<n_sources; j++){
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
    for (int j = 0; j < thetaminusmean.ncol(); j++) 
    {
      p1(i,j)=thetaminusmean(i,j) * tcholdecomp(i,j) ;
    }
  }

  double p1tp1;
  for (int i = 0; i < p1.nrow(); i++) 
  {
    p1tp1 += pow(p1(1,i),2);
  }
  
  
    

    double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) * log(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i)) 
    - log(tgamma(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)))) 
    +(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) - 1) * theta(0,(i+n_sources)) - 
    lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i) * theta(0,(i+n_sources));
  }
 
 return -0.5 * n_sources * log(2 * M_PI) - 0.5 * sum(log(diag(choldecomp))) 
   - 0.5 * p1tp1 +gamman;
    
 
}
  
  
  
  
  
  