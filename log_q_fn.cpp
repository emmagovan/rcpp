#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
List log_q_cpp(NumericMatrix theta, NumericVector lambda, 
                        int n_sources, int n_tracers){
  

  
  NumericVector thetavec = theta(0,_);

  
  NumericMatrix thetaminusmean(1, n_sources);
  for(int i = 0; i <n_sources; i++){
    thetaminusmean(0,i) = thetavec(i) - lambda(i);
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
    for (int j = 0; j < tcholdecomp.ncol(); j++) 
      for (int k = 0; k < tcholdecomp.nrow(); k++){
    {
      p1(i,j) =p1(i,j) + thetaminusmean(i,k) * tcholdecomp(k,j) ;
    }
  }
  }
  



  double p1tp1 = 0;
  for (int i = 0; i < n_sources; i++) 
  {
    p1tp1 += pow(p1(0,i),2);
  }
  
  
    

    double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) * log(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i)) 
    - log(tgamma(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)))) 
    +(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) - 1) * theta(0,(i+n_sources)) - 
    lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i) * theta(0,(i+n_sources));
  }
  
  double sumlogdiag = 0;
  for(int i = 0; i<choldecomp.nrow(); i++){
    for(int j = 0; j<choldecomp.nrow(); j++){
      if((i = j)){
        sumlogdiag += log(choldecomp(i,j));
        }
      
      }}
  
 
 double x = -0.5 * n_sources * log(2 * M_PI) - 0.5 * sumlogdiag 
   - 0.5 * p1tp1 +gamman;
    
    return Rcpp::List::create(Rcpp::Named("total") = x,
                        Rcpp::Named("p1") = p1,
                        Rcpp::Named("p1tp1") = p1tp1,
                        Rcpp::Named("gammam") = gamman,
                        Rcpp::Named("thetavec") = thetavec,
                        Rcpp::Named("sumulogdiag") = sumlogdiag,
                        Rcpp::Named("thminusmean") = thetaminusmean,
                        Rcpp::Named("lambda") = lambda,
                        Rcpp::Named("choldecomp") = choldecomp,
                        Rcpp::Named("tcholdecomp") = tcholdecomp
    )
                        ;
 
}
  
  
  
  
  
  