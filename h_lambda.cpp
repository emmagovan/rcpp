#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;


// This function takes theta and calculates the proportions
//[[Rcpp::export]]
NumericVector hfn(NumericVector theta, int n_sources){
  NumericVector p(n_sources);
  NumericVector exptheta(n_sources);
  double sumexptheta =0;
  
  // gets exp of each theta
  for(int i = 0; i<n_sources; i++){
    exptheta(i) = exp(theta(i));
    
  }
  
  // calculates sum of all exp thetas
  for(int i =0; i<n_sources; i++){
    sumexptheta +=exptheta[i];
    
  }
  // calculates p
  for(int i = 0; i<n_sources; i++){
    p[i] = exptheta[i]/sumexptheta;
    
  }
  
  return p;
  
}





//[[Rcpp::export]]
double hcpp(NumericVector p, int n_sources, int n_isotopes,
            NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
            NumericMatrix correctionmeans,
            NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix theta, NumericMatrix y ){
  double x =0;
  
  
  double ly = y.rows();
  
  // Setting prior values for hyper parameters
  NumericVector prior_means(n_sources); 
  NumericVector prior_sd(n_sources);
  NumericVector c_0(n_isotopes);
  NumericVector d_0(n_isotopes);
  
  // Setting up prior values
  for(int i=0; i<n_sources; i++){
    // mean was zero - this caused problems below when I needed to take the log of the mean
    prior_means(i) = 1;
    prior_sd(i) = 1;
  }
  
  for (int i = 0; i<n_isotopes; i++){
    c_0(i) = 1;
    d_0(i) = 1;
  }
  
  
  if(n_isotopes == 2){
    
    
    // This is to get dnorm(y[,1], sum(p*q) etc) 
    double mutop1 = 0;
    double mubtm1 = 0;
    double mu1 = 0;
    double sigmasq1 = 0;
    double sigmatopsq1 = 0;
    double sigmabtmsq1 = 0;
    
    // Calculate numerator and denominator of mu
    for(int i=0; i<n_sources; i++){
      mutop1 += p(i) * concentrationmeans(i,0) * (sourcemeans(i,0) + correctionmeans(i,0));
      mubtm1 += p(i) * concentrationmeans(i,0);
    }
    
    // Same for sigma
    for(int i=0; i<n_sources; i++){ // sds
      sigmatopsq1 += pow(p(i),2) * pow(concentrationmeans(i,0),2) * (pow(sourcesds(i,0),2) + 
        pow(corrsds(i,0),2));
      sigmabtmsq1 += pow(p(i),2) * pow(concentrationmeans(i,0),2);
    }
    
    //Calculate mu and sd
    mu1 = mutop1/mubtm1;
    sigmasq1 = sigmatopsq1/sigmabtmsq1;
    double sigma1 = pow(sigmasq1 + 1/theta(0, (n_sources)), 0.5);
    
    
    // This is to get dnorm(y[,2], sum(p*q) etc) 
    double mutop2 = 0;
    double mubtm2 = 0;
    double mu2 = 0;
    double sigmasq2 = 0;
    double sigmatopsq2 = 0;
    double sigmabtmsq2 = 0;
    for(int i=0; i<n_sources; i++){
      mutop2 += p(i) * concentrationmeans(i,1) * (sourcemeans(i,1) + correctionmeans(i,1));
      mubtm2 += p(i) * concentrationmeans(i,1);
    }
    
    
    for(int i=0; i<n_sources; i++){
      sigmatopsq2 += pow(p(i),2) * pow(concentrationmeans(i,1),2) * (pow(sourcesds(i,1),2) + 
        pow(corrsds(i,1),2));
      sigmabtmsq2 += pow(p(i),2) * pow(concentrationmeans(i,1),2);
    }
    mu2 = mutop2/mubtm2;
    sigmasq2 = sigmatopsq2/sigmabtmsq2;
    
    double sigma2 = pow(sigmasq2 + 1/theta(0, (1+n_sources)), 0.5);
    
    double yminusmu1 = 0;
    double yminusmu2 = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu1 += pow((y(i,0) - mu1),2);
      yminusmu2 +=  pow((y(i,1) - mu2),2);
    }
    
    
    
    // This is log(dnorm(y, p*q, p^2*q^2 etc) for y1 and y2
    
    x = - ly * log(sigma1) - 0.5 * ly * log(2 * M_PI)
      - 0.5 * yminusmu1 * 1/(pow(sigma1,2))
      - ly * log(sigma2) - 0.5 * ly * log(2 * M_PI)
      - 0.5 * yminusmu2 * 1/(pow(sigma2,2));
      
      
      
      
      
      
      
      
  }  else{
    //This is for just one isotope!!
    
    
    
    
    
    // This is to get dnorm(y[,1], sum(p*q) etc)
    double mutop = 0;
    double mubtm = 0;
    double mu = 0;
    double sigmasq = 0;
    double sigmatopsq = 0;
    double sigmabtmsq = 0;
    
    // calculate mu numerator and denominator
    for(int i=0; i<n_sources; i++){
      mutop += p(i) * concentrationmeans(i,0) * (sourcemeans(i,0) + correctionmeans(i,0));
      mubtm += p(i) * concentrationmeans(i,0);
    }
    
    
    
    
    // same for sigma
    for(int i=0; i<n_sources; i++){
      sigmatopsq += pow(p(i),2) * pow(concentrationmeans(i,0),2) * (pow(sourcesds(i,0),2) + 
        pow(corrsds(i,0),2));
      sigmabtmsq += pow(p(i),2) * pow(concentrationmeans(i,0),2);
    }
    
    
    
    //Calculate mu and sd
    mu = mutop/mubtm;
    sigmasq = sigmatopsq/sigmabtmsq;
    double sigma = pow(sigmasq + 1/theta(0, n_sources), 0.5);
    
    
    
    
    
    double yminusmu = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu += pow((y(i,0) - mu),2);
    }
    
    
    
    
    
    // This has y
    
    x = - ly * log(sigma) - 0.5 * ly * log(2 * M_PI) - 0.5 * yminusmu* 1/(pow(sigma,2));
    
    
    
    
  }
  
  double thetanorm = 0;
  
  
  
  
  for(int i = 0; i<n_sources; i++){
    thetanorm +=  - n_sources * log(prior_means(i)) - 0.5 * log(2 * M_PI) - (pow((theta(0,i) - prior_means(i)), 2)
                                                                               * 1/(2 * pow(prior_sd(i), 2)));
  }
  
  double gammaprior = 0;
  for (int i=0; i <(n_isotopes); i++){
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +(c_0(i) - 1) * theta(0,(i+n_sources)) - 
      d_0(i) * theta(0,(i+n_sources));
    
  }
  
  double totx = x + gammaprior + thetanorm;
  
  return totx;
  
  
}












//log q function

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
  
  return x;
  
}






// [[Rcpp::export]]
double h_lambdacpp(NumericVector p, int n_sources, int n_isotopes,
                          NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                          NumericMatrix correctionmeans,
                          NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix theta, NumericMatrix y,
                          NumericVector lambda) {
  
 
  
  
  
  
  
  
  return hcpp(p, n_sources, n_isotopes, concentrationmeans, sourcemeans, correctionmeans,
              corrsds, sourcesds, theta, y) - log_q_cpp(theta, lambda, n_sources, n_isotopes);
}












// h_lambda(lambda, theta[1,], y)