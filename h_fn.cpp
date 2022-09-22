#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector hfn(NumericVector theta, int n_sources){
  NumericVector p(n_sources);
  NumericVector exptheta(n_sources);
  double sumexptheta =0;
  
  
  for(int i = 0; i<n_sources; i++){
    exptheta(i) = exp(theta(i));
    
  }
  
  
  for(int i =0; i<n_sources; i++){
    sumexptheta +=exptheta[i];
    
  }
  
  for(int i = 0; i<n_sources; i++){
    p[i] = exptheta[i]/sumexptheta;
    
  }
  
  
  
  
  
  
  return p;
  
}


//[[Rcpp::export]]
List hcpp(NumericVector p, int n_sources, int n_isotopes,
         NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
         NumericMatrix correctionmeans,
         NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix theta, NumericMatrix y ){
  double x =0;

  
  double ly = y.rows();
  NumericVector prior_means(n_sources); //length???
  NumericVector prior_sd(n_sources);
  NumericVector c_0(n_isotopes);
  NumericVector d_0(n_isotopes);
  
  // Setting up prior values
  for(int i=0; i<n_sources; i++){
    prior_means(i) = 1;
    prior_sd(i) = 1;
    }
  
  for (int i = 0; i<n_isotopes; i++){
    c_0(i) = 1;
    d_0(i) = 1;
    }
    
  
  if(n_isotopes == 2){
    //double pi = 3.14159;
    
    // This is to get dnorm(y[,1], sum(p*q) etc) and for y[,2] as well
    double mutop1 = 0;
    double mubtm1 = 0;
    double mu1 = 0;
    double sigma1 = 0;
    double sigmatop1 = 0;
    double sigmabtm1 = 0;
    for(int i=0; i<n_sources; i++){
      mutop1 += p(i) * concentrationmeans(i,0) * (sourcemeans(i,0) + correctionmeans(i,0));
      mubtm1 += p(i) * concentrationmeans(i,0);
    }
    
    
    for(int i=0; i<n_sources; i++){ // sds
      sigmatop1 += pow(p(i),2) * pow(concentrationmeans(i,0),2) * (pow(sourcesds(i,0),2) + 
        pow(corrsds(i,0),2));
      sigmabtm1 += pow(p(i),2) * pow(concentrationmeans(i,0),2);
    }
    mu1 = mutop1/mubtm1;
    sigma1 = sigmatop1/sigmabtm1;
    
    
    // This is to get dnorm(y[,2], sum(p*q) etc) 
    double mutop2 = 0;
    double mubtm2 = 0;
    double mu2 = 0;
    double sigma2 = 0;
    double sigmatop2 = 0;
    double sigmabtm2 = 0;
    for(int i=0; i<n_sources; i++){
      mutop2 += p(i) * concentrationmeans(i,1) * (sourcemeans(i,1) + correctionmeans(i,1));
      mubtm2 += p(i) * concentrationmeans(i,1);
    }
    
    
    for(int i=0; i<n_sources; i++){
      sigmatop2 += pow(p(i),2) * pow(concentrationmeans(i,1),2) * (pow(sourcesds(i,1),2) + 
        pow(corrsds(i,1),2));
      sigmabtm2 += pow(p(i),2) * pow(concentrationmeans(i,1),2);
    }
    mu2 = mutop2/mubtm2;
    sigma2 = sigmatop2/sigmabtm2;
    
    double yminusmu1 = 0;
    double yminusmu2 = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu1 += y(i,0) - mu1;
      yminusmu2 += y(i,1) - mu2;
    }
    
    
    
    // This has y1 and y2 so far, just need to add the rest
    
    x = - ly * log(sigma1) - 0.5 * log(pow(M_PI,2)) - (pow((yminusmu1),2)) * 1/(2*(pow(sigma1,2))) 
      - ly * log(sigma2) - 0.5 * log(pow(M_PI,2)) - (pow((yminusmu2),2)) * 1/(2*(pow(sigma2,2)));
    

      
      
  }  else{
    //This is for just one isotope!!
    
    
    
    //double pi = 3.14159;
    
    // This is to get dnorm(y[,1], sum(p*q) etc) and for y[,2] as well
    double mutop = 0;
    double mubtm = 0;
    double mu = 0;
    double sigma = 0;
    double sigmatop = 0;
    double sigmabtm = 0;
    for(int i=0; i<n_sources; i++){
      mutop += p(i) * concentrationmeans(i,0) * (sourcemeans(i,0) + correctionmeans(i,0));
      mubtm += p(i) * concentrationmeans(i,0);
    }
    
    
    for(int i=0; i<n_sources; i++){
      sigmatop += pow(p(i),2) * pow(concentrationmeans(i,0),2) * (pow(sourcesds(i,0),2) + 
        pow(corrsds(i,0),2));
      sigmabtm += pow(p(i),2) * pow(concentrationmeans(i,0),2);
    }
    mu = mutop/mubtm;
    sigma = sigmatop/sigmabtm;
    double sigma_sq = pow(sigma,2);
    
    double yminusmu = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu += y(i,0) - mu;
      }
    
    
    
    
    
    // This has y
    
    x = - ly * log(mu) - 0.5 * log(pow(M_PI,2)) - (pow((yminusmu),2))* 1/(2*sigma_sq);
    
  }
  
  double thetanorm = 0;

  for(int i = 0; i<n_sources; i++){
    thetanorm +=  - n_sources * log(prior_means(i)) - 0.5 * log(pow(M_PI,2)) - (pow((theta(0,i) - prior_means(i)), 2)
    * 1/(2 * pow(prior_sd(i), 2)));
    }

  double gammaprior = 0;
  for (int i=0; i <(n_isotopes); i++){
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +(c_0(i) - 1) * theta(0,(i+n_sources)) - 
      d_0(i) * theta(0,(i+n_sources));

    }

  x = x + gammaprior + thetanorm;

  return Rcpp::List::create(Rcpp::Named("total") = x,
                            Rcpp::Named("gamma") = gammaprior,
                            Rcpp::Named("theta") = thetanorm);

  
}













