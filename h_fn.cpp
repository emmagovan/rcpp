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
List hcpp(NumericVector p, int n_sources, int n_isotopes,
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
    double sigma1 = pow(sigmasq1, 0.5) + 1/theta(0, (n_isotopes));
    
    
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
    
    double sigma2 = pow(sigmasq2, 0.5) + theta(0, 1+n_isotopes);
    
    double yminusmu1 = 0;
    double yminusmu2 = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu1 += y(i,0) - mu1;
      yminusmu2 += y(i,1) - mu2;
    }
    
    
    
    // This is log(dnorm(y, p*q, p^2*q^2 etc) for y1 and y2
    
    x = - ly * log(sigma1) - 0.5 * ly * log(pow(M_PI,2)) 
        - 0.5 * (pow((yminusmu1),2)) * 1/(pow(sigma1,2)) 
        - ly * log(sigma2) - 0.5 * ly * log(pow(M_PI,2)) 
        - 0.5 * (pow((yminusmu2),2)) * 1/(pow(sigma2,2));
    

      
      
  }  else{
    //This is for just one isotope!!
    
    
    
   
    
    // This is to get dnorm(y[,1], sum(p*q) etc)
    double mutop = 0;
    double mubtm = 0;
    double mu = 0;
    double sigma = 0;
    double sigmatop = 0;
    double sigmabtm = 0;
    
    // calculate mu numerator and denominator
    for(int i=0; i<n_sources; i++){
      mutop += p(i) * concentrationmeans(i,0) * (sourcemeans(i,0) + correctionmeans(i,0));
      mubtm += p(i) * concentrationmeans(i,0);
    }
    
    // same for sigma
    for(int i=0; i<n_sources; i++){
      sigmatop += pow(p(i),2) * pow(concentrationmeans(i,0),2) * (pow(sourcesds(i,0),2) + 
        pow(corrsds(i,0),2));
      sigmabtm += pow(p(i),2) * pow(concentrationmeans(i,0),2);
    }
    
    //Calculate mu and sd
    mu = mutop/mubtm;
    sigma = sigmatop/sigmabtm;
    double sigma_sq = pow(sigma,2);
    
    double yminusmu = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu += y(i,0) - mu;
      }
    
    
    
    
    
    // This has y
    
    x = - ly * log(sigma) - 0.5 * ly * log(pow(M_PI,2)) - (pow((yminusmu),2))* 1/(2*sigma_sq);
    
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













