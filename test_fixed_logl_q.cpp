#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}



// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = true) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
NumericMatrix crossprod(NumericMatrix X){
  NumericMatrix ans(X.nrow(), X.ncol());
  
  for(int i = 0; i<X.ncol(); i++){
    for(int j=0; j<X.ncol(); j++){
      for(int n =0; n<X.nrow(); n++){
        
        ans(i,j) += X(n,i) * X(n,j);
      }
    }}
  return(ans);
}

// [[Rcpp::export]]
NumericMatrix rMVNormCpp(const double n,
                         const arma::vec mu,
                         const NumericMatrix U) {
  
  
  // Dimension of MVN
  int p = mu.size();
  
  // Simulate iid standard normals
  arma::mat Z(p, n);
  Z.imbue(norm_rand);
  
  // Now backsolve and add back on the means
  arma::mat X = solve(as<arma::mat>(U), Z);
  for ( int i = 0; i < n; ++i ) {
    X.col(i) += mu;
  }
  
  return Rcpp::wrap(X.t());
}


// [[Rcpp::export]]
NumericMatrix solvearma(const NumericMatrix X) {
  
  arma::mat b = arma::eye(X.nrow(), X.ncol());
  
  
  // Now backsolve and add back on the means
  arma::mat ans = solve(as<arma::mat>(X), b);
  
  
  return Rcpp::wrap(ans.t());
}




// x <- t(rMVNormC(S, mu = mean, U = chol_prec))
// i think insert r code to do first part here
// then call it into the function and output
// r code and maths all bound together
//should be able to do with a loop

//[[Rcpp::export]]
NumericMatrix sim_thetacpp(int S, NumericVector lambda, int n_sources, 
                           int n_tracers){
  NumericMatrix theta(S, (n_sources + n_tracers));
  
  NumericVector mean(n_sources);
  
  for(int i = 0; i<n_sources; i++){
    mean(i) = lambda(i);
  }
  
  
  NumericMatrix chol_prec(n_sources, n_sources);
  int count = 0;
  for(int j = 0; j< n_sources; j++){ 
    for(int i = 0; i<n_sources; i++){
      if (i <= j){
        count +=1;
        chol_prec((i),(j)) = lambda(n_sources -1 +count);
        
        
      }
      else{
        chol_prec(i,j) = 0;
      }
      
    }
  }
  //NumericMatrix tchol_prec = transpose(chol_prec);
  
  
  NumericMatrix normmat(S, n_sources);
  
  
  //for (int i = 0; i<n_sources; i++){
  
  // stop("Hello!");
  normmat = rMVNormCpp(S, mean, chol_prec);
  
  //}
  NumericMatrix gammam(S, n_tracers);
  
  for(int i = 0; i<n_tracers; i++){
    gammam(_,i) = Rcpp::rgamma(S,  lambda((n_sources + (n_sources * (n_sources + 1)) / 2) + i),
           1/lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i));
  }
  
  
  for(int i=0; i<n_sources; i++){
    theta(_,i) = normmat(_,i);
  }
  
  
  for(int i = 0; i < (n_tracers); i++){
    theta(_,i+n_sources) = gammam(_,i);
  }
  
  return theta;
  
}





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
double hcpp(int n_sources, int n_isotopes,
            NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
            NumericMatrix correctionmeans,
            NumericMatrix corrsds, NumericMatrix sourcesds, NumericVector theta, NumericMatrix y ){
  
  double x =0;
  
  NumericVector p(n_sources);
  
  p = hfn(theta, n_sources);
  
  double ly = y.rows();
  
  // Setting prior values for hyper parameters
  NumericVector prior_means(n_sources);
  NumericVector prior_sd(n_sources);
  NumericVector c_0(n_isotopes);
  NumericVector d_0(n_isotopes);
  
  // Setting up prior values
  for(int i=0; i<n_sources; i++){
    // mean was zero - this caused problems below when I needed to take the log of the mean
    prior_means(i) = 0;
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
    
    // NumericMatrix pxconc1(p.length(), n_isotopes);
    // for(int i=0; i<n_sources; i++){
    //   for(int j = 0; j<n_isotopes; j++){
    //     pxconc1(i,j) = p(i)*concentrationmeans(i,j);
    //   }}
    
    // Calculate numerator and denominator of mu
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_isotopes; j++){
        mutop1 +=  p(i)*concentrationmeans(i,j) * (sourcemeans(i,0) + correctionmeans(i,0));
        mubtm1 += p(i) * concentrationmeans(i,j);
      }
    }
    
    // Same for sigma
    for(int i=0; i<n_sources; i++){ // sds
      for(int j =0; j<n_isotopes; j++){
        sigmatopsq1 += pow(p(i),2) * pow(concentrationmeans(i,j),2) * (pow(sourcesds(i,0),2) +
          pow(corrsds(i,0),2));
        sigmabtmsq1 += pow(p(i),2) * pow(concentrationmeans(i,j),2);
      }
    }
    
    //Calculate mu and sd
    mu1 = mutop1/mubtm1;
    sigmasq1 = sigmatopsq1/sigmabtmsq1;
    double sigma1 = pow(sigmasq1 + 1/theta((n_sources)), 0.5);
    
    
    // This is to get dnorm(y[,2], sum(p*q) etc)
    double mutop2 = 0;
    double mubtm2 = 0;
    double mu2 = 0;
    double sigmasq2 = 0;
    double sigmatopsq2 = 0;
    double sigmabtmsq2 = 0;
    for(int i=0; i<n_sources; i++){
      for(int j =0; j<n_isotopes; j++){
        mutop2 += p(i) * concentrationmeans(i,j) * (sourcemeans(i,1) + correctionmeans(i,1));
        mubtm2 += p(i) * concentrationmeans(i,j);
      }
    }
    
    
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_isotopes; j++){
        sigmatopsq2 += pow(p(i),2) * pow(concentrationmeans(i,j),2) * (pow(sourcesds(i,1),2) +
          pow(corrsds(i,1),2));
        sigmabtmsq2 += pow(p(i),2) * pow(concentrationmeans(i,j),2);
      }
    }
    
    mu2 = mutop2/mubtm2;
    sigmasq2 = sigmatopsq2/sigmabtmsq2;
    
    double sigma2 = pow(sigmasq2 + 1/theta((1+n_sources)), 0.5);
    
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
    double sigma = pow(sigmasq + 1/theta(n_sources), 0.5);
    
    
    
    
    
    double yminusmu = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu += pow((y(i,0) - mu),2);
    }
    
    
    
    
    
    // This has y
    
    x = - ly * log(sigma) - 0.5 * ly * log(2 * M_PI) - 0.5 * yminusmu* 1/(pow(sigma,2));
    
    
    
    
  }
  
  double thetanorm = 0;
  
  
  
  
  for(int i = 0; i<n_sources; i++){
    thetanorm +=  - n_sources * log(prior_sd(i)) - 0.5 * log(2 * M_PI) - (pow((theta(i) - prior_means(i)), 2)
                                                                            * 1/(2 * pow(prior_sd(i), 2)));
  }
  
  double gammaprior = 0;
  for (int i=0; i <(n_isotopes); i++){
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +(c_0(i) - 1) * theta((i+n_sources)) -
      d_0(i) * theta((i+n_sources));
    
  }
  
  double totx = x + gammaprior + thetanorm;
  
  return totx;
  
}



//[[Rcpp::export]]
double log_q_cpp(NumericVector theta, NumericVector lambda, 
                 int n_sources, int n_tracers){
  
  NumericMatrix thetaminusmean(1, n_sources);
  
  for(int i = 0; i <n_sources; i++){
    thetaminusmean(0,i) = theta(i) - lambda(i);
  }
  
  NumericMatrix chol_prec(n_sources, n_sources);
  int count = 0;
  for(int j = 0; j< n_sources; j++){ 
    for(int i = 0; i<n_sources; i++){
      if (i <= j){
        count +=1;
        chol_prec((i),(j)) = lambda(n_sources -1 +count);
        
        
      }
      else{
        chol_prec(i,j) = 0;
      }
      
    }
  }
  NumericMatrix prec(n_sources, n_sources);
  prec = crossprod(chol_prec);
  
  NumericMatrix solve_prec(n_sources, n_sources);
  solve_prec = solvearma(prec);
  
  NumericMatrix y(1, n_sources);
  NumericVector mean(n_sources);
  
  for(int i = 0; i<n_sources; i++){
    y(0,i) = theta(i);
    mean(i) = lambda(i);
  }
  
  
  
  
  double thetanorm = 0;
  
  thetanorm = 
    *REAL(Rcpp::wrap(dmvnrm_arma_fast(as<arma::mat>(y), mean, as<arma::mat>(solve_prec))));
    
    
    
    
    
    double gamman = 0;
    for (int i=0; i <(n_tracers); i++){
      gamman += lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) * 
        log(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i)) 
      - log(tgamma(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)))) 
      +(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) - 1) * log(theta((i+n_sources))) - 
      lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i) * theta((i+n_sources));
    }
    
    
    
    
    double x = thetanorm + gamman;
    
    return x;
   

      
    
}