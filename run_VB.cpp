#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
NumericMatrix rMVNormC(const double n,
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
  NumericMatrix tchol_prec = transpose(chol_prec);
  
  // 
  // NumericMatrix rMVNorm(theta.nrow(), n_sources);
  // 
  // // now need to make rMVNormC(S, mu = mean, U = chol_preg)
  // double p = (mean.length());
  // NumericMatrix Z(p, S);
  // for(int i = 0; i<(p*S); i++){
  //   Z(_,i) = Rcpp::rnorm(p);
  //   }
  // 
  // x = U^-1 Z
  
  
  
  for (int i = 0; i<n_sources; i++){
    
    // stop("Hello!");
    theta(_,i) = rMVNormC(S, mean, tchol_prec);
    
  }
  NumericMatrix gammam(S, n_tracers);
  
  for(int i = 0; i<n_tracers; i++){
    gammam(_,i) = Rcpp::rgamma(S,  lambda((n_sources + (n_sources * (n_sources + 1)) / 2) + i),
           1/lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i));
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
    
    NumericMatrix pxconc1(p.length(), concentrationmeans.ncol());
    for(int i=0; i<n_sources; i++){
      for(int j = 0; j<n_sources; j++){
        pxconc1(i,j) = p(i)*concentrationmeans(i,j);
      }}
    
    // Calculate numerator and denominator of mu
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_sources; j++){
        mutop1 += pxconc1(i,j) * (sourcemeans(i,0) + correctionmeans(i,0));
        mubtm1 += p(i) * concentrationmeans(i,j);
      }
    }
    
    // Same for sigma
    for(int i=0; i<n_sources; i++){ // sds
      for(int j =0; j<n_sources; j++){
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
      for(int j =0; j<n_sources; j++){
        mutop2 += p(i) * concentrationmeans(i,j) * (sourcemeans(i,1) + correctionmeans(i,1));
        mubtm2 += p(i) * concentrationmeans(i,j);
      }
    }
    
    
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_sources; j++){
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
  
  NumericMatrix choldecomp(n_sources, n_sources);
  int count = 0;
  for(int j = 0; j< n_sources; j++){ 
    for(int i = 0; i<n_sources; i++){
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
  
  
  NumericMatrix pt = transpose(p1);
  
  double p1tp1 = 0;
  for (int i = 0; i < p1.nrow(); i++) 
  {
    for (int j = 0; j < pt.ncol(); j++) 
      for (int k = 0; k < pt.nrow(); k++){
        {
          p1tp1 += p1(i,k) * pt(k,j) ;
        }
      }
  }
  
  
  
  
  double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) * 
      log(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i)) 
    - log(tgamma(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)))) 
    +(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) - 1) * log(theta((i+n_sources))) - 
    lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i) * theta((i+n_sources));
  }
  
  double sumlogdiag = 0;
  for(int i = 0; i<n_sources; i++){
    for(int j = 0; j<n_sources; j++){
      if(i == j){
        sumlogdiag += log(choldecomp(i,j));
      }
    }}
  
  
  double x = -0.5 * n_sources * log(2 * M_PI) - 0.5 * sumlogdiag 
    - 0.5 * p1tp1 + gamman;
  
  return x
    ;
  
}


// [[Rcpp::export]]
NumericVector delta_lqltcpp(NumericVector lambda, NumericVector theta, 
                            double eps, int n_sources, int n_tracers) {
  eps = 0.001;
  double k = lambda.length();
  NumericVector ans(k);
  NumericVector d(k);
  NumericVector lambdaplusd(k);
  NumericVector lambdaminusd(k);
  
  
  
  for(int i = 0; i<k; i++){
    
    for (int j = 0; j<k; j++){
      d(j) = 0;
    }
    d(i) = eps;
    
    
    for (int j = 0; j<k; j++){
      lambdaplusd(j) = lambda(j) + d(j);
      lambdaminusd(j) = lambda(j) - d(j);
    }
    ans(i) = (log_q_cpp(theta, lambdaplusd, n_sources, n_tracers) -  
      log_q_cpp(theta, lambdaminusd, n_sources, n_tracers))/(2 * eps);
  }
  return  ans;
}


// [[Rcpp::export]]
double h_lambdacpp(int n_sources, int n_isotopes,
                   NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                   NumericMatrix correctionmeans,
                   NumericMatrix corrsds, NumericMatrix sourcesds,
                   NumericVector theta, NumericMatrix y,
                   NumericVector lambda) {
  
  
  
  
  
  
  
  
  return hcpp(n_sources, n_isotopes, concentrationmeans, sourcemeans, correctionmeans,
              corrsds, sourcesds, theta, y) - log_q_cpp(theta, lambda, n_sources, n_isotopes);
}


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

// [[Rcpp::export]]
NumericVector nabla_LB_cpp(NumericVector lambda, NumericMatrix theta, int n_sources, int n_tracers,
                           NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                           NumericMatrix correctionmeans,
                           NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix y,
                           NumericVector c){
  
  NumericMatrix big_c(theta.nrow(), c.length());
  
  //working
  NumericMatrix big_delta_lqlt(theta.nrow(), lambda.length()); 
  NumericMatrix big_h_lambda_rep(lambda.length(), theta.nrow());
  NumericMatrix big_h_lambda_rep_transpose(theta.nrow(), lambda.length());
  //NumericMatrix big_delta_lqlt_transpose(lambda.length(), theta.nrow());
  
  //not working - size is right but numbers wrong
  NumericVector big_h_lambda(theta.nrow());
  NumericVector big_h_lambda_transpose(theta.nrow());
  
  
  for(int i = 0; i <theta.nrow(); i++){
    big_delta_lqlt(i,_) = delta_lqltcpp(lambda, theta(i,_), 0.01, n_sources, n_tracers);
  }
  
  for(int i =0; i<theta.nrow(); i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta(i,_), y,
                 lambda);
  }
  
  //big_delta_lqlt_transpose = transpose(big_delta_lqlt);
  
  for(int i =0; i<lambda.length(); i++){
    big_h_lambda_rep(i,_) = big_h_lambda;
  }
  big_h_lambda_rep_transpose = transpose(big_h_lambda_rep);
  
  // for(int i = 0; i<lambda.length(); i++){
  //   c(i) = 0;
  //  }
  
  for(int i =0; i<theta.nrow(); i++){
    big_c(i,_) = c;
  }
  
  NumericMatrix big_h_minus_c(theta.nrow(), lambda.length());
  //NumericMatrix big_h_minus_c_t(lambda.length(), theta.nrow());
  
  for (int i = 0; i<theta.nrow(); i++){
    for(int j = 0; j<lambda.length(); j++){
      big_h_minus_c(i,j) = big_h_lambda_rep_transpose(i,j) - c(j);
    }
  }
  
  //big_h_minus_c_t = transpose(big_h_minus_c);
  
  NumericMatrix ansmat(big_delta_lqlt.nrow(), big_h_minus_c.ncol());
  
  for (int i = 0; i < big_delta_lqlt.nrow(); i++) 
  {
    for (int j = 0; j < big_delta_lqlt.ncol(); j++) {
      
      
      ansmat(i,j) = big_delta_lqlt(i,j) * big_h_minus_c(i,j);
      
      
    }
  }
  
  NumericVector ans(ansmat.ncol());
  for(int i = 0; i<ansmat.ncol(); i++){
    
    ans(i) = mean(ansmat(_,i));
    
  }
  
  return ans;
}







// [[Rcpp::export]]
NumericVector control_var_cpp(NumericVector lambda, NumericMatrix theta, int n_sources, int n_tracers,
                              NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                              NumericMatrix correctionmeans,
                              NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix y){
  NumericMatrix big_delta_lqlt(theta.nrow(), lambda.length()); 
  NumericMatrix big_h_lambda_rep(lambda.length(), theta.nrow());
  NumericMatrix big_h_lambda_rep_transpose(theta.nrow(), lambda.length());
  NumericVector big_h_lambda(theta.nrow());
  NumericVector big_h_lambda_transpose(theta.nrow());
  
  for(int i = 0; i <theta.nrow(); i++){
    big_delta_lqlt(i,_) = delta_lqltcpp(lambda, theta(i,_), 0.01, n_sources, n_tracers);
  }
  
  for(int i =0; i<theta.nrow(); i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta(i,_), y,
                 lambda);
  }
  
  
  
  for(int i =0; i<lambda.length(); i++){
    big_h_lambda_rep(i,_) = big_h_lambda;
  }
  big_h_lambda_rep_transpose = transpose(big_h_lambda_rep);
  
  NumericMatrix big_nabla(big_delta_lqlt.nrow(), big_h_lambda_rep_transpose.ncol());
  
  for (int i = 0; i < big_delta_lqlt.nrow(); i++) 
  {
    for (int j = 0; j < big_delta_lqlt.ncol(); j++) {
      
      
      big_nabla(i,j) = big_delta_lqlt(i,j) * big_h_lambda_rep_transpose(i,j);
      
      
    }
  }
  
  NumericVector var_big_delta_lqlt(big_delta_lqlt.ncol());
  
  for(int i = 0; i<big_delta_lqlt.ncol(); i++){
    var_big_delta_lqlt(i) = var(big_delta_lqlt(_,i));
  }
  
  NumericMatrix covmat(big_nabla.ncol(), big_delta_lqlt.ncol());
  
  covmat = cov_mat_cpp(big_nabla, big_delta_lqlt);
  
  NumericVector diag(covmat.ncol());
  for(int i =0; i<covmat.ncol(); i++){
    for(int j =0; j<covmat.ncol(); j++){
      if(i == j){
        diag(i) = covmat(i,j);
      }
    }}
  
  NumericVector ans(diag.length());
  for(int i =0; i<diag.length(); i++){
    ans(i) = diag(i)/var_big_delta_lqlt(i);
  }
  
  //divide
  return ans;
}

// [[Rcpp::export]]
double LB_lambda_cpp(NumericMatrix theta, NumericVector lambda, NumericVector p, int n_sources, int n_isotopes,
                     NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                     NumericMatrix correctionmeans,
                     NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix y){
  
  NumericVector hlambdaapply(theta.nrow());
  
  for(int i = 0; i <theta.nrow(); i++){
    hlambdaapply(i) = h_lambdacpp(n_sources, n_isotopes, concentrationmeans, sourcemeans,
                 correctionmeans, corrsds, sourcesds, theta(i,_), y, lambda);
  }
  
  double ans = mean(hlambdaapply);
  
  return ans;
  
  
}



// [[Rcpp::export]]
List run_VB_cpp(NumericVector lambdastart,
                         int n_sources,
                         int n_tracers,
                         NumericMatrix concentrationmeans,
                         NumericMatrix sourcemeans,
                         NumericMatrix correctionmeans,
                         NumericMatrix corrsds,
                         NumericMatrix sourcesds,
                         NumericMatrix y
                             ){

  // starting values
  int S = 100;
    int P = 10;
    double beta_1 = 0.9;
    double beta_2 = 0.9;
    int tau = 1000;
    double eps_0 = 0.1;
    int t_W = 49;

  NumericMatrix theta(S, (n_sources + n_tracers));

  theta = sim_thetacpp(S, lambdastart, n_sources, n_tracers);

  NumericVector c(lambdastart.length());
  c = control_var_cpp(lambdastart, theta, n_sources, n_tracers, 
                      concentrationmeans,
                      sourcemeans, correctionmeans,
                      corrsds, sourcesds, y);

  NumericVector g_0(lambdastart.length());
  g_0 = nabla_LB_cpp(lambdastart, theta, 
                     n_sources, n_tracers, 
                     concentrationmeans, sourcemeans,
                     correctionmeans, corrsds, 
                     sourcesds, y, c);

  NumericVector nu_0(lambdastart.length());
  
  for(int i = 0; i<lambdastart.length(); i++){
    nu_0(i) = pow(g_0(i),2);
    }


  NumericVector g_bar(lambdastart.length());
  g_bar = g_0;

  NumericVector nu_bar(lambdastart.length());
  nu_bar = nu_0;

  NumericVector g_t(lambdastart.length());
  NumericVector nu_t(lambdastart.length());

  double patience = 0;
  bool stop = FALSE;
  double max_LB_bar = R_NegInf;
  double alpha_t = 0;
  double t = 0;
  NumericVector LB(t_W+1);

  for(int i=0; i<t_W; i++){
    LB(i) = 0;
    }
  NumericVector lambda = lambdastart;

  while(stop == FALSE){

  theta = sim_thetacpp(S, lambda, n_sources, n_tracers);

  g_t = nabla_LB_cpp(lambda, theta, n_sources, n_tracers, concentrationmeans,
                     sourcemeans, correctionmeans, corrsds, sourcesds,
                     y, c);

  c = control_var_cpp(lambda, theta,n_sources,n_tracers,
                      concentrationmeans, sourcemeans,
                      correctionmeans,
                      corrsds,sourcesds, y);

  for(int i=0; i<nu_t.length(); i++){
    nu_t(i) = pow(g_t(i),2);
    }

  for(int i=0; i<g_bar.length(); i++){
    g_bar(i) = beta_1 * g_bar(i) + (1-beta_1) * g_t(i);
    nu_bar(i) = beta_2 * nu_bar(i) + (1-beta_2) * nu_t(i);
    }

  NumericVector alpha_min(2);

  for(int i=0; i<2; i++){
    alpha_min(0) = eps_0;
    alpha_min(1) = eps_0 * tau/(t+1);
    }

  alpha_t = Rcpp::min(alpha_min);



 // # Update lambda
  for(int i = 0; i<lambda.length(); i++){
    lambda(i) = lambda(i) + alpha_t * g_bar(i)/pow(nu_bar(i), 0.5);
    }
  // # Compute the moving average LB if out of warm-up
  if(t<=t_W){

  // # Compute a new lower bound estimate
  NumericVector p = hfn(theta, n_sources);

  LB(t) = LB_lambda_cpp(theta,lambda, p, n_sources, n_tracers,
    concentrationmeans, sourcemeans,
    correctionmeans,
    corrsds,sourcesds, y);
  }
  else{
    for (int i = 0; i<(t_W-1); i++){
      LB(i) = LB(i+1);
      }

    NumericVector p = hfn(theta, n_sources);

    LB(t_W) = LB_lambda_cpp(theta, lambda, p, n_sources, n_tracers,
       concentrationmeans, sourcemeans,
       correctionmeans,
       corrsds,sourcesds, y);


  double LB_bar = mean(LB);

  NumericVector maxbar(2);
  for(int i=0; i<2; i++){
    maxbar(0) = max_LB_bar;
    maxbar(1) = LB_bar;
  }

  max_LB_bar = Rcpp::max(maxbar);

    if(LB_bar>= max_LB_bar){
      patience = 0;
    } else{
      patience = patience +1;
      }

    }

//////////////////////
    if(patience>P){
      stop = TRUE;
      }
    t = t + 1;
    }

  return Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("patience") = patience,
                            Rcpp::Named("t_W") = t_W,
                            Rcpp::Named("alpha_t") = alpha_t,
                            Rcpp::Named("LB") = LB,
                            Rcpp::Named("lambdastart") = lambdastart,
                            Rcpp::Named("gbar") = g_bar,
                            Rcpp::Named("g0") = g_0,
                            Rcpp::Named("nu0") = nu_0,
                            Rcpp::Named("nubar") = nu_bar
                              );
  }

// alpha_t <- min(eps_0, eps_0 * tau / t)
//
// # Update lambda
//   lambda <- lambda + alpha_t * g_bar / sqrt(nu_bar)
//
// # Compute the moving average LB if out of warm-up


