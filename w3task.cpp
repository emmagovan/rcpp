#include <Rcpp.h>
using namespace Rcpp;

//What we want to do is start with a vector x and specify i iterations
//We get the mean and sd of x
//We use those to generate a new x
//Repeat i times, return vector of mu, vector of sd, and last x



//[[Rcpp::export]]
List my_fn(NumericVector xvec, int iter){
  //This has created a way to store mu and sd
  //x is what we give the function to start
  //we can just overwrite x every time
  NumericVector muvec(iter);
  NumericVector sdvec(iter);
  int xs = xvec.length();
  
  NumericVector sdnum(xs);
  double sdnum2;
  double sdnum3;
  
 
  int i;
  int iteration = iter;
  
for(i = 0; i<iteration; i++){
  double xsum = 0;
    //this gets the sum of xvalues
    for(int j = 0; j<xs; j++){
     xsum += xvec(j);
      //return xsum;
      }

    //This calculates the mean
    muvec(i) = xsum/xs;

  //   Now get sd

  for (int j = 0; j<xs; j++){
    sdnum(j) = pow((xvec(j) - muvec(i)), 2);
    
    }
  
  sdnum2 = 0;
  for (int j=0; j<xs; j++){
    sdnum2 += sdnum(j);
    }
  
  sdnum3 = 0;
  sdnum3= sdnum2/(xs-1);
  sdvec(i) = pow(sdnum3, 0.5);
  
  double currmu = muvec(i);
  double currsd = sdvec(i);
  
  for(int j = 0; j<xs; j++){
    xvec(j) = 0;
    xvec(j) = R::rnorm(currmu, currsd);
  }

}

   
return Rcpp::List::create(Rcpp::Named("mean") = muvec,
                          Rcpp::Named("sd") = sdvec,
                          Rcpp::Named("sdnum") = sdnum,
                          Rcpp::Named("sdnum2") = sdnum2,
                          Rcpp::Named("sdnum3") = sdnum3,
                          Rcpp::Named("xvector") = xvec);

  

  }




// So far this returns the mean for x every time
// Now need it to return sd for x every time
// Then create new x for each loop




