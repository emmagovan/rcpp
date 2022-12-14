# See if we can re-create the crash

library(simmr)

# Souce in all the cpp files
Rcpp::sourceCpp("run_VB.cpp")

# Create the data required
K <- 4
lambda<-c(rep(0,4), rep(1,14))
n_isotopes <- 2

# Load in sim_theta
sim_theta <- function(S, lambda) {
  
  # For K parameters you will have
  # lambda is of length K+K*(K+1)/2 +n_isotopes*2
  # mean <- lambda[1:K]
  # chol_prec is made up of lambda[(K + 1):(K+(K*(K+1))/2)]
  # Tau is made up of lambda[((K+(K*(K+1))/2)+1):((K+(K*(K+1))/2)+n_isotopes*2)]
  # (f) ~ MVN(lambda[1:K], solve(crossprod(chol_prec)))
  
  mean <- lambda[1:K]
  # K*(K-1) precision terms
  chol_prec <- matrix(0, nrow = K, ncol = K)
  chol_prec[upper.tri(chol_prec, diag = TRUE)] <- lambda[(K + 1):(K + (K * (K + 1)) / 2)]
  
  # This is a placeholder for more advanced code using chol_prec directly
  theta <- cbind(
    rMVNormC(S, mu = mean, U = chol_prec),
    matrix(rgamma(S * n_isotopes,
                  shape = lambda[((K + (K * (K + 1)) / 2) + 1):(((K + (K * (K + 1)) / 2)) + n_isotopes)],
                  rate = lambda[(((K + (K * (K + 1)) / 2)) + n_isotopes + 1):(((K + (K * (K + 1)) / 2)) + n_isotopes * 2)]
    ),
    nrow = S,
    ncol = n_isotopes,
    byrow = TRUE
    )
  )
  
  return(theta)
}
theta<-sim_theta(100, lambda)

# Load in the simmr object
mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54, 
               -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89, 
               11.73, 10.89, 11.05, 12.3), ncol = 2, nrow = 10)
colnames(mix) = c('d13C','d15N')
s_names = c("Zostera", "Grass", "U.lactuca", "Enteromorpha")
s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol = 2, nrow = 4)
s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol = 2, nrow = 4)
c_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol = 2, nrow = 4)
c_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol = 2, nrow = 4)
conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol = 2, nrow = 4)
simmr_in = simmr_load(mixtures = mix,
                      source_names = s_names,
                      source_means = s_means,
                      source_sds = s_sds,
                      correction_means = c_means,
                      correction_sds = c_sds,
                      concentration_means = conc)

# This is the line that causes everything to crash
c <-control_var_cpp(lambda, theta, 4, 2,
                    simmr_in$concentration_means,
                    simmr_in$source_means,
                    simmr_in$correction_means,
                    simmr_in$correction_sds,
                    simmr_in$source_sds,
                    simmr_in$mixtures)

hcpp_ans <- with(simmr_in,
                 hcpp(n_sources, n_isotopes,
                      concentration_means,
                      source_means,
                      correction_means,
                      correction_sds,
                      source_sds,
                      theta[1,],
                      mixtures))



