# Script to run a fixed form VB model on some stable isotope data with 2 isotopes
# Note that this doesn't include much maths as it's explained elsewhere
# The rough structure of this file is:
# 1. Simulate some data or use a real data set
# 2. Run JAGS on the data
# 3. Run FFVB on the data
# 4. Perform some analysis comparing the results

# Set up
rm(list = ls())
library(R2jags)
library(tidyverse)
library(mvnfast)


# Or use simmr data -------------------------------------------------------
# cons <- read.csv("wolves_consumer.csv")[c(1:21),]
# disc <- read.csv("wolves_discrimination.csv")
# sources <- read.csv("wolves_sources.csv")[c(1,4,7),]

cons <- read.csv("mantis_consumer.csv")
disc <- read.csv("mantis_discrimination.csv")
sources <- read.csv("mantis_source.csv")[c(1,3,5,7,9,11),]

mix <- cons[,c(1:2)]
s_names <- sources[,1]
s_means <- sources[,c(3,5)]
s_sds <- sources[,c(4,6)]
c_means <- disc[,c(2,4)]
c_sds <- disc[,c(3,5)]
conc <- sources[,c(7:8)]

y <- as.data.frame(mix)
n <- nrow(y)
c_0 <- c(1, 1)
d_0 <- c(1, 1)
n_isotopes <- 2

mu_kj <- s_means 
colnames(mu_kj) <- c("Meand13C", "Meand15N")
K <- nrow(mu_kj)

# Prior parameters
fmean_0 <- c(rep(0, K))
fsd_0 <- c(rep(1, K))

# # Create an iso-space plot
ggplot() +
  geom_point(data = y, aes(x = d15N, y = d13C)) +
  geom_point(data = as.data.frame(mu_kj), aes(x = Meand15N, y = Meand13C), colour = "red", shape = 18, size = 3)

# Run JAGS ----------------------------------------------------------------

JAGS_data <- list(
  N = length(y),
  y = y,
  s_mean = mu_kj,
  c_mean = c_means,
  s_sd = s_sds,
  c_sd = c_sds,
  q = conc,
  K = K,
  J = n_isotopes,
  mu_f_mean = fmean_0,
  sigma_f_sd = fsd_0
)

# Choose which parameters to save
JAGS_parameters <- c("p", "var_y", "f")

# Run the model
JAGS_run <- jags(
  data = JAGS_data,
  parameters.to.save = JAGS_parameters,
  model.file = "~/OneDrive/Documents/PhD/2nd Year/FFVB Models/Emma_code/Emma_code/complex.jags"
)

# plot(JAGS_run)
# print(JAGS_run)


# Run VB ------------------------------------------------------------------

# Source in all the generic functions
source("FF_VB_generic_functions_correct.R")

S <- 100

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
    t(rMVNormC(S, mu = mean, U = chol_prec)),
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
# K <- 4; lambda <- c(rnorm(K+K*(K+1)/2), rgamma(n_isotopes*2, 1,1))
# theta <- sim_theta(50, lambda)

# Log of likelihood added to prior
h <- function(theta) {
  p <- exp(theta[1:K]) / (sum((exp(theta[1:K]))))
  return(sum(dnorm(y[, 1],
                   mean = sum(p * conc * (mu_kj[, 1]+c_means[,1]))/sum(p*conc),
                   sd = sqrt(sum(p^2 * conc^2 * (s_sds[, 1]^2+c_sds[,1]^2))/sum(p^2*conc^2) + 1 / theta[K + 1]),
                   log = TRUE
  )) +
    sum(dnorm(y[, 2],
              mean = sum(p * conc * (mu_kj[, 2]+c_means[,2]))/sum(p*conc),
              sd = sqrt(sum(p^2 * conc^2 * (s_sds[, 2]^2+c_sds[,2]^2))/sum(p^2*conc^2) + 1 / theta[K + 2]),
              log = TRUE
    )) +
    sum(dnorm(theta[1:K], fmean_0, fsd_0, log = TRUE)) +
    sum(dgamma(theta[(K + 1):(K + 2)], shape = c_0, rate = d_0, log = TRUE)))
}
# h(theta[1,])

log_q <- function(lambda, theta) {
  
   mean <- lambda[1:K]
  # K*(K-1) precision terms
  chol_prec <- matrix(0, nrow = K, ncol = K)
  chol_prec[upper.tri(chol_prec, diag = TRUE)] <- lambda[(K + 1):(K + (K * (K + 1)) / 2)]
  
  # This version uses chol_prec directly
  # p1 <- matrix(theta[1:K] - mean, nrow = 1) %*% t(chol_prec)
  # # log_det <- unlist(determinant(prec, logarithm = TRUE))["modulus"]
  # return(-0.5 * K * log(2 * pi) - 0.5 * sum(log(diag(chol_prec))) - 0.5 * p1%*%t(p1)
  #        + sum(dgamma(theta[(K + 1)],
  #                     shape = lambda[((K + (K * (K + 1)) / 2) + 1)],
  #                     rate = lambda[(((K + (K * (K + 1)) / 2)) + 2)],
  #                     log = TRUE
  #        )))
  
  # Alternatively use the full dmvn function
  prec <- crossprod(chol_prec)
  return(dmvn(theta[1:K], mean, solve(prec), log = TRUE) + 
           + sum(dgamma(theta[(K + 1):(K + 2)],
                        shape = lambda[((K + (K * (K + 1)) / 2) + 1):(((K + (K * (K + 1)) / 2)) + n_isotopes)],
                        rate = lambda[(((K + (K * (K + 1)) / 2)) + n_isotopes + 1):(((K + (K * (K + 1)) / 2)) + n_isotopes * 2)],
                        log = TRUE
           )))
} 

#log_q(lambda, theta[1,])

# Now run it!
lambda_out <- run_VB(lambda = c(rep(0, K), rep(1, (((K * (K + 1)) / 2) + n_isotopes * 2)))) # Starting value of lambda

#----------------------------
#ffvb rcpp
lambdastart = c(rep(0, K), rep(1, (((K * (K + 1)) / 2) + n_isotopes * 2)))
Rcpp::sourceCpp("run_VB.cpp")
lambda_out_rcpp <- run_VB_cpp(as.vector(lambdastart), as.integer(K), as.integer(n_isotopes),
                              as.matrix(cons), as.matrix(s_means), as.matrix(c_means), as.matrix(c_sds), as.matrix(s_sds), as.matrix(mix))
# Perform comparisons -----------------------------------------------------

# Get all the parameters from JAGS and VB
f_JAGS <- JAGS_run$BUGSoutput$sims.list$f
n <- nrow(f_JAGS)
all_vb <- sim_theta(n, lambda_out)
f_VB <- all_vb[,1:K]

for (i in 1:K) {
  print(data.frame(
    f = c(f_JAGS[,i], f_VB[,i]),
    Fit = c(rep('JAGS', n), c(rep("VB", n)))
  ) %>% ggplot(aes(x = f, fill = Fit)) + geom_density(alpha = 0.5) + 
    ggtitle(paste("f: Iso",i)))
}

# Try the p 
p_fun <- function(x) exp(x)/sum(exp(x))
p_JAGS <- JAGS_run$BUGSoutput$sims.list$p
p_VB <- t(apply(all_vb[,1:K], 1, p_fun))

for (i in 1:K) {
  print(data.frame(
    p = c(p_JAGS[,i], p_VB[,i]),
    Fit = c(rep('JAGS', n), c(rep("VB", n)))
  ) %>% ggplot(aes(x = p, fill = Fit)) + geom_density(alpha = 0.5) +
    ggtitle(paste("p: source",i)))
}

# Try tau
tau_JAGS <- 1/JAGS_run$BUGSoutput$sims.list$var_y
tau_VB <- all_vb[,(K+1):(K + n_isotopes)]
for (i in 1:n_isotopes) {
  print(data.frame(
    sigma = c(1/sqrt(tau_JAGS[,i]), 1/sqrt(tau_VB[,i])),
    Fit = c(rep('JAGS', n), c(rep("VB", n)))
  ) %>% ggplot(aes(x = sigma, fill = Fit)) + geom_density(alpha = 0.5) +
    ggtitle(paste("sigma: Iso",i)))
}

# Now have a look at the values of h over the JAGS and VB output
h_VB <- apply(all_vb, 1, h)
JAGS_pars <- cbind(JAGS_run$BUGSoutput$sims.list$f,
                   JAGS_run$BUGSoutput$sims.list$var_y)
h_JAGS <- apply(JAGS_pars, 1, h)

data.frame(
  h = c(h_JAGS, h_VB),
  Fit = c(rep('JAGS', n), c(rep("VB", n)))
) %>% ggplot(aes(x = h, fill = Fit)) + geom_density(alpha = 0.5) +
  ggtitle("log posterior")
# So VB is getting consistently higher values of the log posterior

# Finally have a look at the two posterior predictives match the data
mu_JAGS <- p_JAGS%*%mu_kj
mu_VB <- p_VB%*%mu_kj
ypp_JAGS <- ypp_VB <- matrix(NA, ncol = n_isotopes, nrow = n)
colnames(ypp_JAGS) <- colnames(y)
colnames(ypp_VB) <- colnames(y)
for(i in 1:n_isotopes) {
  ypp_JAGS[,i] <- rnorm(n, mean = mu_JAGS[,i], sd = 1/sqrt(tau_JAGS[,i]))
  ypp_VB[,i] <- rnorm(n, mean = mu_VB[,i], sd = 1/sqrt(tau_VB[,i]))
}

ggplot() +
  geom_point(data = as.data.frame(ypp_JAGS), aes(x = d15N, y = d13C), 
             colour = 'blue', alpha = 0.5) + 
  geom_point(data = as.data.frame(ypp_VB), aes(x = d15N, y = d13C),
             colour = 'green', alpha = 0.5) + 
  geom_point(data = y, aes(x = d15N, y = d13C)) +
  geom_point(data = as.data.frame(mu_kj), aes(x = Meand15N, y = Meand13C), colour = "red", shape = 18, size = 3)

