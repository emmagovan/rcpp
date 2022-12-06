# Fit a 1 isotope model using VB

# Clear workspace and call in packages
rm(list = ls())
library(mvnfast)

# Call in the external VB functions
source("FF_VB_generic_functions_correct.R")

# Simulate some data
set.seed(123)
N <- 30
K <- 3
mu_s <- seq(-10, 10, length = K)
sigma_s <- rep(1, K)
s <- rnorm(K, mu_s, sigma_s)
tau <- 1
p <- c(0.8, 0.15, 0.05)
y <- rnorm(N, sum(p * s), 1 / sqrt(tau))

# Prior distributions

# Prior for tau
c_0 <- 1
d_0 <- 1

# Prior for p <- exp(f)/sum(exp(f))
fmean_0 <- c(rep(0, K))
fsd_0 <- diag(K)

# Parameters for VB algorithm
S <- 100
lambda_start <- c(rep(0, K), rep(1, 2 + K + K*(K-1)/2)) # Don't think it's K^2

# Sim_theta function
sim_theta <- function(S, lambda) {
  mean <- lambda[1:K]
  # K*(K-1) precision terms
  chol_prec <- matrix(0, nrow = K, ncol = K)
  chol_prec[upper.tri(chol_prec, diag = TRUE)] <- lambda[(K + 1):(K + (K * (K + 1)) / 2)]

  # This is a placeholder for more advanced code using chol_prec directly
  theta <- cbind(
    t(rMVNormC(S, mu = mean, U = chol_prec)),
    rgamma(S,
      shape = lambda[((K + (K * (K + 1)) / 2) + 1)],
      rate = lambda[(((K + (K * (K + 1)) / 2)) + 2)]
    )
  )
}
# theta <- sim_theta(S, lambda_start)

# Log of likelihood added to prior
h <- function(theta) {
  p <- exp(theta[1:K]) / sum(exp(theta[1:K]))
  sum(dnorm(y,
    mean = sum(p * s),
    sd = sqrt(sum((p^2) * (sigma_s^2)) / (sum(p^2)) + (1 / theta[(K + 1)])), log = TRUE
  )) +
    dmvn(theta[1:K], fmean_0, fsd_0, log = TRUE) +
    dgamma(theta[(K + 1)], shape = c_0, rate = d_0, log = TRUE)
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
           sum(dgamma(theta[(K + 1)],
                      shape = lambda[((K + (K * (K + 1)) / 2) + 1)],
                      rate = lambda[(((K + (K * (K + 1)) / 2)) + 2)],
                      log = TRUE
           )))
  
}
# log_q(lambda_start, theta[1,])

# Run the model
lambda <- run_VB(lambda = lambda_start)

# Check the results -------------------------------------------------------

# Generate theta
theta <- sim_theta(10000, lambda)

# Convert to p
fcalc <- function(x) exp(x) / sum(exp(x))
prop <- t(apply(theta[, 1:K], 1, fcalc))

# Create plots
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, .7, 0), tck = -0.01, las = 1)
hist(prop[, 1], breaks = 30)
abline(v = p[1])
hist(prop[, 2], breaks = 30)
abline(v = p[2])
hist(prop[, 3], breaks = 30)
abline(v = p[3])
hist(theta[, 4], breaks = 30)
abline(v = tau)
par(mfrow = c(1, 1))
