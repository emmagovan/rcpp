# Fit a 1 isotope model using VB

# Clear workspace and call in packages
rm(list = ls())
library(mvnfast)

# Call in the external VB functions
source("FF_VB_generic_functions_correct.R")


# Simulate some data
col1<-rnorm(10, mean = 8, sd =0.2)
col2<-rnorm(10, mean = 8, sd=0.2)
y<-matrix(c(col1, col2), ncol = 2, byrow = FALSE)
colnames(y)<-c("d15N", "d13C")
y<-as.data.frame(y)
mu_kj <- matrix(c(1,1,10,1,10,10), nrow = 3)
colnames(mu_kj)<-c("Meand15N", "Meand13C")
K<-nrow(mu_kj)

# Prior distributions
n <- nrow(y)
c_0 <- c(1, 1)
d_0 <- c(1, 1)
n_isotopes <- 2

# Prior for p <- exp(f)/sum(exp(f))
fmean_0 <- c(rep(0, K))
fsd_0 <- diag(K)


# # Create an iso-space plot
ggplot() +
  geom_point(data = y, aes(x = d15N, y = d13C)) +
  geom_point(data = as.data.frame(mu_kj), aes(x = Meand15N, y = Meand13C), colour = "red", shape = 18, size = 3)

#--------------------------------

# Parameters for VB algorithm
S <- 100
lambda_start <- c(rep(0, K), rep(1, 2 + K + K*(K-1)/2)) # Don't think it's K^2

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
                   mean = sum(p * (mu_kj[, 1]))/sum(p),
                   sd = sqrt(sum(p^2 * (s_sds[, 1]^2+c_sds[,1]^2))/sum(p^2*conc^2) + 1 / theta[K + 1]),
                   log = TRUE
  )) +
    sum(dnorm(y[, 2],
              mean = sum(p * (mu_kj[, 2]))/sum(p),
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

# Now run it!
lambda_out <- run_VB(lambda = c(rep(0, K), rep(1, (((K * (K + 1)) / 2) + n_isotopes * 2)))) # Starting value of lambda


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
