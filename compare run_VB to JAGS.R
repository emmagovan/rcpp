library(R2jags)
library(tidyverse)
library(Rcpp)

#use simmr data -------------------------------------------------------

mix <- matrix(c(
  -10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54,
  -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89,
  11.73, 10.89, 11.05, 12.3
), ncol = 2, nrow = 10)
colnames(mix) <- c("d13C", "d15N")
s_names <- c("Zostera", "Grass", "U.lactuca", "Enteromorpha")
s_means <- matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol = 2, nrow = 4)
s_sds <- matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol = 2, nrow = 4)
c_means <- matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol = 2, nrow = 4)
c_sds <- matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol = 2, nrow = 4)
conc <- matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol = 2, nrow = 4)

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
sourceCpp("run_VB.cpp")
lambdastart = c(rep(0, K), rep(1, (((K * (K + 1)) / 2) + n_isotopes * 2)))
#lambdastart<-c(0,0,0,0,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,0,1)
res<-run_VB_cpp(lambdastart, 4,2, conc, s_means, c_means, c_sds,
                s_sds, mix)

#Check time --------------------------------------------------------------
library(simmr)
simmr_in = simmr_load(mixtures=mix,
                      source_names=s_names,
                      source_means=s_means,
                      source_sds=s_sds,
                      correction_means=c_means,
                      correction_sds=c_sds,
                      concentration_means = conc)

library(microbenchmark)
microbenchmark(simmr_out = simmr_mcmc(simmr_in),
VB = run_VB_cpp(lambdastart, 4,2, conc, s_means, c_means, c_sds,
           s_sds, mix), times = 10L
)


# Perform comparisons -----------------------------------------------------

# Get all the parameters from JAGS and VB
f_JAGS <- JAGS_run$BUGSoutput$sims.list$f
n <- nrow(f_JAGS)
all_vb <- sim_thetacpp(n, res$lambda, K, n_isotopes)
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
    ggtitle(paste("p: Iso",i)))
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


