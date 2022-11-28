#compare time for rcpp and r

library(microbenchmark)
library(simmr)

#read in data

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



lambda = c(0,0,0,0,rep(1,14))


microbenchmark(rcpp = sim_thetacpp(100, lambda, 4, 2),
              r = sim_theta(100, lambda))

# rcpp much faster


theta<-sim_thetacpp(100, lambdastart, 4, 2)

microbenchmark(rcpp = hcpp(4,2, simmr_in$concentration_means, 
                           simmr_in$source_means, 
                           simmr_in$correction_means, 
                           simmr_in$correction_sds, 
                           simmr_in$source_sds, theta[1,], 
                           simmr_in$mixtures),
               r = h(theta[1,]))
#rcpp much faster again

microbenchmark(rcpp = log_q_cpp(theta[1,], lambda, 4, 2),
               r = log_q(lambda, theta[1,]))

microbenchmark(rcpp = delta_lqltcpp(lambda, theta[1,], 0.001, 4, 2),
               r = delta_lqlt(lambda, theta[1,]))


microbenchmark(rcpp = h_lambdacpp(4,2, simmr_in$concentration_means, simmr_in$source_means, simmr_in$correction_means, simmr_in$correction_sds,
                                  simmr_in$source_sds, theta[1,], simmr_in$mixtures, lambda),
               r = h_lambda(lambda, theta[1,], simmr_in$mixtures))



microbenchmark(rcpp = cov_mat_cpp(simmr_in$source_sds, simmr_in$source_means),
               r = cov(simmr_in$source_sds, simmr_in$source_means))

microbenchmark(rcpp = nabla_LB_cpp(lambda, theta, 4,2, simmr_in$concentration_means, simmr_in$source_means, simmr_in$correction_means, simmr_in$correction_sds,
                                   simmr_in$source_sds, simmr_in$mixtures, c = lambda),
               r = nabla_LB(lambda, theta, lambda))

microbenchmark(rcpp = control_var_cpp(lambda, theta, 4, 2, simmr_in$concentration_means, simmr_in$source_means, simmr_in$correction_means, simmr_in$correction_sds,
                                      simmr_in$source_sds, simmr_in$mixtures),
               r = control_var(lambda, theta))

microbenchmark(rcpp = LB_lambda_cpp(theta, lambda, hfn(theta, 4),4, 2, simmr_in$concentration_means, simmr_in$source_means, simmr_in$correction_means, simmr_in$correction_sds,
                                    simmr_in$source_sds, simmr_in$mixtures),
               r = LB_lambda(lambda, theta))

microbenchmark(rcpp = run_VB_cpp(lambdastart = c(0,0,0,0, rep(1,14)), simmr_in$n_sources,
                                 simmr_in$n_tracers,
                                 simmr_in$concentration_means,
                                 simmr_in$source_means,
                                 simmr_in$correction_means,
                                 simmr_in$correction_sds,
                                 simmr_in$source_sds,
                                 simmr_in$mixtures),
               run_VB(lambda))









