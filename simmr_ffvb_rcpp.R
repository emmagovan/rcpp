simmr_ffvb<-function(simmr_in,
                     prior_control = list(
                       means = rep(
                         0,
                         simmr_in$n_sources),
                       sd = rep(
                         1,
                         simmr_in$n_sources),
                       c_0 = rep(
                         1,
                         simmr_in$n_tracers),
                       d_0 = rep(
                         1,
                         simmr_in$n_tracers))
){
  
  # Throw a warning if less than 4 observations in a group - 1 is ok as it wil do a solo run
  if (min(table(simmr_in$group)) > 1 & min(table(simmr_in$group)) < 4) warning("At least 1 group has less than 4 observations - either put each observation in an individual group or use informative prior information")
  
  
  output <- vector("list", length = simmr_in$n_groups)
  names(output) <- levels(simmr_in$group)
  K<-simmr_in$n_sources
  n_tracers <- simmr_in$n_tracers
  n_output<-3600
  lambdares<-matrix(rep(NA, ((( K + (K * (K + 1)) / 2)) + n_tracers * 2) * simmr_in$n_groups),
                    nrow =((( K + (K * (K + 1)) / 2)) + n_tracers * 2),
                    ncol = simmr_in$n_groups)
  thetares<-matrix(rep(NA, ((K+n_tracers) * n_output*simmr_in$n_groups)),
                   ncol = (K+n_tracers),
                   nrow = n_output*simmr_in$n_groups)
  
  # Loop through all the groups
  for (i in 1:simmr_in$n_groups) {
    if (simmr_in$n_groups > 1) cat(paste("\nRunning for group", levels(simmr_in$group)[i], "\n\n"))
    
    curr_rows <- which(simmr_in$group_int == i)
    curr_mix <- simmr_in$mixtures[curr_rows, , drop = FALSE]
    
    # Determine if a single observation or not
    if (nrow(curr_mix) == 1) {
      cat("Only 1 mixture value, performing a simmr solo run...\n")
      solo <- TRUE
    } else {
      solo <- FALSE
    }
    
    
    n_tracers <- simmr_in$n_tracers
    n_sources <- simmr_in$n_sources
    s_names <- simmr_in$source_names
    K<-simmr_in$n_sources
    S <- 100
    source_means = simmr_in$source_means
    source_sds = simmr_in$source_sds
    correction_means = simmr_in$correction_means
    correction_sds = simmr_in$correction_sds
    concentration_means = simmr_in$concentration_means
    y = curr_mix
    
    Rcpp::sourceCpp("run_VB.cpp")

    lambdastart = c(rep(0, K), rep(1, (((K * (K + 1)) / 2) + n_tracers * 2)))
   
    lambdares[,i]<-run_VB_cpp(lambdastart, K, n_tracers, concentration_means, 
                          source_means, correction_means, correction_sds,
                    source_sds, y)
    
    thetares[(1+3600*(i-1)):(3600*i),] = sim_thetacpp(n_output, lambdares[,i], K, n_tracers)
  }
  

  mylist<-list(lambda = lambdares,
               n_sources = simmr_in$n_sources, 
               n_tracers = simmr_in$n_tracers,
               group = simmr_in$n_groups, 
               source_names = simmr_in$source_names,
               theta = thetares
  )
  
  class(mylist) <- "simmr_ffvb_output"
  return(mylist)
}

#Check
#simmr_out1 = simmr_ffvb(simmr_in)










