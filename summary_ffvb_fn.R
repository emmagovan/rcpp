summary_ffvb<-function(simmr_out, type = c("quantiles", "statistics", "correlations"), group = 1, ...){
  
  # Get the specified type
  type <- match.arg(type, several.ok = TRUE)
  K = simmr_out$n_sources
  n_tracers = simmr_out$n_tracers
  lambda = simmr_out$lambda
  n_groups = simmr_out$group
  names = simmr_out$source_names
  theta = simmr_out$theta
  
  
  
  
  for (i in 1:n_groups) {
    cat(paste("\nSummary for group", i))
    
    all_VB <- theta# Starting value of lambda
    
    p_fun <- function(x) exp(x)/sum(exp(x))
    p_VB <- t(apply(all_VB[,1:K], 1, p_fun))
    out_all <- matrix(cbind(p_VB, 1/sqrt(all_VB[,(K+1):(K+n_tracers)])), ncol = (K+n_tracers))
    
    out_quantiles <- array(rep(NA, 5*ncol(all_VB)*n_groups), dim = c(ncol(all_VB), 5, n_groups))
    #5 because 5 different quantiles taken
    
    colnames(out_quantiles) = c("0.025", "0.25", "0.5", "0.75", "0.975")
    rownames(out_quantiles)=c(simmr_out$source_names, rep("sd", n_tracers))
    
    
    out_statistics <- array(rep(NA,n_groups*2*ncol(all_VB)), dim = c(ncol(all_VB), 2, n_groups))
    colnames(out_statistics) = c("mean", "sd")
    rownames(out_statistics)=c(simmr_out$source_names, rep("sd", n_tracers))
    
    
    out_cor <- matrix(rep(NA,n_groups^2), ncol = n_groups, nrow = n_groups)
    
    out_quantiles[,,i] <- t(apply(out_all, 2, "quantile", probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    #  coda:::summary.mcmc.list(object$output)$quantiles
    out_statistics[,,i] <- t(apply(out_all, 2, function(x) {
      return(c(mean = mean(x), sd = stats::sd(x)))
    }))
    # coda:::summary.mcmc.list(object$output)$statistics[,1:2]
    out_cor <- stats::cor(out_all)
    rownames(out_cor)=c(simmr_out$source_names, rep("sd", n_tracers))
    colnames(out_cor)=c(simmr_out$source_names, rep("sd", n_tracers))
    
    
    if ("quantiles" %in% type) {
      # Print out quantiles argument
      cat(paste0("\nquantiles\n"))
      print(round(out_quantiles[,,i], 3))
    }
    
    if ("statistics" %in% type) {
      cat(paste0("\nstatistics\n"))
      
      # Print out quantiles argument
      print(round(out_statistics[,,i], 3))
    }
    
    
  }
  
  if ("correlations" %in% type) {
    cat(paste0("\ncorrelations\n"))
    # Print out quantiles argument
    print(round(out_cor, 3))
  }
}
