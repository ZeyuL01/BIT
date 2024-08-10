#Show the demo

#' Demo to show BIMTR, which can be used without loading the ChIP-seq data.
#'
#' @param N number of iterations for gibbs sampler, recommended for >= 5000.
#' @param dat pre-processed Demo data, including CTCF, ZBTB7A. For details please check the manual for the data.
#'
#' @return A list object contains results of Gibbs sampler.
#' @export
Demo = function(N,dat=c("CTCF","ZBTB7A")){
  if(dat=="CTCF"){
    Demo_dat = CTCF_Demo
  }else if(dat=="ZBTB7A"){
    Demo_dat = ZBTB7A_Demo
  }

  tf_labs <- as.numeric(factor(Demo_dat$TF))
  xct <- Demo_dat$GOOD
  nct <- Demo_dat$TOTAL

  start_time <- Sys.time()

  results <- Main_Sampling(N,xct, nct, tf_labs)

  end_time <- Sys.time()
  time_used <- difftime(end_time,start_time)
  print(paste0("Time used: ", round(time_used,2)," ", units(time_used)," for ", N," rounds."))

  results[["TF_names"]] <- Demo_dat$TF
  return(results)
}


#' gelman-rubin demo
#'
#' @description Gelman-Rubin convergence diagnostic for Demo dataset, recommend for N >= 5000
#' @param N number of rounds
#' @param dat pre-processed Demo data, including CTCF, ZBTB7A. For details please check the manual for the data.
#' @param nchains number of chains.
#'
#' @return A list object contains gelman-rubin diagnostic results.
#' @export
gelman_rubin_Demo <- function(N = 5000 , dat=c("CTCF","ZBTB7A"), nchains = parallel::detectCores()-1){
  if(nchains > parallel::detectCores()){
    stop("number of chains computed should be less or equal to number of cores.")
  }
  local_cl <- parallel::makeCluster(nchains)
  parallel::clusterExport(cl=local_cl, c("N","dat"),envir = environment())
  results <- parallel::clusterEvalQ(local_cl, {
    library(BIMTR)
    Demo(N = N, dat = dat)
  })
  parallel::stopCluster(local_cl)

  mu0_nchains_list <- mcmc_list_transformer(results,"mu0",nchains)
  tau0S_nchains_list <- mcmc_list_transformer(results,"tau0S",nchains)
  tau1S_nchains_list <- mcmc_list_transformer(results,"tau1S",nchains)
  theta_ij_nchains_list <- mcmc_list_transformer(results,"theta_ij",nchains)
  theta_i_nchains_list <- mcmc_list_transformer(results,"theta_i",nchains)
  sigmaS_nchains_list <- mcmc_list_transformer(results,"sigmaS",nchains)

  results_list<-list()
  results_list[["mu0"]] <- coda::gelman.diag(mu0_nchains_list)
  results_list[["tau0S"]] <- coda::gelman.diag(tau0S_nchains_list)
  results_list[["tau1S"]] <- coda::gelman.diag(tau1S_nchains_list)
  results_list[["theta_ij"]] <- list()
  for(j in 1:length(theta_ij_nchains_list)){
    results_list[["theta_ij"]][[j]] <- coda::gelman.diag(theta_ij_nchains_list[[j]])
  }

  results_list[["theta_i"]] <- list()
  for(j in 1:length(theta_i_nchains_list)){
    results_list[["theta_i"]][[j]] <- coda::gelman.diag(theta_i_nchains_list[[j]])
  }

  results_list[["sigmaS"]] <- list()
  for(j in 1:length(sigmaS_nchains_list)){
    results_list[["sigmaS"]][[j]] <- coda::gelman.diag(sigmaS_nchains_list[[j]])
  }

  return(results_list)
}
