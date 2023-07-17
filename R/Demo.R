#Show the demo

#' Demo to show BayesIMTR, which can be used without loading the ChIP-seq data.
#'
#' @param N number of rounds for gibbs sampler, recommended for > 2000.
#' @param dat pre-processed Demo data, including CTCF, ZBTB7A and KDM1A. For details please check the manual for the data.
#'
#' @return A list object contains results of Gibbs sampler.
#' @export
Demo = function(N,dat=c("CTCF","ZBTB7A","KDM1A")){
  if(dat=="CTCF"){
    Demo_dat = CTCF_Demo
  }else if(dat=="ZBTB7A"){
    Demo_dat = ZBTB7A_Demo
  }else if(dat=="KDM1A"){
    Demo_dat = KDM1A_Demo
  }

  tf_labs<-as.numeric(factor(Demo_dat$TF))
  xct<-Demo_dat$GOOD
  nct<-Demo_dat$TOTAL

  start_time <- Sys.time()

  results<-Main_Sampling(N,xct,nct,tf_labs)

  end_time <- Sys.time()
  time_used<-difftime(end_time,start_time)
  print(paste0("Time used: ", round(time_used,2)," ",units(time_used)," for ",N," rounds."))

  results[["TF_names"]]<-Demo_dat$TF
  return(results)
}

#' To show the ranking table by inspecting the results of Gibbs sampler.
#'
#' @param dat results from Main_sampling.
#' @param burnin number of samples used for burn-in.
#'
#' @return a list object contains two data.frame. one by ranking the theta_ij, the other one by ranking theta_i.
#' @export
Show_Results<-function(dat,burnin){
  TF_names <- dat[["TF_names"]]
  theta_ij_mat<-dat$theta_ij
  theta_i_mat<-dat$theta_i
  tf_results_ij<-rowMeans(theta_ij_mat[,(dim(theta_ij_mat)[2]-burnin):dim(theta_ij_mat)[2]])
  tf_results_i<-rowMeans(theta_i_mat[,(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]])

  results_theta_ij=data.frame(TF=TF_names,Theta_ij=tf_results_ij,Rank_ij=rank(-tf_results_ij))
  results_theta_ij=results_theta_ij[order(-results_theta_ij$Theta_ij),]
  row.names(results_theta_ij) <- NULL

  results_theta_i=data.frame(TF=TF_names,Theta_i=tf_results_i)
  results_theta_i=results_theta_i[!duplicated(results_theta_i$TF),]
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank_i=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  return_list<-list()
  return_list[["Theta_ij"]]<-results_theta_ij
  return_list[["Theta_i"]]<-results_theta_i
  return(return_list)
}



#' gelman-rubin demo
#'
#' @description Gelman-Rubin convergence diagnostic for Demo dataset, recommend for N > 2000
#' @param N number of rounds
#' @param dat pre-processed Demo data, including CTCF, ZBTB7A and KDM1A. For details please check the manual for the data.
#' @param nchains number of chains.
#'
#' @return A list object contains gelman-rubin diagnostic results.
#' @export
gelman_rubin_Demo <- function(N = 2000 , dat=c("CTCF","ZBTB7A","KDM1A"), nchains = parallel::detectCores()-1){
  if(nchains > parallel::detectCores()){
    stop("number of chains computed should be less or equal to number of cores.")
  }
  local_cl <- parallel::makeCluster(nchains)
  parallel::clusterExport(cl=local_cl, c("N","dat"),envir = environment())
  results <- parallel::clusterEvalQ(local_cl, {
    library(BayesIMTR)
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
