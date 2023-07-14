##Used for common covergence diagnostic, such as gelman-rubin convergence diagnostic.
##requires parallel computation, recommend to have cores > 10.

#' gelman-rubin convergence diagnostic
#'
#' @description gelman-rubin convergence diagnostic based on parallel computation, used chains can be varied but must less or equal to the usable cores.
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param N number of rounds for gibbs sampler, recommended for > 1000.
#' @param bin_width desired width of bin, should be in 100/500/1000.
#' @param nchains number of chains used in diagnostic.
#'
#' @return
#' @export
#'
#' @examples
gelman_rubin <- function(file, format=c("bed","narrowPeak","broadPeak","bigNarrowPeak"), N = 1000 ,
             bin_width = 1000, nchains = parallel::detectCores()){
  if(nchains > parallel::detectCores()){
    stop("number of chains computed should be less or equal to number of cores.")
  }
  results <- BayesIMTR_multi(file, format, N, bin_width, numCores = nchains)

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

#' mcmc_list_transformer
#'
#' @description function to transform parallel chains into a mcmc list.
#'
#' @param dat input results from parallel chains.
#' @param label label of TFs.
#' @param nchains number of chains.
#'
#' @return
#'
#' @examples
mcmc_list_transformer<-function(dat, label, nchains){
  dat_dims <- dim(dat[[1]][[label]])
  mcmc_trans_list <- list()
  if(dat_dims[2]==1){
    for(i in 1:nchains){
      mcmc_trans_list[[i]] <- coda::mcmc(dat[[i]][[label]])
    }
    return(mcmc_trans_list)
  }else{
    for(j in 1:dat_dims[1]){
      mcmc_trans_list[[j]] <- list()
    }

    for(i in 1:nchains){
      for(j in 1:dat_dims[1]){
        mcmc_trans_list[[j]][[i]] <- coda::mcmc(dat[[i]][[label]][j,])
      }
    }
    return(mcmc_trans_list)
  }
}
