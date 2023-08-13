##Main function

#' BayesIMTR
#' @description Main interface to run BayesIMTR method, user need to give the input file path, input file format, number of rounds and bin width.
#' users can also change default parameters used in Gibbs sampler.
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param N number of rounds for gibbs sampler, recommended for > 2000.
#' @param bin_width desired width of bin, should be in 100/500/1000.
#'
#' @return A list object contains results of Gibbs sampler.
#' @export
BayesIMTR <- function(file, format=c("bed","narrowPeak","broadPeak","bigNarrowPeak"), N = 2000 ,bin_width = 1000, version = c(1,2), centerized = FALSE){
  print("Load and map peaks to bins...")
  input_peak_inds <- import_peaks(file = file, format = format, bin_width = bin_width)
  print("Done.")

  print(paste0("Align the loaded peaks with the pre-compiled reference ChIP-seq data, bin width used: ",bin_width," bps"))
  alignment_results <- alignment_wrapper(input_peak_inds, bin_width = bin_width)
  print("Done.")

  xct <- alignment_results$GOOD
  zct <- alignment_results$BAD
  nct <- alignment_results$TOTAL
  tf_labels <- as.numeric(factor(alignment_results$TF))

  print(paste0("Start BayesIMTR core, rounds: ",N))
  if(version == 1){
    gibbs_sampler_results <- Main_Sampling(N, xct, nct, tf_labels)
  }else{
    n2ct <- xct + zct
    gibbs_sampler_results <- Main_Sampling(N, xct, n2ct, tf_labels)
  }

  gibbs_sampler_results[["TF_names"]] <- alignment_results$TF
  print("Done.")
  if(centerized==TRUE){
    gibbs_sampler_results <- centerized(gibbs_sampler_results, tf_labels)
  }
  return(gibbs_sampler_results)
}


#' BayesIMTR parallel computation.
#' @description the multi-cores parallel version of BayesIMTR to acquire parallel MCMC chains for gelman-rubin diagnostic.
#'
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param N number of rounds for gibbs sampler, recommended for > 2000.
#' @param bin_width desired width of bin, should be in 100/500/1000.
#' @param numCores number of cores used in parallel computation.
#'
#' @return A list object of list objects that contains results of Gibbs sampler.
BayesIMTR_multi <- function(file, format=c("bed","narrowPeak","broadPeak","bigNarrowPeak"), N = 2000 ,bin_width = 1000, version = c(1,2), numCores){
    local_cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl=local_cl, c("file","format","N","bin_width","version"),envir = environment())
    multi_results <- parallel::clusterEvalQ(local_cl, {
      library(BayesIMTR)
      BayesIMTR(file = file,format = format,N = N, bin_width = bin_width, version = version)
    })
    parallel::stopCluster(local_cl)

    return(multi_results)
}






