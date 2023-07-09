##Main function

BayesIMTR <- function(file, format, N=1000 ,window_length = 1000){
  input_peak_inds <- import_peaks(file = file, format = format, window_length = window_length)

  alignment_results <- Alignment_wrapper(input_peak_inds)

  xct <- alignment_results$x
  nct <- alignment_results$n

  gibbs_sampler_results <- Main_Sampling(N, xct, nct, tf_labels)

  return(Show_Results(gibbs_sampler_results))
}
