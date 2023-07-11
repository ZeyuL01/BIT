##Main function

#' Title
#'
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param N number of rounds for gibbs sampler, recommended for > 1000.
#' @param bin_width desired width of bin, should be in 100/500/1000.
#'
#' @return A list object contains results of Gibbs sampler.
#' @export
#'
#' @examples
#' input_peak <- "./CTCF.bed"
#' gibbs_results <- BayesIMTR(input_peak, format = "bed", N = 1000, bin_width = 1000)
BayesIMTR <- function(file, format=c("bed","narrowPeak","broadPeak","bigNarrowPeak"), N = 1000 ,bin_width = 1000){
  print("Load and map peaks to bins...")
  input_peak_inds <- import_peaks(file = file, format = format, bin_width = bin_width)
  print("Done.")

  print(paste0("Align the loaded peaks with the pre-compiled reference ChIP-seq data, bin width used: ",bin_width," bps"))
  alignment_results <- alignment_wrapper(input_peak_inds, bin_width = bin_width)
  print("Done.")

  xct <- alignment_results$GOOD
  nct <- alignment_results$TOTAL
  tf_labels <- as.numeric(factor(alignment_results$TF))

  print(paste0("Start BayesIMTR core, rounds: ",N))
  gibbs_sampler_results <- Main_Sampling(N, xct, nct, tf_labels)
  gibbs_sampler_results[["TF_names"]] <- alignment_results$TF
  print("Done.")
  return(gibbs_sampler_results)
}
