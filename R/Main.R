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
  input_peak_inds <- import_peaks(file = file, format = format, bin_width = bin_width)

  alignment_results <- alignment_wrapper(input_peak_inds, bin_width = bin_width)

  xct <- alignment_results$GOOD
  nct <- alignment_results$TOTAL
  tf_labels <- alignment_results$TF

  gibbs_sampler_results <- Main_Sampling(N, xct, nct, tf_labels)

  return(gibbs_sampler_results)
}
