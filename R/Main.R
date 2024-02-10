##Main function

#' BIT
#' @description Main interface to run BIT method, please set the input file path, input file format, number of iterations and bin width.
#' users can also change default parameters used in Gibbs sampler.
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param N number of iterations for gibbs sampler, recommended for >= 5000.
#' @param bin_width desired width of bin, should be in 100/500/1000.
#' @param option option to filter peaks with candidate cis-regulatory elements from ENCODE.
#'
#' @return A list object contains results of Gibbs sampler.
#' @export
BIT <- function(file, format=c("bed","narrowPeak","broadPeak","bigNarrowPeak","csv"), N = 5000 ,bin_width = c(100,500,1000), option=c("ALL","CREs","PLS","ELS")){
  print("Load and map peaks to bins...")
  input_peak_inds <- import_peaks(file = file, format = format, bin_width = bin_width)
  filtered_peak_inds <- filter_peaks(input_peak_inds,option = option, bin_width = bin_width)
  print("Done.")

  print(paste0("Align the loaded peaks with the pre-compiled reference ChIP-seq data, bin width used: ",bin_width," bps"))
  alignment_results <- alignment_wrapper(filtered_peak_inds, bin_width = bin_width, option = option)
  print("Done.")

  xct <- alignment_results$GOOD
  zct <- alignment_results$BAD
  nct <- alignment_results$TOTAL

  tf_labels <- as.numeric(factor(alignment_results$TF))

  print(paste0("Start BIMTR Gibbs sampler, iterations: ",N))

  gibbs_sampler_results <- Main_Sampling(N, xct, nct, tf_labels)

  gibbs_sampler_results[["TF_names"]] <- alignment_results$TF

  print("Done.")

  return(gibbs_sampler_results)

}






