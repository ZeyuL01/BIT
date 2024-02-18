##Main function

#' BIT
#' @description Main interface to run BIT method, please set the input file path, input file format, number of iterations and bin width.
#' @param file file path to the user-input.
#' @param output_path absoluate or relative directory to store the Gibbs sampler data
#' @param format if specify as NULL, BIT will automatically read and judge the file type based on extension.
#' @param N number of iterations for gibbs sampler, recommended for >= 5000, default: 5000.
#' @param bin_width desired width of bin, should be in 100/500/1000, default: 1000.
#' @param option option to filter peaks with candidate cis-regulatory elements from ENCODE,
#' default use "ALL", can be "PLS" for all promoters, "ELS" for all enhancers.
#'
#' @return NULL
#' @export
BIT <- function(file, output_path, show=TRUE, plot.bar=TRUE, format=NULL, N = 5000 ,bin_width = 1000, option="ALL",burnin=NULL){
  print("Load and map peaks to bins...")

  output_path = R.utils::getAbsolutePath(output_path)

  input_peak_inds <- import_input_regions(file = file, format = format, bin_width = bin_width)
  filtered_peak_inds <- filter_peaks(input_peak_inds, option = option, bin_width = bin_width)

  print("Done.")
  print(paste0("compare the input regions with the pre-compiled reference ChIP-seq data, bin width used: ",bin_width," bps"))

  alignment_results <- alignment_wrapper(filtered_peak_inds, bin_width = bin_width, option = option)

  print("Done.")

  xct <- alignment_results$GOOD
  nct <- alignment_results$TOTAL

  tr_labels <- as.numeric(factor(alignment_results$TR))

  print(paste0("Start BIT Gibbs sampler, iterations: ",N))

  gibbs_sampler_results <- Main_Sampling(N, xct, nct, tr_labels)

  gibbs_sampler_results[["TR_names"]] <- alignment_results$TR

  print("Done.")

  file_name<-paste0(output_path,"/",tools::file_path_sans_ext(basename(file)),".rds")
  saveRDS(gibbs_sampler_results,file_name)

  print(paste0("Output data saved as ",file_name))

  if(show==TRUE){
    display_tables(file_path = file_name, output_path = output_path, burnin = burnin)
  }

  if(plot.bar==TRUE){
    rank_plot(file_path = file_name, output_path = output_path, burnin = burnin)
  }

  return()
}

#' BIT_compare
#' @description compare BIT identifid TRs for two user input epigenomic region sets.
#' @param file1 file path to the user-input file 1.
#' @param file2 file path to the user-input file 2.
#' @param output_path absoluate or relative directory to store the Gibbs sampler data.
#' @param format if specify as NULL, BIT will automatically read and judge the file type based on extension, default: NULL.
#' @param N number of iterations for gibbs sampler, recommended for >= 5000, default: 5000.
#' @param bin_width desired width of bin, should be in 100/500/1000, default: 1000.
#' @param option option to filter peaks with candidate cis-regulatory elements from ENCODE,
#' default use "ALL", can be "PLS" for all promoters, "ELS" for all enhancers.
#'
#' @return NULL
#' @export
BIT_compare <- function(file1, file2, output_path, show=TRUE, plot.scatter=TRUE, format=c(NULL,NULL), N = 5000, bin_width = 1000, option="ALL", burnin=NULL){
  print("Load and map peaks to bins...")

  output_path = R.utils::getAbsolutePath(output_path)

  input_peak_inds_file1 <- import_input_regions(file = file1, format = format[1], bin_width = bin_width)
  input_peak_inds_file2 <- import_input_regions(file = file2, format = format[2], bin_width = bin_width)

  filtered_peak_inds_file1 <- filter_peaks(input_peak_inds_file1, option = option, bin_width = bin_width)
  filtered_peak_inds_file2 <- filter_peaks(input_peak_inds_file2, option = option, bin_width = bin_width)

  print("Done.")
  print(paste0("compare the input regions with the pre-compiled reference ChIP-seq data, bin width used: ",bin_width," bps"))

  alignment_results_file1 <- alignment_wrapper(filtered_peak_inds_file1, bin_width = bin_width, option = option)
  alignment_results_file2 <- alignment_wrapper(filtered_peak_inds_file2, bin_width = bin_width, option = option)

  print("Done.")

  xct_file1 <- alignment_results_file1$GOOD
  nct_file1 <- alignment_results_file1$TOTAL

  xct_file2 <- alignment_results_file2$GOOD
  nct_file2 <- alignment_results_file2$TOTAL

  tr_labels_file1 <- as.numeric(factor(alignment_results_file1$TR))
  tr_labels_file2 <- as.numeric(factor(alignment_results_file2$TR))

  print(paste0("Start BIT Gibbs sampler for file 1, iterations: ",N))

  gibbs_sampler_results_file1 <- Main_Sampling(N, xct_file1, nct_file1, tr_labels_file1)
  gibbs_sampler_results_file1[["TR_names"]] <- alignment_results_file1$TR
  file1_name<-paste0(output_path,"/",tools::file_path_sans_ext(basename(file1)),".rds")
  saveRDS(gibbs_sampler_results_file1,file1_name)

  print("Done.")
  print(paste0("file1 saved as ",file1_name))

  print(paste0("Start BIT Gibbs sampler for file 2, iterations: ",N))

  gibbs_sampler_results_file2 <- Main_Sampling(N, xct_file2, nct_file2, tr_labels_file2)
  gibbs_sampler_results_file2[["TR_names"]] <- alignment_results_file2$TR
  file2_name<-paste0(output_path,"/",tools::file_path_sans_ext(basename(file2)),".rds")
  saveRDS(gibbs_sampler_results_file2,file2_name)

  print("Done.")
  print(paste0("file2 saved as ",file2_name))

  if(show==TRUE){
    display_tables(file_path = file1_name, output_path = output_path, burnin = burnin)
    display_tables(file_path = file2_name, output_path = output_path, burnin = burnin)
  }

  if(plot.scatter==TRUE){
    compare_scatter_plot(file1_name,file2_name)
  }

  return()
}



BIT_parallel <- function(files_paths, output_path, show=TRUE, plot.bar=TRUE, format=NULL, N = 5000 ,bin_width = 1000, option="ALL",burnin=NULL){
  mclapply(files_paths, function(files_paths) BIT(files_paths, output_path=output_path, show=show, plot.bar=plot.bar, format=format,N=N,bin_width=bin_width,option=option, burnin=burnin),
           mc.cores = detectCores() - 1)
}

