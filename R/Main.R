##Main function

#' BIT
#' @description Main interface to run BIT method, please set the input file path, input file format, number of iterations and bin width.
#' @param file file path to the user-input.
#' @param show Whether to show the results table, default: TRUE.
#' @param plot_bar Whether to plot the top 10 TRs BIT score by a horizontal barplot, default: TRUE.
#' @param output_path absoluate or relative directory to store the Gibbs sampler data.
#' @param format if specify as NULL, BIT will automatically read and judge the file type based on extension.
#' @param N number of iterations for gibbs sampler, recommended for >= 5000, default: 5000.
#' @param bin_width desired width of bin, default: 1000.
#' @param genome Genome of the reference TR ChIP-seq data, must be "hg38" or "mm10".
#' @param seed Set random seed for reproducibility.
#'
#' @return NULL
#' @export
BIT <- function(file,
                output_path,
                show=TRUE,
                plot_bar=TRUE,
                format=NULL,
                N=5000,
                bin_width=1000,
                burnin=NULL,
                genome=c("hg38", "mm10"),
                seed=42) {

  # Define log file path
  output_path <- R.utils::getAbsolutePath(output_path)
  log_file <- file.path(output_path, paste0(tools::file_path_sans_ext(basename(file)), "_log.txt"))

  status <- tryCatch({

  # Start logging
  sink(log_file, append=TRUE)
  cat("BIT process started.\n")

  # Validate inputs
  if (!file.exists(file)) stop("Input file does not exist.")
  if (!dir.exists(output_path)) stop("Output path does not exist.")
  if (!genome %in% c("hg38", "mm10")) stop("Invalid genome specified.")

  # Standardize genome input
  genome <- match.arg(genome)

  cat("Loading and mapping peaks to bins...\n")

  # Resolve output path to an absolute path
  output_path <- R.utils::getAbsolutePath(output_path)

  # Import and map peaks
  input_peak_inds <- import_input_regions(file = file, format = format, bin_width = bin_width, genome = genome)

  cat("Done loading.\n")
  cat(paste0("Comparing the input regions with the pre-compiled reference ChIP-seq data, using a bin width of ",
             bin_width, " bps...\n"))

  # Perform alignment
  alignment_results <- alignment_wrapper(input_peak_inds, bin_width = bin_width, genome = genome)

  cat("Alignment complete.\n")

  # Extract relevant information for the Gibbs sampler
  xct <- alignment_results$GOOD
  nct <- alignment_results$TOTAL
  tr_labels <- as.numeric(factor(alignment_results$TR))

  # Set the seed for reproducibility
  set.seed(seed)  # Set the seed before running the Gibbs sampler

  # Run Gibbs sampling
  gibbs_sampler_results <- Main_Sampling(N, xct, nct, tr_labels, log_file)
  gibbs_sampler_results[["TR_names"]] <- alignment_results$TR

  # Save results to an RDS file
  file_name <- file.path(output_path, paste0(tools::file_path_sans_ext(basename(file)), ".rds"))
  saveRDS(gibbs_sampler_results, file_name)

  cat(paste0("Output data saved as ", file_name, "\n"))

  # Optionally display tables
  if (show) {
    display_tables(file_path = file_name, output_path = output_path, burnin = burnin)
  }

  # Optionally plot bar chart
  if (plot_bar) {
    rank_plot(file_path = file_name, output_path = output_path, burnin = burnin)
  }

  cat("BIT process completed.\n")
  cat("status: 0")

  sink()

  }, error = function(e) {
    # Error handler
    cat("An error occurred: ", conditionMessage(e), "\n")
    cat("BIT process failed.\n")
    cat("status: -1")

    # Stop logging and write the error message to the log file
    sink()

    # Return status code -1 for failure
    return(-1)
  })

  return(status)
}


#' BIT_compare
#' @description compare BIT identifid TRs for two user input epigenomic region sets.
#' @param file1 file path to the user-input file 1.
#' @param file2 file path to the user-input file 2.
#' @param show Whether to show the results table, default: TRUE
#' @param plot_scatter Whether to plot the BIT scores of two region sets in one scatter plot, default: TRUE
#' @param output_path absoluate or relative directory to store the Gibbs sampler data.
#' @param format if specify as NULL, BIT will automatically read and judge the file type based on extension, default: NULL.
#' @param N number of iterations for gibbs sampler, recommended for >= 5000, default: 5000.
#' @param bin_width desired width of bin, default: 1000.
#' @param genome Genome of the reference TR ChIP-seq data, must be "hg38" or "mm10"
#' @param seed Set random seed for reproducibility.
#'
#' @return NULL
#' @export
BIT_compare <- function(file1,
                        file2,
                        output_path,
                        show=TRUE,
                        plot_scatter=TRUE,
                        format=c(NULL, NULL),
                        N=5000,
                        bin_width=1000,
                        burnin=NULL,
                        genome=c("hg38", "mm10"),
                        seed=42) {

  # Validate inputs
  if (!file.exists(file1)) stop("Input file1 does not exist.")
  if (!file.exists(file2)) stop("Input file2 does not exist.")
  if (!dir.exists(output_path)) stop("Output path does not exist.")
  if (!genome %in% c("hg38", "mm10")) stop("Invalid genome specified.")

  # Standardize genome input
  genome <- match.arg(genome)

  cat("Loading and mapping peaks to bins for both files...\n")

  # Resolve output path to an absolute path
  output_path <- R.utils::getAbsolutePath(output_path)

  # Import and map peaks for both files
  input_peak_inds_file1 <- import_input_regions(file = file1, format = format[1], bin_width = bin_width, genome = genome)
  input_peak_inds_file2 <- import_input_regions(file = file2, format = format[2], bin_width = bin_width, genome = genome)

  cat("Done loading.\n")
  cat(paste0("Comparing input regions with pre-compiled reference ChIP-seq data, bin width: ", bin_width, " bps...\n"))

  # Perform alignment for both files
  alignment_results_file1 <- alignment_wrapper(input_peak_inds_file1, bin_width = bin_width, genome = genome)
  alignment_results_file2 <- alignment_wrapper(input_peak_inds_file2, bin_width = bin_width, genome = genome)

  cat("Alignment complete for both files.\n")

  # Extract relevant information for the Gibbs sampler for both files
  xct_file1 <- alignment_results_file1$GOOD
  nct_file1 <- alignment_results_file1$TOTAL
  xct_file2 <- alignment_results_file2$GOOD
  nct_file2 <- alignment_results_file2$TOTAL

  tr_labels_file1 <- as.numeric(factor(alignment_results_file1$TR))
  tr_labels_file2 <- as.numeric(factor(alignment_results_file2$TR))

  # Gibbs Sampling for file 1
  cat(paste0("Starting BIT Gibbs sampler for file 1, iterations: ", N, "...\n"))

  gibbs_sampler_results_file1 <- Main_Sampling(N, xct_file1, nct_file1, tr_labels_file1)
  gibbs_sampler_results_file1[["TR_names"]] <- alignment_results_file1$TR

  # Save file 1 results
  file1_name <- file.path(output_path, paste0(tools::file_path_sans_ext(basename(file1)), ".rds"))
  saveRDS(gibbs_sampler_results_file1, file1_name)

  cat(paste0("File 1 results saved as ", file1_name, "\n"))

  # Gibbs Sampling for file 2
  cat(paste0("Starting BIT Gibbs sampler for file 2, iterations: ", N, "...\n"))

  gibbs_sampler_results_file2 <- Main_Sampling(N, xct_file2, nct_file2, tr_labels_file2)
  gibbs_sampler_results_file2[["TR_names"]] <- alignment_results_file2$TR

  # Save file 2 results
  file2_name <- file.path(output_path, paste0(tools::file_path_sans_ext(basename(file2)), ".rds"))
  saveRDS(gibbs_sampler_results_file2, file2_name)

  cat(paste0("File 2 results saved as ", file2_name, "\n"))

  # Optionally display tables
  if (show) {
    display_tables(file_path = file1_name, output_path = output_path, burnin = burnin)
    display_tables(file_path = file2_name, output_path = output_path, burnin = burnin)
  }

  # Optionally plot scatter comparison
  if (plot_scatter) {
    compare_scatter_plot(file1_name, file2_name, output_path = output_path, burnin = burnin)
  }

  cat("BIT comparison process completed.\n")

}







