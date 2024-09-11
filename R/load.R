##functions to transform peaks to windows.

#' Count the 'good' and 'total' informative cases between input peak set with reference database.
#'
#' @param input_vec A input vector contains index of peaks transformed by applying import_peaks.
#' @param bin_width width of bin, should be in 100/500/1000.
#'
#' @return a data frame has three columns, TR labels, number of 'good' windows, number of 'total' informative cases.
alignment_wrapper <- function(input_vec, bin_width, genome=c("hg38", "mm10")) {

  # Validate genome input
  genome <- match.arg(genome)

  # Load the meta table
  cat("Loading meta table...\n")
  meta_table_path <- system.file("meta_table.rds", package = "BIT")

  if (!file.exists(meta_table_path)) {
    stop("Meta table not found in package directory. Please ensure the package is installed correctly.")
  }

  meta_table <- readRDS(meta_table_path)

  # Get the appropriate ChIP-seq reference data
  file_table <- meta_table[[paste0("meta_", genome, "_", bin_width)]]

  if (is.null(file_table)) {
    stop("ChIP-seq files not found. Please download and load the ChIP-seq data first.\n",
         "You may follow the tutorial on: https://github.com/ZeyuL01/BIT")
  }

  # Initialize the results table
  chip_table <- data.frame(TR = file_table$TR, GOOD = numeric(nrow(file_table)), TOTAL = numeric(nrow(file_table)))

  # Progress bar setup
  cat("Starting alignment process...\n")
  pb <- txtProgressBar(min = 0, max = nrow(chip_table), style = 3, width = 50, char = "=")

  # Preallocate vectors for storing results
  good_vec <- numeric(nrow(chip_table))
  total_vec <- numeric(nrow(chip_table))

  # Loop over each ChIP-seq file and perform alignment
  for (i in seq_len(nrow(chip_table))) {
    ref_vec <- data.table::fread(file_table$File_Path[i])[[1]]

    # Perform alignment
    alignment_result <- Alignment(input_vec, ref_vec)

    # Store results in pre-allocated vectors
    good_vec[i] <- alignment_result$Xi_GOOD
    total_vec[i] <- alignment_result$Ni_TOTAL

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  # Close the progress bar
  close(pb)

  # Update chip_table with results
  chip_table$GOOD <- good_vec
  chip_table$TOTAL <- total_vec

  return(chip_table)
}


#' Transform the input file to vector
#'
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param bin_width desired width of bin, should be in 100/500/1000.
#'
#' @return A numeric vector contains the index of peaks with pre-specified number of bins in each chromosome.
#' @export
import_input_regions <- function(file, format = NULL, bin_width = 1000, genome=c("hg38", "mm10")) {
  # Determine the format if not provided
  if (is.null(format)) {
    extensions <- tools::file_ext(file)
    format <- switch(extensions,
                     bed = "bed",
                     bb = "bigNarrowPeak",
                     narrowPeak = "narrowPeak",
                     broadPeak = "broadPeak",
                     csv = "csv",
                     stop("File type not in (bed, bigNarrowPeak, narrowPeak, broadPeak, csv), please check the file type.")
    )
  }

  # Import data based on the format
  peak_dat <- switch(format,
                     bed = rtracklayer::import(file, format = "BED"),
                     narrowPeak = rtracklayer::import(file, format = "BED", extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")),
                     broadPeak = rtracklayer::import(file, format = "BED", extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")),
                     bigNarrowPeak = rtracklayer::import(file),
                     csv = read.csv(file),
                     stop("Unsupported format")
  )

  # Define fixed window numbers
  if (genome == "hg38") {
    N <- 3031030 * (1000 / bin_width)
    chr_windows <- c(0, 248956, 242193, 198295, 190214, 181538, 170805, 159345, 145138, 138394, 133797, 135086, 133275, 114364, 107043, 101991, 90338, 83257, 80373, 58617, 64444, 46709, 50818, 156040)
    chr_labels <- paste0("chr", c(1:22, "X"))
  } else if (genome == "mm10") {
    N <- 2631636 * (1000 / bin_width)
    chr_windows <- c(0, 195372, 182014, 159940, 156359, 151735, 149587, 145342, 129302, 124496, 130595, 121983, 120029, 120322, 124803, 103944, 98108, 94888, 90603, 61332, 170882)
    chr_labels <- paste0("chr", c(1:19, "X"))
  } else {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }

  chr_windows <- chr_windows * (1000 / bin_width)
  chr_windows_cs <- cumsum(chr_windows)

  bin_inds <- c()

  for (i in seq_along(chr_labels)) {
    chr_lab <- chr_labels[i]
    chr_win_num <- chr_windows[i + 1]

    inds <- switch(format,
                   bed = {
                     start <- rtracklayer::start(peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab])
                     end <- rtracklayer::end(peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab])
                     (end + start) %/% (2 * bin_width) + 1
                   },
                   narrowPeak = peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]$peak %/% bin_width + 1,
                   broadPeak = peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]$peak %/% bin_width + 1,
                   bigNarrowPeak = peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]$abs_summit %/% bin_width + 1,
                   csv = {
                     start <- peak_dat$Start[peak_dat$Chrom == chr_lab]
                     end <- peak_dat$End[peak_dat$Chrom == chr_lab]
                     (end + start) %/% (2 * bin_width) + 1
                   }
    )

    inds <- inds[inds <= chr_win_num]
    inds <- inds[!duplicated(inds)]
    bin_inds <- c(bin_inds, (chr_windows_cs[i] + inds))
  }

  return(sort(bin_inds))
}



##functions to load chip-seq datasets
##just need to run once.

#' load the pre-compiled chip-seq data.
#' @description load the pre-compiled chip-seq data. Please follow the tutorial on: https://github.com/ZeyuL01/BIT.
#' @param data_path path to the ChIP-seq data folder, can be absolute or relative path.
#' @param bin_width width of bin, which should be in 100/500/1000 and map with your ChIP-seq data.
#'
#' @export
load_chip_data <- function(data_path, bin_width, genome=c("hg38","mm10")){
  data_path = R.utils::getAbsolutePath(data_path)
  #Check the parameters
  if (!genome %in% c("hg38", "mm10")) {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }

  if(!bin_width %in% c(100,500,1000)){
    stop("bin width should be 100/500/1000!")
  }

  if(dir.exists(data_path)){
    if(!file.exists(paste0(system.file(package = "BIT"),"/meta_table.rds"))){

      data_list<-list()
      data_list[["path"]]=data_path

      ChIP_seq_files<-list.files(data_path)
      TR_labels<-sapply(strsplit(ChIP_seq_files,"_",fixed=TRUE),function(x){return(x[[1]])})
      meta_table<-data.frame(matrix(ncol=2,nrow=length(ChIP_seq_files)))
      colnames(meta_table)<-c("TR","File_Path")

      meta_table$TR <- TR_labels
      meta_table$File_Path <- paste0(data_path,"/",ChIP_seq_files)

      data_list[[paste0("meta_",genome,"_",bin_width)]] = meta_table

      saveRDS(data_list,paste0(system.file(package = "BIT"),"/meta_table.rds"))

    }else{
      data_list <- readRDS(paste0(system.file(package = "BIT"),"/meta_table.rds"))
      if(!is.null(data_list[[paste0("meta_",genome,"_",bin_width)]])){
        warning("Overwriting previous loaded meta-table for bin width of ", bin_width)
      }
      ChIP_seq_files<-list.files(data_path)
      TR_labels<-sapply(strsplit(ChIP_seq_files,"_",fixed=TRUE),function(x){return(x[[1]])})
      meta_table<-data.frame(matrix(ncol=2,nrow=length(ChIP_seq_files)))
      colnames(meta_table)<-c("TR","File_Path")

      meta_table$TR <- TR_labels
      meta_table$File_Path <- paste0(data_path,"/",ChIP_seq_files)

      data_list[[paste0("meta_",genome,"_",bin_width)]] = meta_table

      saveRDS(data_list,paste0(system.file(package = "BIT"),"/meta_table.rds"))
    }
  }else{

    stop("ChIP-seq data directory does not exist.")

  }
  print("ChIP-seq data successfully loaded, please run BIT with input to check!")
  return()
}

