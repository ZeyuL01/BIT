##functions to transform peaks to windows.

#' Count the 'good' and 'total' informative cases between input peak set with reference database.
#'
#' @param input_vec A input vector contains index of peaks transformed by applying import_peaks.
#' @param bin_width width of bin, should be in 100/500/1000.
#'
#' @return a data frame has three columns, TF labels, number of 'good' windows, number of 'total' informative cases.
alignment_wrapper <- function(input_vec, bin_width){
  meta_table <- readRDS(paste0(system.file(package = "BayesIMTR"),"/meta_table.rds"))

  file_table <- meta_table[[paste0("meta_",bin_width)]]

  if(is.null(file_table)){
    stop("ChIP-seq files not found, please download and load the ChIP-seq data first.
         You may follow the tutorial on: https://github.com/ZeyuL01/BayesIMTR")
  }

  chip_table<-data.frame(matrix(ncol=4,nrow=nrow(file_table)))
  colnames(chip_table) <- c("TF","GOOD","BAD","TOTAL")

  chip_table$TF <- file_table$TF

  good_vec <- c()
  bad_vec <- c()
  total_vec <- c()

  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = nrow(chip_table), # Maximum value of the progress bar
                       style = 1,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar

  for(i in 1:nrow(chip_table)){
    ref_vec <- data.table::fread(file_table$File_Path[i])[[1]]
    alignment_result<-Alignment(input_vec,ref_vec)

    good_vec <- c(good_vec, alignment_result$Xi_GOOD)
    bad_vec <- c(bad_vec, alignment_result$Xi_BAD)
    total_vec <- c(total_vec, alignment_result$Ni_TOTAL)
    setTxtProgressBar(pb, i)
  }

  chip_table$GOOD <- good_vec
  chip_table$BAD <- bad_vec
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
import_peaks<-function(file,format=c("bed","narrowPeak","broadPeak","bigNarrowPeak"), bin_width = 1000){
  bin_inds<-c()

  if(format=="bed"){
    peak_dat<-rtracklayer::import(file,format="BED")
  }else if(format=="narrowPeak"){
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", peak = "integer")
    peak_dat <- rtracklayer::import(file, format = "BED",
                            extraCols = extraCols_narrowPeak)

  }else if(format=="broadPeak"){
    extraCols_broadPeak <- c(singnalValue = "numeric", pValue = "numeric",
                             qValue = "numeric", peak = "integer")
    peak_dat <- rtracklayer::import(file, format = "BED",
                                    extraCols = extraCols_broadPeak)
  }else if(format=="bigNarrowPeak"){
    peak_dat <- rtracklayer::import(file)
  }

#Fixed window numbers with 3031030 windows in total.
  N <- 3031030 * (1000 / bin_width)
#Fixed window numbers for each chromosome.
  chr_windows <- c(0,248956,242193,198295,190214,181538,170805,
                   159345,145138,138394,133797,135086,133275,
                   114364,107043,101991,90338,83257,80373,
                   58617,64444,46709,50818,156040)

  chr_windows <- chr_windows * (1000 / bin_width)

  chr_windows_cs <- cumsum(chr_windows)

  chr_labels <- paste0("chr",c(1:22,"X"))
  for(i in 1:23){
    chr_lab <- chr_labels[i]
    chr_win_num <- chr_windows[i+1]

    if(format=="bed"){
      start <- rtracklayer::start(peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab])
      end <- rtracklayer::end(peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab])
      inds <- (end + start) %/% (2 * bin_width) + 1
    }else if(format=="narrowPeak" | format=="broadPeak"){
      inds <- peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab]$peak %/% bin_width + 1
    }else if(format=="bigNarrowPeak"){
      inds <- peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab]$abs_summit %/% bin_width + 1
    }

    inds <- inds[inds<=chr_windows[i+1]]
    inds <- inds[!duplicated(inds)]
    bin_inds <- c(bin_inds,(chr_windows_cs[i]+inds))
  }
  return(bin_inds)
}



