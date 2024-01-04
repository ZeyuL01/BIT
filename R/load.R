##functions to transform peaks to windows.

#' Count the 'good' and 'total' informative cases between input peak set with reference database.
#'
#' @param input_vec A input vector contains index of peaks transformed by applying import_peaks.
#' @param bin_width width of bin, should be in 100/500/1000.
#'
#' @return a data frame has three columns, TR labels, number of 'good' windows, number of 'total' informative cases.
alignment_wrapper <- function(input_vec, bin_width, option){
  meta_table <- readRDS(paste0(system.file(package = "BIMTR"),"/meta_table.rds"))

  file_table <- meta_table[[paste0("meta_",bin_width)]]

  if(is.null(file_table)){
    stop("ChIP-seq files not found, please download and load the ChIP-seq data first.
         You may follow the tutorial on: https://github.com/ZeyuL01/BIMTR")
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

    filtered_ref_vec <- filter_peaks(ref_vec,option=option,bin_width=bin_width)

    alignment_result<-Alignment(input_vec,filtered_ref_vec)

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
import_peaks<-function(file,format=c("bed","narrowPeak","broadPeak","bigNarrowPeak","csv"), bin_width = 1000){
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
  }else if(format=="csv"){
    peak_dat <- read.csv(file)
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
    }else if(format=="csv"){
      start <- peak_dat$Start[which(peak_dat$Chrom==chr_lab)]
      end <- peak_dat$End[which(peak_dat$Chrom==chr_lab)]
      inds <- (end + start) %/% (2 * bin_width) + 1
    }

    inds <- inds[inds<=chr_windows[i+1]]
    inds <- inds[!duplicated(inds)]
    bin_inds <- c(bin_inds,(chr_windows_cs[i]+inds))
  }
  return(sort(bin_inds))
}


##functions to load chip-seq datasets
##just need to run once.

#' load the pre-compiled chip-seq data.
#' @description load the pre-compiled chip-seq data. Please follow the tutorial on: https://github.com/ZeyuL01/BIMTR.
#' @param data_path path to the ChIP-seq data folder.
#' @param bin_width width of bin, which should be in 100/500/1000 and map with your ChIP-seq data.
#'
#' @export
load_chip_data <- function(data_path, bin_width){
  data_path = R.utils::getAbsolutePath(data_path)

  if(dir.exists(data_path)){
    if(!file.exists(paste0(system.file(package = "BIMTR"),"/meta_table.rds"))){

      data_list<-list()
      data_list[["path"]]=data_path

      if(!bin_width %in% c(100,500,1000)){
        stop("bin width should be 100/500/1000!")
      }

      ChIP_seq_files<-list.files(data_path)
      TF_labels<-sapply(strsplit(ChIP_seq_files,"_",fixed=TRUE),function(x){return(x[[1]])})
      meta_table<-data.frame(matrix(ncol=2,nrow=length(ChIP_seq_files)))
      colnames(meta_table)<-c("TF","File_Path")

      meta_table$TF <- TF_labels
      meta_table$File_Path <- paste0(data_path,"/",ChIP_seq_files)

      data_list[[paste0("meta_",bin_width)]] = meta_table

      saveRDS(data_list,paste0(system.file(package = "BIMTR"),"/meta_table.rds"))

    }else{
      data_list <- readRDS(paste0(system.file(package = "BIMTR"),"/meta_table.rds"))
      if(!is.null(data_list[[paste0("meta_",bin_width)]])){
        warning("Overwriting previous loaded meta-table for bin width of ", bin_width)
      }
      ChIP_seq_files<-list.files(data_path)
      TF_labels<-sapply(strsplit(ChIP_seq_files,"_",fixed=TRUE),function(x){return(x[[1]])})
      meta_table<-data.frame(matrix(ncol=2,nrow=length(ChIP_seq_files)))
      colnames(meta_table)<-c("TF","File_Path")

      meta_table$TF <- TF_labels
      meta_table$File_Path <- paste0(data_path,"/",ChIP_seq_files)

      data_list[[paste0("meta_",bin_width)]] = meta_table

      saveRDS(data_list,paste0(system.file(package = "BIMTR"),"/meta_table.rds"))
    }
  }else{

    stop("ChIP-seq data directory does not exist.")

  }
  print("ChIP-seq data successfully loaded, please run BIMTR with input to check!")
  return()
}


#' filter peaks based on CREs, Promoter, Enhancer
#'
#' @param input_vec input vector that contains indices from import_peaks
#' @param option option, ALL: Non-filtering, CREs: All cis-regulatory elements, PLS: Promoters, ELS: Enhancers.
#' @param bin_width bin width
#'
#' @return a vector contains filtered bin indices.
filter_peaks <- function(input_vec,option=c("ALL","CREs","CLS","ELS"),bin_width=c(100,500,1000)){
  if(option == "ALL"){
    return_vec<-input_vec
  }else if(option == "CREs"){
    filter_vec <- CREs_list[[paste0(option,"_",bin_width)]]
    return_vec<-intersect(input_vec,filter_vec)
  }else if(option == "CLS"){
    filter_vec <- PLS_list[[paste0(option,"_",bin_width)]]
    return_vec<-intersect(input_vec,filter_vec)
  }else if(option == "ELS"){
    filter_vec <- ELS_list[[paste0(option,"_",bin_width)]]
    return_vec<-intersect(input_vec,filter_vec)
  }
  return(return_vec)
}
