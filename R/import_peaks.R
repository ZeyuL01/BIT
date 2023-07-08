##functions to transform peaks to windows.

import_peaks<-function(file,format=c("bed","narrowPeak","broadPeak"), window_length = 1000){
  window_inds<-c()

  if(format=="bed"){
    peak_dat<-rtracklayer::import(file,format="BED")
  }else if(format=="narrowPeak"){
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", peak = "integer")
    peak_dat <- import(file, format = "BED",
                            extraCols = extraCols_narrowPeak)

  }else if(format=="broadPeak"){
    extraCols_broadPeak <- c(singnalValue = "numeric", pValue = "numeric",
                             qValue = "numeric", peak = "integer")
    peak_dat <- rtracklayer::import(broadPeak_file, format = "BED",
                                    extraCols = extraCols_broadPeak)
  }

#Fixed window numbers with 3031030 windows in total.
  N <- 3031030
#Fixed windown numbers for each choromsome.
  chr_windows <- c(0,248956,242193,198295,190214,181538,170805,
                   159345,145138,138394,133797,135086,133275,
                   114364,107043,101991,90338,83257,80373,
                   58617,64444,46709,50818,156040)

  chr_windows_cs <- cumsum(chr_windows)

  chr_labels <- paste0("chr",c(1:22,"X"))
  for(i in 1:23){
    chr_lab <- chr_labels[i]
    chr_win_num <- chr_windows[i+1]

    if(format=="bed"){
      start <- rtracklayer::start(peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab])
      end <- rtracklayer::end(peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab])
      inds <- (end + start) %/% (2 * window_length) + 1
    }else{
      inds <- peak_dat$peak %/% window_length + 1
    }

    inds <- inds[inds<=chr_windows[i+1]]

    window_inds <- c(window_inds,(chr_windows_cs[i]+inds))
  }
  return(window_inds)
}



