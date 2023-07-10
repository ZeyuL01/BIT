##functions to load chip-seq datasets
##just need to run once.

#' load the pre-compiled chip-seq data.
#' @description load the pre-compiled chip-seq data. Please follow the tutorial on: https://github.com/ZeyuL01/BayesIMTR.
#' @param data_path path to the ChIP-seq data folder.
#' @param bin_width width of bin, which should be in 100/500/1000 and map with your ChIP-seq data.
#'
#' @export
#'
#' @examples
#' load_chip_data("data_path",bin_width = 1000)
load_chip_data <- function(data_path, bin_width){
  data_path = R.utils::getAbsolutePath(data_path)

  if(dir.exists(data_path)){
    if(!file.exists(paste0(system.file(package = "BayesIMTR"),"/meta_table.rds"))){

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
      meta_table$File_Path <- paste0(data_path,ChIP_seq_files)

      data_list[[paste0("meta_",bin_width)]] = meta_table

      saveRDS(data_list,paste0(system.file(package = "BayesIMTR"),"/meta_table.rds"))

    }else{
      data_list <- readRDS(paste0(system.file(package = "BayesIMTR"),"/meta_table.rds"))
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

      saveRDS(data_list,paste0(system.file(package = "BayesIMTR"),"/meta_table.rds"))
    }
  }else{

    stop("ChIP-seq data directory does not exist.")

  }
  print("ChIP-seq data successfully loaded, please run BayesIMTR with input to check!")
  return()
}
