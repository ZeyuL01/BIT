#Show the demo

#' Demo to show BayesIMTR, which can be used without loading the ChIP-seq data.
#'
#' @param N number of rounds for gibbs sampler, recommended for > 1000.
#' @param dat pre-processed Demo data, including CTCF, ZBTB7A and KDM1A. For details please check the manual for the data.
#'
#' @return A list object contains results of Gibbs sampler.
#' @export
#'
#' @examples
#' CTCF_Demo <- Demo(1000, "CTCF")
#' CTCF_Results <- Show_Results(CTCF_Demo, burnin = 500)
#' head(CTCF_Results[["Theta_ij"]],10)
#' head(CTCF_Results[["Theta_i"]],10)
Demo = function(N,dat=c("CTCF","ZBTB7A","KDM1A")){
  if(dat=="CTCF"){
    Demo_dat = CTCF_Demo
  }else if(dat=="ZBTB7A"){
    Demo_dat = ZBTB7A_Demo
  }else if(dat=="KDM1A"){
    Demo_dat = KDM1A_Demo
  }

  tf_labs<-as.numeric(factor(Demo_dat$TF))
  xct<-Demo_dat$GOOD
  nct<-Demo_dat$TOTAL

  start_time <- Sys.time()

  results<-Main_Sampling(N,xct,nct,tf_labs)

  end_time <- Sys.time()
  time_used<-difftime(end_time,start_time)
  print(paste0("Time used: ", round(time_used,2)," ",units(time_used)," for ",N," rounds."))

  results[["TF_names"]]<-Demo_dat$TF
  return(results)
}

#' To show the ranking table by inspecting the results of Gibbs sampler.
#'
#' @param dat results from Main_sampling.
#' @param burnin number of samples used for burn-in.
#'
#' @return a list object contains two data.frame. one by ranking the theta_ij, the other one by ranking theta_i.
#' @export
#'
#' @examples
#' CTCF_Demo <- Demo(1000, "CTCF")
#' CTCF_Results <- Show_Results(CTCF_Demo, burnin = 500)
#' head(CTCF_Results[["Theta_ij"]],10)
#' head(CTCF_Results[["Theta_i"]],10)
Show_Results<-function(dat,burnin){
  TF_names <- dat[["TF_names"]]
  theta_ij_mat<-dat$theta_ij
  theta_i_mat<-dat$theta_i
  tf_results_ij<-rowMeans(theta_ij_mat[,(dim(theta_ij_mat)[2]-burnin):dim(theta_ij_mat)[2]])
  tf_results_i<-rowMeans(theta_i_mat[,(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]])

  results_theta_ij=data.frame(TF=TF_names,Theta_ij=tf_results_ij,Rank_ij=rank(-tf_results_ij))
  results_theta_ij=results_theta_ij[order(-results_theta_ij$Theta_ij),]
  row.names(results_theta_ij) <- NULL

  results_theta_i=data.frame(TF=TF_names,Theta_i=tf_results_i)
  results_theta_i=results_theta_i[!duplicated(results_theta_i$TF),]
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank_i=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  return_list<-list()
  return_list[["Theta_ij"]]<-results_theta_ij
  return_list[["Theta_i"]]<-results_theta_i
  return(return_list)
}


