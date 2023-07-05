#Show the demo

Demo = function(N,dat=c("CTCF","ZBTB7A","KDM1A")){
  if(dat=="CTCF"){
    dat = CTCF_Demo
  }else if(dat=="ZBTB7A"){
    dat = ZBTB7A_Demo
  }else if(dat=="KDM1A"){
    dat = KDM1A_Demo
  }

  tf_labs<-as.numeric(factor(dat$TF))
  xct<-dat$GOOD
  nct<-dat$TOTAL

  start_time <- Sys.time()

  results<-Main_Sampling(N,xct,nct,tf_labs)

  end_time <- Sys.time()
  time_used<-difftime(end_time,start_time)
  print(paste0("Time used: ", round(time_used,2)," ",units(time_used)," for ",N," rounds."))

  return(results)
}

Show_Results<-function(dat,burnin){
  TF_names=CTCF_Demo$TF
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


