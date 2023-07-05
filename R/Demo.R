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

Show_Results<-function(dat,TF_names=CTCF_Demo$TF,burnin){
  theta_ij_mat<-dat$theta_ij
  theta_i_mat<-dat$theta_i
  tf_results<-theta_ij_mat[,(dim(theta_ij_mat)[2]-burnin):dim(theta_ij_mat)[2]]
  tf_results_i<-theta_i_mat[,(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]]
  results_tab<-data.frame(TF=TF_names,res=rowMeans(tf_results),res_i=rowMeans(tf_results_i))
  results_tab$rank=rank(-results_tab$res)
  return(results_tab)
}


