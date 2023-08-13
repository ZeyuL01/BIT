##Compute centerized theta.

centerized<-function(dat, tf_labels){
  Theta_ij <- dat[["theta_ij"]]
  Theta_i <- dat[["theta_i"]]
  TF_nondup_id <- !duplicated(dat[["TF_names"]])

  sigmaS <- dat[["sigmaS"]]
  tau0S <- dat[["tau0S"]]

  centerized_theta_ij <- matrix(ncol=ncol(Theta_ij),nrow=nrow(Theta_ij))

  for(i in 1:length(tf_labels)){
    sigmaS_vec <- sigmaS[tf_labels[i]]
    centerized_theta_ij[i,] = Theta_ij[i] / sqrt(tau0S + sigmaS_vec)
  }

  dat[["centerized_theta_ij"]] <- centerized_theta_ij
  return(dat)
}


#' To show the ranking table by inspecting the results of Gibbs sampler.
#'
#' @param dat results from Main_sampling.
#' @param burnin number of samples used for burn-in.
#'
#' @return a list object contains two data.frame. one by ranking the theta_ij, the other one by ranking theta_i.
#' @export
Show_Results<-function(dat, burnin, centerized = FALSE){
  TF_names <- dat[["TF_names"]]
  theta_ij_mat<-dat$theta_ij
  theta_i_mat<-dat$theta_i

  tf_results_ij<-rowMeans(theta_ij_mat[,(dim(theta_ij_mat)[2]-burnin):dim(theta_ij_mat)[2]])
  tf_results_i<-rowMeans(theta_i_mat[,(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]])

  if(centerized==TRUE){
    theta_ij_centerized_mat<-dat$centerized_theta_ij
    tf_results_ij_centerized<-rowMeans(theta_ij_centerized_mat[,(dim(theta_ij_centerized_mat)[2]-burnin):dim(theta_ij_centerized_mat)[2]])
  }

  results_theta_ij=data.frame(TF=TF_names,Theta_ij=tf_results_ij,Rank_ij=rank(-tf_results_ij))
  results_theta_ij=results_theta_ij[order(-results_theta_ij$Theta_ij),]
  results_theta_ij=results_theta_ij[!duplicated(results_theta_ij$TF),]
  row.names(results_theta_ij) <- NULL

  results_theta_i=data.frame(TF=TF_names,Theta_i=tf_results_i)
  results_theta_i=results_theta_i[!duplicated(results_theta_i$TF),]
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank_i=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  if(centerized==TRUE){
    results_centerized_theta_ij = data.frame(TF=TF_names,Theta_ij_centerized=tf_results_ij_centerized,Rank_ij=rank(-tf_results_ij_centerized))
    results_centerized_theta_ij = results_centerized_theta_ij[order(-results_centerized_theta_ij$Theta_ij_centerized),]
    row.names(results_centerized_theta_ij) <- NULL
  }

  return_list<-list()
  return_list[["Theta_ij"]]<-results_theta_ij
  return_list[["Theta_i"]]<-results_theta_i
  if(centerized==TRUE){
    return_list[["Theta_ij_centerized"]]<-results_centerized_theta_ij
  }

  return(return_list)
}

