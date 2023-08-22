#' To show the ranking table by inspecting the results of Gibbs sampler.
#'
#' @param dat results from Main_sampling.
#' @param burnin number of samples used for burn-in.
#'
#' @return a list object contains two data.frame. one by ranking the theta_ij, the other one by ranking theta_i.
#' @export
Show_Results<-function(dat, burnin){
  TF_names <- dat[["TF_names"]]
  theta_ij_mat<-dat$theta_ij
  theta_i_mat<-dat$theta_i

  tf_results_ij<-rowMeans(theta_ij_mat[,(dim(theta_ij_mat)[2]-burnin):dim(theta_ij_mat)[2]])
  tf_results_i<-rowMeans(theta_i_mat[,(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]])

  results_theta_ij=data.frame(TF=TF_names,Theta_ij=tf_results_ij,Rank_ij=rank(-tf_results_ij))
  results_theta_ij=results_theta_ij[order(-results_theta_ij$Theta_ij),]
  results_theta_ij=results_theta_ij[!duplicated(results_theta_ij$TF),]
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

