#' To show the ranking table by inspecting the results of Gibbs sampler.
#'
#' @param dat results from Main_sampling.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#'
#' @return a data.frame object contains TR names, BIT scores and 95 CI for each TR.
#' @export
display_tables<-function(dat, burnin=NULL){
  TR_names <- dat[["TR_names"]]
  theta_i_mat<-dat$theta_i

  if(is.null(burnin)){
    burnin=dim(theta_i_mat)[2]%/%2
  }

  tr_results_i<-rowMeans(theta_i_mat[,(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]])

  results_theta_i=data.frame(TR=TR_names,Theta_i=tr_results_i)
  results_theta_i=results_theta_i[!duplicated(results_theta_i$TR),]
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank_i=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  return_list<-list()
  return_list[["Theta_ij"]]<-results_theta_ij
  return_list[["Theta_i"]]<-results_theta_i

  return(return_list)
}

